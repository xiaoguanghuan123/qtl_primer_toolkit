"""
工作流管理器
整合所有模块：变异解析 → 引物设计 → BLAST验证 → 可视化 → 报告生成
"""
import os
import logging
import json
import yaml
import pandas as pd
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, asdict, field
from typing import Dict, List, Optional, Tuple, Union, Any, Callable
from concurrent.futures import ThreadPoolExecutor, as_completed

from ..core.variant_parser import VariantParser, parse_variants
from ..core.primer_designer import PrimerDesigner, design_primers_from_dataframe
from ..core.specificity_checker import SpecificityChecker, check_primers_specificity
from ..visualization.gel_electrophoresis import GelElectrophoresisVisualizer
from ..visualization.recombination_plots import RecombinationPlotter
from ..config.defaults import DEFAULT_CONFIG, get_config, validate_config
from ..utils.file_utils import FileUtils
from ..utils.logging_utils import setup_logging

logger = logging.getLogger(__name__)


@dataclass
class WorkflowStep:
    """工作流步骤定义"""
    name: str
    description: str
    function: Callable
    enabled: bool = True
    depends_on: List[str] = field(default_factory=list)
    timeout: int = 3600  # 超时时间（秒）
    retry_count: int = 3
    
    def execute(self, *args, **kwargs):
        """执行步骤"""
        if not self.enabled:
            logger.info(f"跳过步骤: {self.name}")
            return None
        
        logger.info(f"开始步骤: {self.name}")
        logger.debug(f"步骤描述: {self.description}")
        
        for attempt in range(self.retry_count):
            try:
                result = self.function(*args, **kwargs)
                logger.info(f"步骤完成: {self.name}")
                return result
            except Exception as e:
                logger.warning(f"步骤 {self.name} 第 {attempt + 1} 次尝试失败: {str(e)}")
                if attempt == self.retry_count - 1:
                    logger.error(f"步骤 {self.name} 最终失败: {str(e)}")
                    raise
        return None


@dataclass
class WorkflowConfig:
    """工作流配置"""
    # 输入文件
    variants_file: str
    reference_fasta: str
    blast_db: str
    
    # 输出设置
    output_dir: str = "./qtl_results"
    sample_prefix: str = "QTL"
    create_subdirectories: bool = True
    
    # 处理参数
    max_workers: int = 8
    batch_size: int = 50
    save_intermediate: bool = True
    
    # 模块配置
    variant_parser_config: Dict = field(default_factory=dict)
    primer_designer_config: Dict = field(default_factory=dict)
    specificity_checker_config: Dict = field(default_factory=dict)
    visualization_config: Dict = field(default_factory=dict)
    
    # 步骤控制
    enabled_steps: List[str] = field(default_factory=lambda: [
        'parse_variants',
        'design_primers',
        'check_specificity',
        'create_visualizations',
        'generate_reports'
    ])
    
    def to_dict(self):
        """转换为字典"""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, config_dict: Dict) -> 'WorkflowConfig':
        """从字典创建配置"""
        return cls(**config_dict)
    
    def save(self, file_path: Union[str, Path]):
        """保存配置到文件"""
        file_path = Path(file_path)
        data = self.to_dict()
        
        if file_path.suffix in ['.yaml', '.yml']:
            with open(file_path, 'w') as f:
                yaml.dump(data, f, default_flow_style=False)
        elif file_path.suffix == '.json':
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
        else:
            raise ValueError(f"不支持的配置文件格式: {file_path.suffix}")
        
        logger.info(f"配置已保存: {file_path}")


@dataclass
class WorkflowResult:
    """工作流结果"""
    # 基本信息
    output_dir: Path
    start_time: datetime
    end_time: datetime
    success: bool
    
    # 结果数据
    config: WorkflowConfig
    indels_df: Optional[pd.DataFrame] = None
    primers_df: Optional[pd.DataFrame] = None
    blast_results_df: Optional[pd.DataFrame] = None
    visualization_results: Dict = field(default_factory=dict)
    report_files: Dict = field(default_factory=dict)
    
    # 统计信息
    statistics: Dict = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    
    # 性能信息
    step_times: Dict = field(default_factory=dict)
    
    @property
    def duration(self):
        """运行时长"""
        return self.end_time - self.start_time
    
    @property
    def indel_count(self):
        """Indel数量"""
        if self.indels_df is not None:
            return len(self.indels_df)
        return 0
    
    @property
    def primer_count(self):
        """引物数量"""
        if self.primers_df is not None:
            return len(self.primers_df)
        return 0
    
    @property
    def successful_primer_count(self):
        """成功设计的引物数量"""
        if self.primers_df is not None:
            return len(self.primers_df[self.primers_df['QUALITY'] != 'FAILED'])
        return 0
    
    def to_dict(self):
        """转换为字典"""
        result_dict = asdict(self)
        
        # 转换Path对象为字符串
        result_dict['output_dir'] = str(self.output_dir)
        
        # 转换DataFrame为字典（仅摘要）
        if self.indels_df is not None:
            result_dict['indels_summary'] = {
                'count': len(self.indels_df),
                'columns': list(self.indels_df.columns),
                'head': self.indels_df.head().to_dict('records')
            }
        
        if self.primers_df is not None:
            result_dict['primers_summary'] = {
                'count': len(self.primers_df),
                'successful': len(self.primers_df[self.primers_df['QUALITY'] != 'FAILED']),
                'columns': list(self.primers_df.columns),
                'top_primers': self.primers_df.head(10).to_dict('records')
            }
        
        return result_dict
    
    def save_summary(self, file_path: Union[str, Path]):
        """保存摘要到文件"""
        file_path = Path(file_path)
        summary = self.to_dict()
        
        if file_path.suffix in ['.yaml', '.yml']:
            with open(file_path, 'w') as f:
                yaml.dump(summary, f, default_flow_style=False)
        elif file_path.suffix == '.json':
            with open(file_path, 'w') as f:
                json.dump(summary, f, indent=2)
        else:
            with open(file_path, 'w') as f:
                f.write(self._generate_text_summary())
        
        logger.info(f"结果摘要已保存: {file_path}")
    
    def _generate_text_summary(self) -> str:
        """生成文本摘要"""
        lines = []
        lines.append("=" * 60)
        lines.append("QTL精细定位引物设计工作流 - 结果摘要")
        lines.append("=" * 60)
        
        lines.append(f"\n基本信息:")
        lines.append(f"  输出目录: {self.output_dir}")
        lines.append(f"  开始时间: {self.start_time}")
        lines.append(f"  结束时间: {self.end_time}")
        lines.append(f"  运行时长: {self.duration}")
        lines.append(f"  是否成功: {'是' if self.success else '否'}")
        
        lines.append(f"\n结果统计:")
        lines.append(f"  候选Indel数量: {self.indel_count}")
        lines.append(f"  设计引物对数: {self.primer_count}")
        lines.append(f"  成功设计引物: {self.successful_primer_count}")
        
        if self.statistics:
            lines.append(f"\n详细统计:")
            for key, value in self.statistics.items():
                lines.append(f"  {key}: {value}")
        
        if self.errors:
            lines.append(f"\n错误信息 ({len(self.errors)} 个):")
            for error in self.errors[:5]:  # 只显示前5个错误
                lines.append(f"  - {error}")
            if len(self.errors) > 5:
                lines.append(f"  ... 还有 {len(self.errors) - 5} 个错误")
        
        if self.warnings:
            lines.append(f"\n警告信息 ({len(self.warnings)} 个):")
            for warning in self.warnings[:5]:  # 只显示前5个警告
                lines.append(f"  - {warning}")
        
        lines.append(f"\n生成文件:")
        for file_type, file_path in self.report_files.items():
            lines.append(f"  {file_type}: {file_path}")
        
        lines.append(f"\n步骤耗时:")
        for step, duration in self.step_times.items():
            lines.append(f"  {step}: {duration:.2f}秒")
        
        lines.append("\n" + "=" * 60)
        return "\n".join(lines)


class QTLFineMappingWorkflow:
    """QTL精细定位完整工作流管理器"""
    
    def __init__(self, config: Optional[Union[Dict, WorkflowConfig]] = None):
        """
        初始化工作流管理器
        
        Args:
            config: 工作流配置
        """
        # 处理配置
        if config is None:
            self.config = WorkflowConfig(
                variants_file="",
                reference_fasta="",
                blast_db=""
            )
        elif isinstance(config, dict):
            self.config = WorkflowConfig.from_dict(config)
        elif isinstance(config, WorkflowConfig):
            self.config = config
        else:
            raise TypeError(f"不支持的配置类型: {type(config)}")
        
        # 工作流步骤
        self.steps = self._initialize_steps()
        self.results = {}
        self.current_step = 0
        self.total_steps = len(self.steps)
        
        # 性能监控
        self.start_time = None
        self.step_start_time = None
        
        logger.info("QTL精细定位工作流管理器初始化完成")
    
    def _initialize_steps(self) -> Dict[str, WorkflowStep]:
        """初始化工作流步骤"""
        steps = {}
        
        # 步骤1: 解析变异
        steps['parse_variants'] = WorkflowStep(
            name="parse_variants",
            description="解析变异文件，提取候选Indel",
            function=self._step_parse_variants,
            enabled='parse_variants' in self.config.enabled_steps
        )
        
        # 步骤2: 设计引物
        steps['design_primers'] = WorkflowStep(
            name="design_primers",
            description="为候选Indel设计PCR引物",
            function=self._step_design_primers,
            enabled='design_primers' in self.config.enabled_steps,
            depends_on=['parse_variants']
        )
        
        # 步骤3: BLAST特异性验证
        steps['check_specificity'] = WorkflowStep(
            name="check_specificity",
            description="使用BLAST验证引物特异性",
            function=self._step_check_specificity,
            enabled='check_specificity' in self.config.enabled_steps,
            depends_on=['design_primers']
        )
        
        # 步骤4: 创建可视化
        steps['create_visualizations'] = WorkflowStep(
            name="create_visualizations",
            description="生成电泳模拟图和重组图",
            function=self._step_create_visualizations,
            enabled='create_visualizations' in self.config.enabled_steps,
            depends_on=['check_specificity']
        )
        
        # 步骤5: 生成报告
        steps['generate_reports'] = WorkflowStep(
            name="generate_reports",
            description="生成综合报告",
            function=self._step_generate_reports,
            enabled='generate_reports' in self.config.enabled_steps,
            depends_on=['create_visualizations']
        )
        
        return steps
    
    def run(self, variants_file: Optional[str] = None,
            reference_fasta: Optional[str] = None,
            blast_db: Optional[str] = None,
            output_dir: Optional[str] = None) -> WorkflowResult:
        """
        运行完整工作流
        
        Args:
            variants_file: 变异文件路径（覆盖配置）
            reference_fasta: 参考基因组路径（覆盖配置）
            blast_db: BLAST数据库路径（覆盖配置）
            output_dir: 输出目录（覆盖配置）
            
        Returns:
            WorkflowResult: 工作流结果
        """
        # 更新配置
        if variants_file:
            self.config.variants_file = variants_file
        if reference_fasta:
            self.config.reference_fasta = reference_fasta
        if blast_db:
            self.config.blast_db = blast_db
        if output_dir:
            self.config.output_dir = output_dir
        
        # 验证必需参数
        if not all([self.config.variants_file, self.config.reference_fasta, self.config.blast_db]):
            raise ValueError("必需参数缺失: variants_file, reference_fasta, blast_db")
        
        # 创建输出目录
        output_dir_path = Path(self.config.output_dir)
        if self.config.create_subdirectories:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_dir_path = output_dir_path / f"{self.config.sample_prefix}_{timestamp}"
        
        output_dir_path.mkdir(parents=True, exist_ok=True)
        self.config.output_dir = str(output_dir_path)
        
        # 保存配置
        config_file = output_dir_path / "workflow_config.yaml"
        self.config.save(config_file)
        
        # 初始化结果
        self.start_time = datetime.now()
        result = WorkflowResult(
            output_dir=output_dir_path,
            start_time=self.start_time,
            end_time=self.start_time,
            success=False,
            config=self.config
        )
        
        # 运行工作流
        try:
            logger.info("=" * 60)
            logger.info("开始QTL精细定位工作流")
            logger.info("=" * 60)
            logger.info(f"变异文件: {self.config.variants_file}")
            logger.info(f"参考基因组: {self.config.reference_fasta}")
            logger.info(f"BLAST数据库: {self.config.blast_db}")
            logger.info(f"输出目录: {self.config.output_dir}")
            logger.info("=" * 60)
            
            # 执行步骤
            for step_name, step in self.steps.items():
                if not step.enabled:
                    logger.info(f"跳过禁用的步骤: {step_name}")
                    continue
                
                # 检查依赖
                for dep in step.depends_on:
                    if dep not in self.results:
                        raise RuntimeError(f"步骤 {step_name} 缺少依赖: {dep}")
                
                # 执行步骤
                step_start = datetime.now()
                try:
                    step_result = step.execute()
                    step_end = datetime.now()
                    step_duration = (step_end - step_start).total_seconds()
                    
                    self.results[step_name] = step_result
                    result.step_times[step_name] = step_duration
                    
                    logger.info(f"步骤 {step_name} 完成，耗时 {step_duration:.2f} 秒")
                    
                except Exception as e:
                    logger.error(f"步骤 {step_name} 执行失败: {e}")
                    result.errors.append(f"{step_name}: {str(e)}")
                    raise
            
            # 收集结果
            self._collect_results(result)
            result.success = True
            
            logger.info("=" * 60)
            logger.info("🎉 工作流成功完成!")
            logger.info("=" * 60)
            
        except Exception as e:
            logger.error(f"工作流执行失败: {e}")
            result.errors.append(f"工作流失败: {str(e)}")
            result.success = False
        
        finally:
            # 设置结束时间
            result.end_time = datetime.now()
            duration = (result.end_time - result.start_time).total_seconds()
            logger.info(f"总运行时间: {duration:.2f} 秒")
            
            # 保存结果摘要
            summary_file = output_dir_path / "workflow_summary.txt"
            result.save_summary(summary_file)
        
        return result
    
    def _step_parse_variants(self):
        """步骤1: 解析变异"""
        logger.info("步骤1: 解析变异文件")
        
        parser = VariantParser(config=self.config.variant_parser_config)
        df = parser.parse(self.config.variants_file)
        
        # 过滤indel
        indel_df = parser.filter_indels(df)
        
        if indel_df.empty:
            logger.warning("未找到符合条件的indel变异")
        
        # 保存结果
        if self.config.save_intermediate:
            output_file = Path(self.config.output_dir) / "01_indels_filtered.csv"
            indel_df.to_csv(output_file, index=False)
            logger.info(f"Indel结果已保存: {output_file}")
        
        return {
            'parser': parser,
            'df': df,
            'indel_df': indel_df,
            'output_file': str(output_file) if self.config.save_intermediate else None
        }
    
    def _step_design_primers(self):
        """步骤2: 设计引物"""
        logger.info("步骤2: 设计引物")
        
        # 获取上一步结果
        prev_result = self.results.get('parse_variants')
        if prev_result is None or prev_result['indel_df'] is None:
            raise ValueError("需要先完成变异解析步骤")
        
        indel_df = prev_result['indel_df']
        
        if indel_df.empty:
            logger.warning("没有Indel可供设计引物")
            return {
                'designer': None,
                'primers_df': pd.DataFrame(),
                'output_file': None
            }
        
        # 设计引物
        designer = PrimerDesigner(
            reference_fasta=self.config.reference_fasta,
            config=self.config.primer_designer_config
        )
        
        primers_df = design_primers_from_dataframe(
            indel_df, self.config.reference_fasta,
            config=self.config.primer_designer_config
        )
        
        # 保存结果
        output_file = None
        if self.config.save_intermediate and not primers_df.empty:
            output_dir = Path(self.config.output_dir)
            output_file = output_dir / "02_primers_designed.csv"
            primers_df.to_csv(output_file, index=False)
            
            # 保存FASTA文件
            fasta_file = output_dir / "02_primers_sequences.fasta"
            self._save_primers_fasta(primers_df, fasta_file)
            
            logger.info(f"引物设计结果已保存: {output_file}")
        
        return {
            'designer': designer,
            'primers_df': primers_df,
            'output_file': output_file
        }
    
    def _step_check_specificity(self):
        """步骤3: BLAST特异性验证"""
        logger.info("步骤3: BLAST特异性验证")
        
        # 获取上一步结果
        prev_result = self.results.get('design_primers')
        if prev_result is None or prev_result['primers_df'] is None:
            raise ValueError("需要先完成引物设计步骤")
        
        primers_df = prev_result['primers_df']
        
        if primers_df.empty:
            logger.warning("没有引物需要验证特异性")
            return {
                'checker': None,
                'blast_results': pd.DataFrame(),
                'output_file': None
            }
        
        # 准备引物数据
        primer_pairs = []
        for _, row in primers_df.iterrows():
            if row['LEFT_PRIMER'] not in ['FAILED', '']:
                primer_pairs.append({
                    'Primer_Name': f"{row['ID']}_P{row['PAIR_INDEX']}_F",
                    'Sequence': row['LEFT_PRIMER'],
                    'PRIMER_PAIR_ID': f"{row['ID']}_P{row['PAIR_INDEX']}",
                    'CHR': row['CHR'],
                    'POS': row['POS'],
                    'PRODUCT_SIZE': row['PRODUCT_SIZE']
                })
            
            if row['RIGHT_PRIMER'] not in ['FAILED', '']:
                primer_pairs.append({
                    'Primer_Name': f"{row['ID']}_P{row['PAIR_INDEX']}_R",
                    'Sequence': row['RIGHT_PRIMER'],
                    'PRIMER_PAIR_ID': f"{row['ID']}_P{row['PAIR_INDEX']}",
                    'CHR': row['CHR'],
                    'POS': row['POS'],
                    'PRODUCT_SIZE': row['PRODUCT_SIZE']
                })
        
        primers_for_blast = pd.DataFrame(primer_pairs)
        
        # BLAST验证
        checker = SpecificityChecker(
            blast_db=self.config.blast_db,
            config=self.config.specificity_checker_config
        )
        
        blast_results = checker.batch_blast_check(
            primer_pairs, 
            batch_size=self.config.batch_size
        )
        
        # 保存结果
        output_file = None
        if self.config.save_intermediate and not blast_results.empty:
            output_dir = Path(self.config.output_dir)
            output_file = output_dir / "03_primers_with_specificity.csv"
            blast_results.to_csv(output_file, index=False)
            logger.info(f"特异性验证结果已保存: {output_file}")
        
        return {
            'checker': checker,
            'primers_for_blast': primers_for_blast,
            'blast_results': blast_results,
            'output_file': output_file
        }
    
    def _step_create_visualizations(self):
        """步骤4: 创建可视化"""
        logger.info("步骤4: 创建可视化")
        
        # 获取上一步结果
        prev_result = self.results.get('check_specificity')
        if prev_result is None:
            raise ValueError("需要先完成特异性验证步骤")
        
        blast_results = prev_result.get('blast_results', pd.DataFrame())
        
        if blast_results.empty:
            logger.warning("没有BLAST结果可用于可视化")
            return {'visualization_dir': None, 'files': {}}
        
        # 创建可视化目录
        output_dir = Path(self.config.output_dir) / "04_visualizations"
        output_dir.mkdir(exist_ok=True)
        
        # 凝胶电泳可视化
        gel_visualizer = GelElectrophoresisVisualizer(
            config=self.config.visualization_config
        )
        
        # 选择前N个最佳引物进行可视化
        top_n = 3
        if 'Specificity_Score' in blast_results.columns:
            top_primers = blast_results.nlargest(top_n, 'Specificity_Score')
        else:
            top_primers = blast_results.head(top_n)
        
        visualization_files = {}
        
        for idx, row in top_primers.iterrows():
            primer_id = row.get('Primer_Name', f'primer_{idx}').replace('_F', '').replace('_R', '')
            
            # 生成三种场景的电泳图
            for scenario, (geno_a, geno_b) in gel_visualizer.genotype_examples.items():
                output_file = output_dir / f"{primer_id}_{scenario}.png"
                
                gel_visualizer.visualize_recombination_test(
                    marker_a_genotype=geno_a,
                    marker_b_genotype=geno_b,
                    marker_a_product_size=350,  # 默认大小
                    marker_b_product_size=550,
                    sample_name=primer_id,
                    output_path=str(output_file)
                )
                
                visualization_files[f"{primer_id}_{scenario}"] = str(output_file)
        
        # 生成教学图
        tutorial_file = output_dir / "gel_electrophoresis_tutorial.png"
        gel_visualizer.create_tutorial_figure(str(tutorial_file))
        visualization_files['tutorial'] = str(tutorial_file)
        
        # 重组图（如果有位置信息）
        plotter = RecombinationPlotter(config=self.config.visualization_config)
        
        # 示例：创建重组区间图
        if 'CHR' in blast_results.columns and 'POS' in blast_results.columns:
            # 获取第一个引物的染色体和位置信息
            first_primer = blast_results.iloc[0]
            chr_name = first_primer['CHR']
            
            # 收集所有引物的位置
            positions = []
            marker_names = []
            
            for _, row in blast_results.iterrows():
                if 'POS' in row and row['POS']:
                    positions.append(int(row['POS']))
                    marker_names.append(row['Primer_Name'])
            
            if len(positions) >= 2:
                plot_file = output_dir / "recombination_plot.png"
                plotter.plot_recombination_interval(
                    chromosome=chr_name,
                    positions=positions[:10],  # 限制数量
                    marker_names=marker_names[:10],
                    output_file=str(plot_file)
                )
                visualization_files['recombination_plot'] = str(plot_file)
        
        logger.info(f"创建了 {len(visualization_files)} 个可视化文件")
        
        return {
            'visualization_dir': str(output_dir),
            'files': visualization_files
        }
    
    def _step_generate_reports(self):
        """步骤5: 生成报告"""
        logger.info("步骤5: 生成报告")
        
        # 收集所有步骤的结果
        parse_result = self.results.get('parse_variants', {})
        design_result = self.results.get('design_primers', {})
        blast_result = self.results.get('check_specificity', {})
        viz_result = self.results.get('create_visualizations', {})
        
        # 生成报告
        from .reporter import ReportGenerator
        
        reporter = ReportGenerator()
        report_dir = Path(self.config.output_dir) / "05_reports"
        report_dir.mkdir(exist_ok=True)
        
        reports = reporter.generate_all_reports(
            indels_df=parse_result.get('indel_df'),
            primers_df=design_result.get('primers_df'),
            blast_results_df=blast_result.get('blast_results'),
            visualization_files=viz_result.get('files', {}),
            output_dir=str(report_dir),
            config=self.config
        )
        
        logger.info(f"生成了 {len(reports)} 个报告文件")
        
        return {
            'report_dir': str(report_dir),
            'reports': reports
        }
    
    def _save_primers_fasta(self, primers_df: pd.DataFrame, output_path: Path):
        """保存引物序列为FASTA格式"""
        with open(output_path, 'w') as f:
            for _, row in primers_df.iterrows():
                if row['LEFT_PRIMER'] not in ['FAILED', '']:
                    primer_id = f"{row['ID']}_P{row['PAIR_INDEX']}_F"
                    f.write(f">{primer_id}\n{row['LEFT_PRIMER']}\n")
                
                if row['RIGHT_PRIMER'] not in ['FAILED', '']:
                    primer_id = f"{row['ID']}_P{row['PAIR_INDEX']}_R"
                    f.write(f">{primer_id}\n{row['RIGHT_PRIMER']}\n")
    
    def _collect_results(self, result: WorkflowResult):
        """收集所有步骤的结果"""
        # 收集数据
        parse_result = self.results.get('parse_variants', {})
        design_result = self.results.get('design_primers', {})
        blast_result = self.results.get('check_specificity', {})
        viz_result = self.results.get('create_visualizations', {})
        report_result = self.results.get('generate_reports', {})
        
        result.indels_df = parse_result.get('indel_df')
        result.primers_df = design_result.get('primers_df')
        result.blast_results_df = blast_result.get('blast_results')
        result.visualization_results = viz_result.get('files', {})
        result.report_files = report_result.get('reports', {})
        
        # 收集统计信息
        result.statistics = {
            'indel_count': result.indel_count,
            'primer_count': result.primer_count,
            'successful_primer_count': result.successful_primer_count,
            'visualization_count': len(result.visualization_results),
            'report_count': len(result.report_files)
        }
        
        if result.blast_results_df is not None and not result.blast_results_df.empty:
            if 'Specificity_Score' in result.blast_results_df.columns:
                avg_score = result.blast_results_df['Specificity_Score'].mean()
                result.statistics['average_specificity_score'] = f"{avg_score:.2f}"
            
            if 'Specificity_Grade' in result.blast_results_df.columns:
                grade_counts = result.blast_results_df['Specificity_Grade'].value_counts().to_dict()
                result.statistics['specificity_grades'] = grade_counts


# 便捷函数
def run_qtl_fine_mapping_pipeline(
    variants_file: str,
    reference_fasta: str,
    blast_db: str,
    output_dir: str = "./qtl_results",
    config: Optional[Dict] = None
) -> WorkflowResult:
    """
    一键运行QTL精细定位完整流程
    
    Args:
        variants_file: 变异文件路径
        reference_fasta: 参考基因组路径
        blast_db: BLAST数据库路径
        output_dir: 输出目录
        config: 配置字典
        
    Returns:
        WorkflowResult: 工作流结果
    """
    # 创建配置
    workflow_config = WorkflowConfig(
        variants_file=variants_file,
        reference_fasta=reference_fasta,
        blast_db=blast_db,
        output_dir=output_dir
    )
    
    # 合并用户配置
    if config:
        for key, value in config.items():
            if hasattr(workflow_config, key):
                setattr(workflow_config, key, value)
    
    # 运行工作流
    workflow = QTLFineMappingWorkflow(workflow_config)
    result = workflow.run()
    
    return result


def run_workflow_from_config(config_file: str) -> WorkflowResult:
    """
    从配置文件运行工作流
    
    Args:
        config_file: 配置文件路径
        
    Returns:
        WorkflowResult: 工作流结果
    """
    config_file = Path(config_file)
    
    if not config_file.exists():
        raise FileNotFoundError(f"配置文件不存在: {config_file}")
    
    # 加载配置
    if config_file.suffix in ['.yaml', '.yml']:
        with open(config_file, 'r') as f:
            config_dict = yaml.safe_load(f)
    elif config_file.suffix == '.json':
        with open(config_file, 'r') as f:
            config_dict = json.load(f)
    else:
        raise ValueError(f"不支持的配置文件格式: {config_file.suffix}")
    
    # 运行工作流
    workflow = QTLFineMappingWorkflow(config_dict)
    result = workflow.run()
    
    return result
