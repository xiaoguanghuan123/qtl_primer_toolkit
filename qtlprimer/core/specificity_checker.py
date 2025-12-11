"""
批量BLAST验证模块
高效验证引物在参考基因组中的唯一性，避免非特异性结合。
"""
import subprocess
import tempfile
import os
import pandas as pd
import numpy as np
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.Blast.Applications import NcbiblastnCommandline
from ..config.defaults import BLAST_DEFAULTS

logger = logging.getLogger(__name__)

class SpecificityChecker:
    """批量BLAST验证引物特异性 - 优化版本"""
    
    def __init__(self, blast_db, config=None):
        """
        初始化BLAST检查器
        
        Args:
            blast_db: BLAST数据库路径（不含后缀）
            config: 配置字典，覆盖默认值
        """
        self.blast_db = blast_db
        self.config = BLAST_DEFAULTS.copy()
        if config:
            self.config.update(config)
        
        # 验证BLAST数据库是否存在
        self._validate_blast_db()
        
        logger.info(f"BLAST检查器初始化，数据库: {blast_db}")
    
    def _validate_blast_db(self):
        """验证BLAST数据库是否完整"""
        required_extensions = ['.nhr', '.nin', '.nsq']
        missing_files = []
        
        for ext in required_extensions:
            db_file = f"{self.blast_db}{ext}"
            if not os.path.exists(db_file):
                missing_files.append(db_file)
        
        if missing_files:
            error_msg = f"BLAST数据库不完整，缺少文件: {missing_files}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        logger.debug("BLAST数据库验证通过")
    
    def batch_blast_check(self, primers_data, batch_size=100):
        """
        批量BLAST检查 - 高性能版本
        
        Args:
            primers_data: 引物数据列表，每个元素为字典，包含：
                          - 'Primer_Name': 引物名称
                          - 'Sequence': 引物序列
                          - 其他引物属性（可选）
            batch_size: 每个BATCH的大小
        
        Returns:
            DataFrame: 包含BLAST验证结果的DataFrame
        """
        if not primers_data:
            logger.warning("没有提供引物数据")
            return pd.DataFrame()
        
        logger.info(f"开始批量BLAST验证，共 {len(primers_data)} 条引物")
        
        # 1. 将引物分组为批次
        batches = self._create_batches(primers_data, batch_size)
        
        # 2. 并行处理每个批次
        all_results = []
        with ProcessPoolExecutor(max_workers=self.config['NUM_THREADS']) as executor:
            future_to_batch = {
                executor.submit(self._process_batch, batch, batch_idx): batch_idx
                for batch_idx, batch in enumerate(batches)
            }
            
            for future in as_completed(future_to_batch):
                batch_idx = future_to_batch[future]
                try:
                    batch_results = future.result()
                    all_results.extend(batch_results)
                    logger.info(f"批次 {batch_idx} 处理完成: {len(batch_results)} 条结果")
                except Exception as e:
                    logger.error(f"批次 {batch_idx} 处理失败: {str(e)}")
        
        # 3. 合并所有结果
        results_df = pd.DataFrame(all_results)
        
        # 4. 计算特异性评分
        results_df = self._calculate_specificity_scores(results_df)
        
        logger.info(f"批量BLAST完成，共处理 {len(results_df)} 条引物")
        return results_df
    
    def _create_batches(self, primers_data, batch_size):
        """将引物数据分成批次"""
        batches = []
        for i in range(0, len(primers_data), batch_size):
            batch = primers_data[i:i + batch_size]
            batches.append(batch)
        
        logger.info(f"创建 {len(batches)} 个批次，每批最多 {batch_size} 条引物")
        return batches
    
    def _process_batch(self, batch, batch_idx):
        """处理单个批次的BLAST验证"""
        # 创建临时目录存放当前批次的文件
        with tempfile.TemporaryDirectory(prefix=f"blast_batch_{batch_idx}_") as temp_dir:
            # 1. 准备输入FASTA文件
            fasta_path = os.path.join(temp_dir, "primers.fasta")
            self._write_batch_fasta(batch, fasta_path)
            
            # 2. 运行BLAST
            blast_output_path = os.path.join(temp_dir, "blast_results.tsv")
            self._run_blastn(fasta_path, blast_output_path)
            
            # 3. 解析BLAST结果
            batch_results = self._parse_batch_results(batch, blast_output_path)
            
            return batch_results
    
    def _write_batch_fasta(self, batch, fasta_path):
        """将批次引物写入FASTA文件"""
        with open(fasta_path, 'w') as f:
            for primer in batch:
                f.write(f">{primer['Primer_Name']}\n")
                f.write(f"{primer['Sequence']}\n")
        
        logger.debug(f"写入FASTA文件: {fasta_path}, {len(batch)} 条序列")
    
    def _run_blastn(self, query_fasta, output_path):
        """运行BLASTN命令（使用Biopython接口）"""
        try:
            blastn_cline = NcbiblastnCommandline(
                task=self.config['TASK'],
                query=query_fasta,
                db=self.blast_db,
                out=output_path,
                outfmt="6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore",
                evalue=str(self.config['EVALUE_THRESHOLD']),
                num_threads=str(min(4, self.config['NUM_THREADS'])),  # 每个BLAST进程的线程数
                word_size=str(self.config['WORD_SIZE']),
                perc_identity=str(self.config['PERC_IDENTITY'])
            )
            
            logger.debug(f"BLAST命令: {blastn_cline}")
            stdout, stderr = blastn_cline()
            
            if stderr:
                logger.warning(f"BLAST标准错误: {stderr}")
                
        except Exception as e:
            logger.error(f"BLAST运行失败: {str(e)}")
            # 创建空输出文件
            with open(output_path, 'w') as f:
                pass
            raise
    
    def _parse_batch_results(self, batch, blast_output_path):
        """解析批次的BLAST结果"""
        results = []
        
        # 检查输出文件是否存在且不为空
        if not os.path.exists(blast_output_path) or os.path.getsize(blast_output_path) == 0:
            logger.warning(f"BLAST输出文件为空: {blast_output_path}")
            # 为批次中的每个引物创建无匹配结果
            for primer in batch:
                results.append(self._create_no_hit_result(primer))
            return results
        
        try:
            # 读取BLAST结果
            blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
            blast_df = pd.read_csv(blast_output_path, sep='\t', names=blast_cols)
            
            # 按引物名称分组处理
            primer_names = {primer['Primer_Name']: primer for primer in batch}
            
            for primer_name, primer_data in primer_names.items():
                primer_results = self._analyze_primer_hits(primer_name, primer_data, blast_df)
                results.append(primer_results)
                
        except Exception as e:
            logger.error(f"解析BLAST结果失败: {str(e)}")
            # 为批次中的每个引物创建错误结果
            for primer in batch:
                results.append(self._create_error_result(primer, str(e)))
        
        return results
    
    def _analyze_primer_hits(self, primer_name, primer_data, blast_df):
        """分析单个引物的BLAST匹配情况"""
        # 筛选该引物的所有匹配
        primer_hits = blast_df[blast_df['qseqid'] == primer_name].copy()
        primer_len = len(primer_data['Sequence'])
        
        # 基础统计
        total_hits = len(primer_hits)
        
        # 完全匹配统计
        exact_hits = primer_hits[
            (primer_hits['pident'] == 100.0) &
            (primer_hits['mismatch'] == 0) &
            (primer_hits['gapopen'] == 0) &
            (primer_hits['length'] == primer_len)
        ]
        exact_hit_count = len(exact_hits)
        
        # 高质量匹配统计（identity > 95%, mismatch <= 1）
        high_quality_hits = primer_hits[
            (primer_hits['pident'] >= 95.0) &
            (primer_hits['mismatch'] <= 1) &
            (primer_hits['gapopen'] == 0) &
            (primer_hits['length'] >= primer_len * 0.95)
        ]
        high_quality_count = len(high_quality_hits)
        
        # 获取最佳匹配
        best_hit = None
        if not primer_hits.empty:
            # 按bitscore排序，取最好的匹配
            best_idx = primer_hits['bitscore'].idxmax()
            best_hit = primer_hits.loc[best_idx]
        
        # 构建结果字典
        result = {
            'Primer_Name': primer_name,
            'Sequence': primer_data['Sequence'],
            'Primer_Length': primer_len,
            'Total_Hits': total_hits,
            'Exact_Hits': exact_hit_count,
            'High_Quality_Hits': high_quality_count,
            'Has_Exact_Match': exact_hit_count > 0,
            'Has_Multiple_Hits': total_hits > 1,
        }
        
        # 添加最佳匹配信息
        if best_hit is not None:
            result.update({
                'Best_Hit_Subject': best_hit['sseqid'],
                'Best_Hit_Identity': best_hit['pident'],
                'Best_Hit_Length': best_hit['length'],
                'Best_Hit_Mismatch': best_hit['mismatch'],
                'Best_Hit_Evalue': best_hit['evalue'],
                'Best_Hit_Bitscore': best_hit['bitscore'],
                'Best_Hit_Location': f"{best_hit['sstart']}-{best_hit['send']}",
            })
        else:
            result.update({
                'Best_Hit_Subject': None,
                'Best_Hit_Identity': None,
                'Best_Hit_Length': None,
                'Best_Hit_Mismatch': None,
                'Best_Hit_Evalue': None,
                'Best_Hit_Bitscore': None,
                'Best_Hit_Location': None,
            })
        
        # 添加原始引物数据中的其他字段
        for key, value in primer_data.items():
            if key not in result:
                result[key] = value
        
        return result
    
    def _calculate_specificity_scores(self, results_df):
        """计算引物特异性评分"""
        if results_df.empty:
            return results_df
        
        df = results_df.copy()
        
        # 特异性评分算法
        def calculate_score(row):
            score = 100.0  # 基础分
            
            # 1. 完全匹配惩罚（我们希望只有1个完全匹配）
            if row['Exact_Hits'] == 0:
                score -= 30  # 没有完全匹配，严重惩罚
            elif row['Exact_Hits'] == 1:
                score += 20  # 完美：只有1个完全匹配
            else:
                score -= 10 * min(row['Exact_Hits'], 10)  # 多个完全匹配，按数量惩罚
            
            # 2. 总匹配数惩罚
            if row['Total_Hits'] > 1:
                penalty = min(40, 5 * np.log1p(row['Total_Hits']))
                score -= penalty
            
            # 3. 高质量匹配惩罚
            if row['High_Quality_Hits'] > 1:
                score -= 5 * min(row['High_Quality_Hits'], 5)
            
            # 4. 基于最佳匹配的质量调整
            if pd.notna(row.get('Best_Hit_Identity')):
                identity = row['Best_Hit_Identity']
                if identity < 80:
                    score += 10  # 低相似度，可能不是真正匹配
                elif identity > 95 and row['Exact_Hits'] == 0:
                    score -= 5  # 高相似度但不是完全匹配
            
            # 确保分数在合理范围内
            return max(0, min(100, score))
        
        # 应用评分
        df['Specificity_Score'] = df.apply(calculate_score, axis=1)
        
        # 特异性等级
        def assign_grade(score):
            if score >= 90:
                return 'A+'
            elif score >= 80:
                return 'A'
            elif score >= 70:
                return 'B'
            elif score >= 60:
                return 'C'
            elif score >= 50:
                return 'D'
            else:
                return 'F'
        
        df['Specificity_Grade'] = df['Specificity_Score'].apply(assign_grade)
        
        # 排序（特异性评分降序）
        df = df.sort_values('Specificity_Score', ascending=False).reset_index(drop=True)
        
        return df
    
    def _create_no_hit_result(self, primer):
        """创建无匹配结果"""
        return {
            'Primer_Name': primer['Primer_Name'],
            'Sequence': primer['Sequence'],
            'Primer_Length': len(primer['Sequence']),
            'Total_Hits': 0,
            'Exact_Hits': 0,
            'High_Quality_Hits': 0,
            'Has_Exact_Match': False,
            'Has_Multiple_Hits': False,
            'Best_Hit_Subject': None,
            'Best_Hit_Identity': None,
            'Best_Hit_Length': None,
            'Best_Hit_Mismatch': None,
            'Best_Hit_Evalue': None,
            'Best_Hit_Bitscore': None,
            'Best_Hit_Location': None,
            'Specificity_Score': 100.0,  # 无匹配是最佳情况
            'Specificity_Grade': 'A+',
            **{k: v for k, v in primer.items() if k not in ['Primer_Name', 'Sequence']}
        }
    
    def _create_error_result(self, primer, error_msg):
        """创建错误结果"""
        result = self._create_no_hit_result(primer)
        result['BLAST_Error'] = error_msg
        result['Specificity_Score'] = 0
        result['Specificity_Grade'] = 'F'
        return result
    
    def check_specificity_for_dataframe(self, primers_df, primer_name_col='Primer_Name', 
                                        sequence_col='Sequence'):
        """
        为DataFrame中的引物运行特异性检查
        
        Args:
            primers_df: 包含引物信息的DataFrame
            primer_name_col: 引物名称列名
            sequence_col: 引物序列列名
            
        Returns:
            DataFrame: 包含特异性检查结果的DataFrame
        """
        # 准备引物数据
        primers_data = []
        for _, row in primers_df.iterrows():
            primer_data = {
                'Primer_Name': row[primer_name_col],
                'Sequence': row[sequence_col],
            }
            
            # 添加其他列作为额外信息
            for col in primers_df.columns:
                if col not in [primer_name_col, sequence_col]:
                    primer_data[col] = row[col]
            
            primers_data.append(primer_data)
        
        # 运行批量BLAST
        results_df = self.batch_blast_check(primers_data)
        
        return results_df


# 便捷函数
def check_primers_specificity(primers_df, blast_db, config=None, **kwargs):
    """
    便捷函数：检查DataFrame中引物的特异性
    
    Args:
        primers_df: 引物DataFrame
        blast_db: BLAST数据库路径
        config: BLAST配置
        **kwargs: 传递给SpecificityChecker的参数
        
    Returns:
        DataFrame: 包含特异性结果的新DataFrame
    """
    checker = SpecificityChecker(blast_db, config)
    results = checker.check_specificity_for_dataframe(primers_df, **kwargs)
    return results
