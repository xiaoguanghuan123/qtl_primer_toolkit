"""
核心引物设计模块
处理模板准备、引物设计、初步质量评估等核心功能。
"""
import primer3
from pyfaidx import Fasta
from ..config.defaults import PRIMER3_DEFAULTS, TOOL_DEFAULTS
import logging

logger = logging.getLogger(__name__)

class PrimerDesigner:
    """引物设计器，封装primer3调用和模板处理逻辑"""
    
    def __init__(self, reference_fasta, config=None):
        """
        初始化设计器
        
        Args:
            reference_fasta: 参考基因组FASTA文件路径
            config: 配置字典，覆盖默认值
        """
        self.reference = Fasta(reference_fasta)
        self.config = PRIMER3_DEFAULTS.copy()
        if config:
            self.config.update(config)
        
        # 工具配置
        self.tool_config = TOOL_DEFAULTS.copy()
        
        logger.info(f"引物设计器初始化，参考基因组: {reference_fasta}")
    
    def prepare_template(self, chrom, pos, ref_len):
        """
        根据indel位置提取模板序列 (基于您脚本中的动态调整逻辑)
        
        Args:
            chrom: 染色体名称
            pos: indel起始位置 (1-based)
            ref_len: 参考等位基因长度
            
        Returns:
            dict: 包含模板序列、目标位置和原始坐标信息
                 如果失败则返回None
        """
        flanking = self.tool_config['FLANKING_SIZE']
        
        # 检查染色体是否存在
        if chrom not in self.reference:
            logger.error(f"染色体 {chrom} 不在参考基因组中")
            return None
        
        chr_length = len(self.reference[chrom])
        
        # 动态计算起始和结束位置 (基于您脚本中的逻辑)
        start = max(0, pos - flanking - 1)  # 转换为0-based
        end = min(chr_length, pos + ref_len + flanking - 1)
        
        # 检查模板长度是否足够
        template_length = end - start
        if template_length < self.tool_config['MIN_TEMPLATE_LENGTH']:
            logger.warning(f"模板序列过短: {template_length}bp (位置: {chrom}:{pos})")
            return None
        
        # 提取模板序列
        try:
            template_seq = str(self.reference[chrom][start:end]).upper()
        except Exception as e:
            logger.error(f"提取模板序列失败: {e}")
            return None
        
        # 计算目标位置在模板中的相对坐标 (0-based)
        target_start = (pos - 1) - start
        target_end = target_start + ref_len
        
        # 验证目标位置在模板范围内
        if target_end > len(template_seq):
            logger.error(f"目标位置超出模板范围: {chrom}:{pos}")
            return None
        
        logger.debug(f"模板准备完成: {chrom}:{pos}, 长度={template_length}bp, 目标位置=[{target_start},{target_end}]")
        
        return {
            'sequence': template_seq,
            'target': [target_start, ref_len],  # primer3需要的格式
            'original_coords': {
                'chrom': chrom,
                'pos': pos,
                'ref_len': ref_len,
                'template_start': start + 1,  # 转换回1-based
                'template_end': end,
                'flanking': flanking
            }
        }
    
    def calculate_product_size_range(self, ref_len):
        """
        动态计算产物大小范围 (基于您脚本中的逻辑)
        
        Args:
            ref_len: 参考等位基因长度
            
        Returns:
            list: [[min_size, max_size]] primer3需要的格式
        """
        min_size = ref_len + self.config['PRODUCT_SIZE_MIN_OFFSET']
        max_size = ref_len + self.config['PRODUCT_SIZE_MAX_OFFSET']
        
        # 确保最小值不超过最大值
        if min_size > max_size:
            min_size, max_size = max_size, min_size
            
        return [[min_size, max_size]]
    
    def design_primers_for_indel(self, chrom, pos, ref, alt, seq_id=None):
        """
        为单个indel设计引物 (封装您脚本中的核心循环逻辑)
        
        Args:
            chrom: 染色体
            pos: 位置
            ref: 参考等位基因
            alt: 替代等位基因
            seq_id: 序列ID
            
        Returns:
            list: 设计的引物对列表，每个元素是一个字典
        """
        if seq_id is None:
            seq_id = f"{chrom}_{pos}"
        
        ref_len = len(ref)
        alt_len = len(alt)
        
        logger.info(f"设计引物: {seq_id}, REF长度={ref_len}, ALT长度={alt_len}")
        
        # 1. 准备模板
        template_info = self.prepare_template(chrom, pos, ref_len)
        if not template_info:
            return [self._create_failed_result(seq_id, chrom, pos, ref, alt, "模板准备失败")]
        
        # 2. 准备primer3参数
        primer3_args = {
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': template_info['sequence'],
            'SEQUENCE_TARGET': template_info['target'],
        }
        
        # 3. 准备设计参数
        design_params = self.config.copy()
        
        # 动态设置产物大小范围
        design_params['PRIMER_PRODUCT_SIZE_RANGE'] = self.calculate_product_size_range(ref_len)
        
        # 移除自定义参数，保留primer3原生参数
        design_params.pop('PRODUCT_SIZE_MIN_OFFSET', None)
        design_params.pop('PRODUCT_SIZE_MAX_OFFSET', None)
        
        # 4. 调用primer3设计引物
        try:
            results = primer3.bindings.design_primers(primer3_args, design_params)
        except Exception as e:
            logger.error(f"调用primer3失败: {e}")
            return [self._create_failed_result(seq_id, chrom, pos, ref, alt, f"primer3错误: {str(e)}")]
        
        # 5. 解析结果
        num_pairs = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
        logger.info(f"{seq_id}: 成功设计 {num_pairs} 对引物")
        
        if num_pairs == 0:
            return [self._create_failed_result(seq_id, chrom, pos, ref, alt, "未设计出引物")]
        
        # 6. 提取和格式化引物对
        primer_pairs = []
        for i in range(num_pairs):
            pair_data = self._extract_primer_pair(results, i, seq_id, chrom, pos, ref, alt)
            
            # 计算额外长度 (您脚本中的Extra_Length)
            if pair_data['PRODUCT_SIZE'] and ref_len:
                pair_data['EXTRA_LENGTH'] = pair_data['PRODUCT_SIZE'] - ref_len
            
            # 质量评估 (基于您脚本中的逻辑)
            pair_data['QUALITY'] = self._assess_quality(pair_data)
            
            primer_pairs.append(pair_data)
        
        return primer_pairs
    
    def _extract_primer_pair(self, results, index, seq_id, chrom, pos, ref, alt):
        """从primer3结果中提取单个引物对信息"""
        prefix = f'PRIMER_PAIR_{index}'
        left_prefix = f'PRIMER_LEFT_{index}'
        right_prefix = f'PRIMER_RIGHT_{index}'
        
        # 提取基本信息
        pair_data = {
            'ID': seq_id,
            'CHR': chrom,
            'POS': pos,
            'REF': ref,
            'ALT': alt,
            'PAIR_INDEX': index,
            
            # 引物序列
            'LEFT_PRIMER': results.get(f'{left_prefix}_SEQUENCE', ''),
            'RIGHT_PRIMER': results.get(f'{right_prefix}_SEQUENCE', ''),
            
            # 位置信息 (0-based)
            'LEFT_POSITION': results.get(f'{left_prefix}', [0, 0])[0],
            'RIGHT_POSITION': results.get(f'{right_prefix}', [0, 0])[0],
            
            # 产物信息
            'PRODUCT_SIZE': results.get(f'{prefix}_PRODUCT_SIZE', 0),
            
            # 热力学参数
            'LEFT_TM': results.get(f'{left_prefix}_TM', 0),
            'RIGHT_TM': results.get(f'{right_prefix}_TM', 0),
            'LEFT_GC': results.get(f'{left_prefix}_GC_PERCENT', 0),
            'RIGHT_GC': results.get(f'{right_prefix}_GC_PERCENT', 0),
            
            # 二级结构 (primer3计算的)
            'SELF_ANY_TH': results.get(f'{prefix}_SELF_ANY_TH', 0),
            'SELF_END_TH': results.get(f'{prefix}_SELF_END_TH', 0),
            'HAIRPIN_TH': results.get(f'{prefix}_HAIRPIN_TH', 0),
        }
        
        # 计算引物自身二聚体和发卡结构 (基于您脚本中的逻辑)
        if pair_data['LEFT_PRIMER']:
            left_seq = pair_data['LEFT_PRIMER']
            homodimer = primer3.calc_homodimer(left_seq)
            hairpin = primer3.calc_hairpin(left_seq)
            
            pair_data['LEFT_SELF_ANY_DG'] = homodimer.dg if homodimer.structure_found else 0
            pair_data['LEFT_HAIRPIN_DG'] = hairpin.dg if hairpin.structure_found else 0
        
        if pair_data['RIGHT_PRIMER']:
            right_seq = pair_data['RIGHT_PRIMER']
            homodimer = primer3.calc_homodimer(right_seq)
            hairpin = primer3.calc_hairpin(right_seq)
            
            pair_data['RIGHT_SELF_ANY_DG'] = homodimer.dg if homodimer.structure_found else 0
            pair_data['RIGHT_HAIRPIN_DG'] = hairpin.dg if hairpin.structure_found else 0
        
        return pair_data
    
    def _assess_quality(self, pair_data):
        """评估引物对质量 (基于您脚本中的逻辑)"""
        quality_filters = self.tool_config['QUALITY_FILTERS']
        
        left_tm = pair_data.get('LEFT_TM', 0)
        right_tm = pair_data.get('RIGHT_TM', 0)
        
        # 检查Tm范围
        if not (quality_filters['MIN_TM'] <= left_tm <= quality_filters['MAX_TM']):
            return "POOR_TM"
        if not (quality_filters['MIN_TM'] <= right_tm <= quality_filters['MAX_TM']):
            return "POOR_TM"
        
        # 检查Tm差异
        tm_diff = abs(left_tm - right_tm)
        if tm_diff > quality_filters['MAX_TM_DIFF']:
            return "POOR_TM_DIFF"
        
        # 检查二级结构
        left_self_dg = pair_data.get('LEFT_SELF_ANY_DG', 0)
        right_self_dg = pair_data.get('RIGHT_SELF_ANY_DG', 0)
        
        if left_self_dg > quality_filters['MAX_SELF_ANY_DG']:
            return "POOR_SELF_DIMER"
        if right_self_dg > quality_filters['MAX_SELF_ANY_DG']:
            return "POOR_SELF_DIMER"
        
        return "GOOD"
    
    def _create_failed_result(self, seq_id, chrom, pos, ref, alt, reason):
        """创建失败结果记录 (保持与您脚本输出的一致性)"""
        return {
            'ID': seq_id,
            'CHR': chrom,
            'POS': pos,
            'REF': ref,
            'ALT': alt,
            'PAIR_INDEX': '',
            'LEFT_PRIMER': 'FAILED',
            'RIGHT_PRIMER': 'FAILED',
            'PRODUCT_SIZE': '',
            'EXTRA_LENGTH': '',
            'LEFT_TM': '',
            'RIGHT_TM': '',
            'LEFT_GC': '',
            'RIGHT_GC': '',
            'SELF_ANY_TH': '',
            'SELF_END_TH': '',
            'HAIRPIN_TH': '',
            'QUALITY': f'FAILED: {reason}'
        }

# 便捷函数
def design_primers_from_dataframe(df, reference_fasta, config=None):
    """
    批量从DataFrame设计引物
    
    Args:
        df: 包含CHR, POS, REF, ALT列的DataFrame
        reference_fasta: 参考基因组路径
        config: 配置字典
        
    Returns:
        list: 所有设计的引物对
    """
    designer = PrimerDesigner(reference_fasta, config)
    all_primers = []
    
    for _, row in df.iterrows():
        primers = designer.design_primers_for_indel(
            chrom=str(row['CHR']),
            pos=int(row['POS']),
            ref=str(row['REF']),
            alt=str(row['ALT']),
            seq_id=f"{row['CHR']}_{row['POS']}" if 'ID' not in row else row['ID']
        )
        all_primers.extend(primers)
    
    return all_primers
