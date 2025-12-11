"""
变异文件解析器
支持多种格式的变异文件解析，专门提取indel用于引物设计。
"""
import pandas as pd
import numpy as np
import logging
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class VariantRecord:
    """变异记录数据类"""
    chrom: str
    pos: int
    ref: str
    alt: str
    id: Optional[str] = None
    qual: Optional[float] = None
    filter: Optional[str] = None
    info: Optional[str] = None
    format: Optional[str] = None
    samples: Optional[Dict] = None
    
    @property
    def is_indel(self) -> bool:
        """判断是否为indel"""
        return len(self.ref) != len(self.alt)
    
    @property
    def indel_type(self) -> str:
        """返回indel类型"""
        if not self.is_indel:
            return "SNP"
        if len(self.ref) > len(self.alt):
            return "DELETION"
        elif len(self.ref) < len(self.alt):
            return "INSERTION"
        else:
            return "COMPLEX"
    
    @property
    def indel_length(self) -> int:
        """返回indel长度（插入为正，删除为负）"""
        return len(self.alt) - len(self.ref)
    
    @property
    def indel_size(self) -> int:
        """返回indel大小（绝对值）"""
        return abs(self.indel_length)

class VariantParser:
    """多格式变异文件解析器"""
    
    SUPPORTED_FORMATS = ['csv', 'tsv', 'vcf', 'vcf.gz', 'txt', 'bed']
    
    # 列名映射：将各种可能的列名映射到标准列名
    COLUMN_MAPPINGS = {
        'CHR': ['chr', 'chrom', 'chromosome', 'contig', 'scaffold', 'seqname', '#CHROM'],
        'POS': ['pos', 'position', 'start', 'coordinate'],
        'REF': ['ref', 'reference', 'ref_allele', 'ref_alleles'],
        'ALT': ['alt', 'alternate', 'alt_allele', 'alt_alleles', 'variant'],
        'ID': ['id', 'snp_id', 'variant_id', 'rsid', 'name'],
        'QUAL': ['qual', 'quality', 'phred'],
        'FILTER': ['filter', 'filters'],
        'INFO': ['info', 'information'],
    }
    
    def __init__(self, config: Optional[Dict] = None):
        """
        初始化解析器
        
        Args:
            config: 配置字典，可包含：
                - min_indel_length: 最小indel长度
                - max_indel_length: 最大indel长度
                - allowed_variant_types: 允许的变异类型列表
                - require_indel: 是否只返回indel
        """
        self.config = config or {}
        
        # 默认配置
        self.default_config = {
            'min_indel_length': 1,
            'max_indel_length': 100,
            'allowed_variant_types': ['INSERTION', 'DELETION'],
            'require_indel': True,
            'strict_mode': False,  # 严格模式，验证所有变异
        }
        
        # 合并配置
        self.default_config.update(self.config)
        
        logger.info("变异文件解析器初始化完成")
    
    def parse(self, file_path: Union[str, Path], 
              file_format: Optional[str] = None,
              **kwargs) -> pd.DataFrame:
        """
        解析变异文件
        
        Args:
            file_path: 文件路径
            file_format: 文件格式，如果为None则自动检测
            **kwargs: 传递给解析函数的参数
            
        Returns:
            DataFrame: 标准化的变异DataFrame
            
        Raises:
            FileNotFoundError: 文件不存在
            ValueError: 不支持的格式或解析错误
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"文件不存在: {file_path}")
        
        # 自动检测格式
        if file_format is None:
            file_format = self._detect_format(file_path)
        
        file_format = file_format.lower()
        
        logger.info(f"解析文件: {file_path}, 格式: {file_format}")
        
        # 根据格式调用相应的解析方法
        try:
            if file_format in ['csv', 'tsv', 'txt']:
                df = self._parse_delimited(file_path, file_format, **kwargs)
            elif file_format in ['vcf', 'vcf.gz']:
                df = self._parse_vcf(file_path, **kwargs)
            elif file_format == 'bed':
                df = self._parse_bed(file_path, **kwargs)
            else:
                raise ValueError(f"不支持的格式: {file_format}. 支持的格式: {self.SUPPORTED_FORMATS}")
            
            # 如果DataFrame为空，返回空DataFrame
            if df.empty:
                logger.warning("文件为空或没有数据行")
                return df
            
            # 标准化DataFrame
            df = self._standardize_dataframe(df)
            
            # 验证必需字段
            self._validate_dataframe(df)
            
            # 添加变异信息列
            df = self._add_variant_info(df)
            
            # 过滤indel
            if self.default_config['require_indel']:
                df = self.filter_indels(df)
            
            logger.info(f"解析完成: {len(df)} 个变异 (其中 {df['Is_Indel'].sum()} 个indel)")
            
            return df
            
        except Exception as e:
            logger.error(f"解析文件失败 {file_path}: {str(e)}")
            raise
    
    def _detect_format(self, file_path: Path) -> str:
        """自动检测文件格式"""
        suffix = file_path.suffix.lower()
        
        if suffix == '.csv':
            return 'csv'
        elif suffix == '.tsv':
            return 'tsv'
        elif suffix == '.txt':
            return 'txt'
        elif suffix == '.vcf':
            return 'vcf'
        elif suffix == '.gz':
            # 检查.gz文件的实际格式
            stem_suffix = file_path.stem.lower().split('.')[-1]
            if stem_suffix == 'vcf':
                return 'vcf.gz'
            elif stem_suffix in ['csv', 'tsv', 'txt']:
                return stem_suffix
        
        # 尝试通过文件内容检测
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                
                # 检查VCF特征
                if first_line.startswith('##fileformat=VCF'):
                    return 'vcf'
                elif first_line.startswith('#'):
                    # 可能是带标题行的CSV/TSV
                    return 'csv'  # 默认尝试CSV
        except:
            pass
        
        # 默认尝试CSV
        logger.warning(f"无法确定文件格式，尝试作为CSV解析: {file_path}")
        return 'csv'
    
    def _parse_delimited(self, file_path: Path, 
                        file_format: str, **kwargs) -> pd.DataFrame:
        """解析分隔符分隔的文件（CSV/TSV/TXT）"""
        # 确定分隔符
        delimiter = kwargs.get('delimiter', None)
        if delimiter is None:
            if file_format == 'csv':
                delimiter = ','
            elif file_format == 'tsv':
                delimiter = '\t'
        
        # 读取参数
        read_kwargs = {
            'sep': delimiter,
            'engine': 'python' if delimiter is None else 'python',
            'comment': kwargs.get('comment', None),
            'skiprows': kwargs.get('skiprows', 0),
            'header': kwargs.get('header', 'infer'),
            'dtype': str,  # 所有列先读为字符串
            'low_memory': kwargs.get('low_memory', False),
        }
        
        # 清理参数，移除None值
        read_kwargs = {k: v for k, v in read_kwargs.items() if v is not None}
        
        try:
            # 处理gzip压缩文件
            if str(file_path).endswith('.gz'):
                import gzip
                with gzip.open(file_path, 'rt') as f:
                    # 对于gzip文件，使用pandas的压缩支持
                    df = pd.read_csv(file_path, compression='gzip', **read_kwargs)
            else:
                df = pd.read_csv(file_path, **read_kwargs)
                
        except Exception as e:
            logger.error(f"解析分隔符文件失败: {e}")
            # 尝试自动检测分隔符
            try:
                df = pd.read_csv(file_path, sep=None, engine='python')
            except Exception as e2:
                logger.error(f"自动检测分隔符也失败: {e2}")
                raise
        
        # 清理列名
        df.columns = df.columns.str.strip()
        
        # 尝试自动识别列
        df = self._auto_detect_columns(df)
        
        return df
    
    def _parse_vcf(self, file_path: Path, **kwargs) -> pd.DataFrame:
        """解析VCF文件（支持gzip压缩）"""
        # 检查是否是gzip压缩
        is_gzipped = str(file_path).endswith('.gz')
        
        header_lines = []
        data_lines = []
        
        try:
            if is_gzipped:
                with gzip.open(file_path, 'rt') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('##'):
                            header_lines.append(line)
                        elif line.startswith('#'):
                            # 列标题行
                            header = line[1:].split('\t')
                        else:
                            data_lines.append(line)
            else:
                with open(file_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('##'):
                            header_lines.append(line)
                        elif line.startswith('#'):
                            # 列标题行
                            header = line[1:].split('\t')
                        else:
                            data_lines.append(line)
        except Exception as e:
            logger.error(f"读取VCF文件失败: {e}")
            raise
        
        # 如果没有数据行，返回空DataFrame
        if not data_lines:
            logger.warning("VCF文件没有数据行")
            return pd.DataFrame()
        
        # 解析数据行
        data = []
        for line in data_lines:
            fields = line.split('\t')
            if len(fields) < 8:  # VCF至少需要8个固定字段
                logger.warning(f"跳过格式不正确的行: {line[:50]}...")
                continue
            
            # 确保有足够的字段
            while len(fields) < len(header):
                fields.append('')
            
            data.append(fields[:len(header)])  # 只取与标题对应的字段
        
        # 创建DataFrame
        df = pd.DataFrame(data, columns=header)
        
        # 重命名标准列
        column_mapping = {}
        for col in df.columns:
            col_upper = col.upper()
            if col_upper == '#CHROM':
                column_mapping[col] = 'CHR'
            elif col_upper in ['CHROM', 'CHROMOSOME']:
                column_mapping[col] = 'CHR'
            elif col_upper == 'POS':
                column_mapping[col] = 'POS'
            elif col_upper == 'REF':
                column_mapping[col] = 'REF'
            elif col_upper == 'ALT':
                column_mapping[col] = 'ALT'
            elif col_upper == 'ID':
                column_mapping[col] = 'ID'
            elif col_upper == 'QUAL':
                column_mapping[col] = 'QUAL'
            elif col_upper == 'FILTER':
                column_mapping[col] = 'FILTER'
            elif col_upper == 'INFO':
                column_mapping[col] = 'INFO'
            elif col_upper == 'FORMAT':
                column_mapping[col] = 'FORMAT'
        
        if column_mapping:
            df = df.rename(columns=column_mapping)
            logger.debug(f"重命名VCF列: {column_mapping}")
        
        return df
    
    def _parse_bed(self, file_path: Path, **kwargs) -> pd.DataFrame:
        """解析BED文件"""
        try:
            # BED文件通常有3-12列，我们至少需要前3列
            df = pd.read_csv(file_path, sep='\t', header=None, 
                           names=['CHR', 'START', 'END'], dtype=str)
            
            # BED是0-based，半开区间，我们转换为1-based位置
            # 假设变异位置在起始位置
            df['POS'] = pd.to_numeric(df['START']) + 1  # 转换为1-based
            
            # BED文件通常没有REF/ALT，我们需要设置默认值
            # 这里可以扩展从其他列提取信息
            df['REF'] = 'N'
            df['ALT'] = 'N'
            
            return df
            
        except Exception as e:
            logger.error(f"解析BED文件失败: {e}")
            raise
    
    def _auto_detect_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """自动识别和重命名列"""
        column_mapping = {}
        
        # 查找匹配的列
        for standard_name, variants in self.COLUMN_MAPPINGS.items():
            for col in df.columns:
                col_lower = col.lower()
                if col_lower in variants or standard_name.lower() == col_lower:
                    column_mapping[col] = standard_name
                    break
        
        # 应用重命名
        if column_mapping:
            df = df.rename(columns=column_mapping)
            logger.debug(f"自动重命名列: {column_mapping}")
        
        return df
    
    def _standardize_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """标准化DataFrame"""
        if df.empty:
            # 返回包含标准列的空的DataFrame
            standard_cols = ['CHR', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO']
            return pd.DataFrame(columns=standard_cols)
        
        df = df.copy()
        
        # 确保必需列存在，如果不存在则创建
        required_cols = ['CHR', 'POS', 'REF', 'ALT']
        
        for col in required_cols:
            if col not in df.columns:
                logger.warning(f"DataFrame缺少列: {col}，创建空列")
                df[col] = ''
        
        # 清理数据
        for col in ['CHR', 'REF', 'ALT']:
            if col in df.columns:
                df[col] = df[col].astype(str).str.strip().str.upper()
        
        # 确保POS是整数
        if 'POS' in df.columns:
            # 尝试转换为数值，错误的值转换为NaN
            df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
            # 将NaN替换为0，然后转换为整数
            df['POS'] = df['POS'].fillna(0).astype(int)
        
        # 创建唯一ID（如果不存在）
        if 'ID' not in df.columns or df['ID'].isnull().all():
            df['ID'] = df.apply(
                lambda row: f"{row['CHR']}_{row['POS']}_{row['REF']}_{row['ALT']}",
                axis=1
            )
        
        # 确保其他列为字符串
        for col in ['ID', 'QUAL', 'FILTER', 'INFO']:
            if col in df.columns:
                df[col] = df[col].astype(str)
        
        return df
    
    def _validate_dataframe(self, df: pd.DataFrame):
        """验证DataFrame的完整性"""
        if df.empty:
            return
        
        required_cols = ['CHR', 'POS', 'REF', 'ALT']
        
        # 检查必需列
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"DataFrame缺少必需列: {missing_cols}")
        
        # 检查缺失值
        missing_data = df[required_cols].isnull().sum()
        if missing_data.any():
            logger.warning(f"数据中存在缺失值: \n{missing_data}")
        
        # 检查数据有效性
        invalid_pos = df[~df['POS'].between(1, 10**9)]  # 合理范围检查
        if not invalid_pos.empty:
            logger.warning(f"发现 {len(invalid_pos)} 个无效的POS值")
        
        # 检查REF/ALT是否为有效碱基
        valid_bases = set('ACGTN-')
        invalid_ref = []
        invalid_alt = []
        
        for idx, row in df.iterrows():
            ref = str(row['REF'])
            alt = str(row['ALT'])
            
            if not all(c in valid_bases for c in ref):
                invalid_ref.append(idx)
            if not all(c in valid_bases for c in alt):
                invalid_alt.append(idx)
        
        if invalid_ref:
            logger.warning(f"发现 {len(invalid_ref)} 个无效的REF值")
        if invalid_alt:
            logger.warning(f"发现 {len(invalid_alt)} 个无效的ALT值")
        
        # 严格模式下的额外检查
        if self.default_config.get('strict_mode', False):
            if invalid_ref or invalid_alt or not invalid_pos.empty:
                error_msgs = []
                if invalid_ref:
                    error_msgs.append(f"无效REF值: {len(invalid_ref)} 个")
                if invalid_alt:
                    error_msgs.append(f"无效ALT值: {len(invalid_alt)} 个")
                if not invalid_pos.empty:
                    error_msgs.append(f"无效POS值: {len(invalid_pos)} 个")
                
                raise ValueError(f"严格模式验证失败: {', '.join(error_msgs)}")
    
    def _add_variant_info(self, df: pd.DataFrame) -> pd.DataFrame:
        """添加变异信息列"""
        if df.empty:
            df['Is_Indel'] = False
            df['Variant_Type'] = ''
            df['Indel_Length'] = 0
            df['Indel_Size'] = 0
            return df
        
        df = df.copy()
        
        # 计算REF和ALT长度
        df['REF_Len'] = df['REF'].apply(lambda x: len(str(x)) if pd.notna(x) else 0)
        df['ALT_Len'] = df['ALT'].apply(lambda x: len(str(x)) if pd.notna(x) else 0)
        
        # 判断是否为indel
        df['Is_Indel'] = df['REF_Len'] != df['ALT_Len']
        
        # 确定变异类型
        def determine_variant_type(row):
            if pd.isna(row['REF']) or pd.isna(row['ALT']):
                return 'UNKNOWN'
            
            ref_len = row['REF_Len']
            alt_len = row['ALT_Len']
            
            if ref_len == alt_len:
                if ref_len == 1:
                    return 'SNP'
                else:
                    return 'MNP'
            elif ref_len > alt_len:
                return 'DELETION'
            elif alt_len > ref_len:
                return 'INSERTION'
            else:
                return 'COMPLEX'
        
        df['Variant_Type'] = df.apply(determine_variant_type, axis=1)
        
        # 计算indel长度（插入为正，删除为负）
        df['Indel_Length'] = df['ALT_Len'] - df['REF_Len']
        
        # 计算indel大小
        df['Indel_Size'] = df.apply(
            lambda row: max(row['REF_Len'], row['ALT_Len']) if row['Is_Indel'] else 0,
            axis=1
        )
        
        # 移除临时列
        df = df.drop(columns=['REF_Len', 'ALT_Len'])
        
        # 记录统计信息
        if not df.empty:
            indel_counts = df['Variant_Type'].value_counts()
            for var_type, count in indel_counts.items():
                logger.info(f"{var_type}: {count} 个")
        
        return df
    
    def filter_indels(self, df: pd.DataFrame, 
                     min_length: Optional[int] = None,
                     max_length: Optional[int] = None,
                     variant_types: Optional[List[str]] = None) -> pd.DataFrame:
        """
        过滤indel变异
        
        Args:
            df: 变异DataFrame
            min_length: 最小indel长度
            max_length: 最大indel长度
            variant_types: 允许的变异类型列表
            
        Returns:
            DataFrame: 过滤后的indel DataFrame
        """
        if df.empty:
            return df
        
        # 确保已经处理过indel
        if 'Is_Indel' not in df.columns:
            df = self._add_variant_info(df)
        
        # 使用配置值或参数值
        min_len = min_length or self.default_config.get('min_indel_length', 1)
        max_len = max_length or self.default_config.get('max_indel_length', 100)
        allowed_types = variant_types or self.default_config.get('allowed_variant_types', ['INSERTION', 'DELETION'])
        
        # 初步过滤：只保留indel
        filtered = df[df['Is_Indel']].copy()
        
        if filtered.empty:
            logger.warning("没有找到indel变异")
            return filtered
        
        # 按长度过滤
        length_filter = filtered['Indel_Size'].between(min_len, max_len)
        filtered = filtered[length_filter]
        
        # 按类型过滤
        if allowed_types:
            if isinstance(allowed_types, str):
                allowed_types = [allowed_types]
            
            type_filter = filtered['Variant_Type'].isin(allowed_types)
            filtered = filtered[type_filter]
        
        logger.info(f"过滤后保留 {len(filtered)} 个indel "
                   f"(长度范围: {min_len}-{max_len}bp)")
        
        return filtered
    
    def to_variant_records(self, df: pd.DataFrame) -> List[VariantRecord]:
        """
        将DataFrame转换为VariantRecord列表
        
        Args:
            df: 变异DataFrame
            
        Returns:
            List[VariantRecord]: VariantRecord对象列表
        """
        records = []
        
        for _, row in df.iterrows():
            record = VariantRecord(
                chrom=str(row.get('CHR', '')),
                pos=int(row.get('POS', 0)),
                ref=str(row.get('REF', '')),
                alt=str(row.get('ALT', '')),
                id=str(row.get('ID', '')) if 'ID' in row and pd.notna(row['ID']) else None,
                qual=float(row['QUAL']) if 'QUAL' in row and pd.notna(row['QUAL']) else None,
                filter=str(row['FILTER']) if 'FILTER' in row and pd.notna(row['FILTER']) else None,
                info=str(row['INFO']) if 'INFO' in row and pd.notna(row['INFO']) else None,
            )
            records.append(record)
        
        return records
    
    def extract_flanking_sequences(self, df: pd.DataFrame, 
                                 reference_fasta: str,
                                 flanking_size: int = 300) -> pd.DataFrame:
        """
        提取indel的侧翼序列
        
        Args:
            df: indel DataFrame
            reference_fasta: 参考基因组FASTA文件路径
            flanking_size: 每侧提取的序列长度
            
        Returns:
            DataFrame: 包含侧翼序列的DataFrame
        """
        if df.empty:
            return df
        
        try:
            from pyfaidx import Fasta
        except ImportError:
            logger.error("需要pyfaidx库来提取序列。请安装: pip install pyfaidx")
            raise
        
        logger.info(f"提取侧翼序列，每侧 {flanking_size}bp")
        
        genome = Fasta(reference_fasta)
        results = []
        
        for _, row in df.iterrows():
            try:
                chrom = str(row['CHR'])
                pos = int(row['POS'])
                ref = str(row['REF'])
                ref_len = len(ref)
                
                # 检查染色体是否存在
                if chrom not in genome:
                    logger.warning(f"染色体 {chrom} 不在参考基因组中")
                    continue
                
                # 计算序列范围
                chr_length = len(genome[chrom])
                start = max(0, pos - flanking_size - 1)  # 转换为0-based
                end = min(chr_length, pos + ref_len + flanking_size - 1)
                
                # 提取序列
                sequence = str(genome[chrom][start:end]).upper()
                
                # 计算indel在序列中的位置
                indel_start = (pos - 1) - start
                indel_end = indel_start + ref_len
                
                results.append({
                    'ID': row.get('ID', f"{chrom}_{pos}"),
                    'CHR': chrom,
                    'POS': pos,
                    'REF': ref,
                    'ALT': row['ALT'],
                    'Flanking_Sequence': sequence,
                    'Indel_Start_In_Seq': indel_start,
                    'Indel_End_In_Seq': indel_end,
                    'Flanking_Size': flanking_size,
                    'Sequence_Length': len(sequence),
                })
                
            except Exception as e:
                logger.error(f"提取序列失败 {row.get('ID', 'Unknown')}: {e}")
                continue
        
        return pd.DataFrame(results)


# 便捷函数
def parse_variants(file_path: Union[str, Path], 
                  file_format: Optional[str] = None,
                  filter_indels: bool = True,
                  **kwargs) -> pd.DataFrame:
    """
    解析变异文件的便捷函数
    
    Args:
        file_path: 文件路径
        file_format: 文件格式
        filter_indels: 是否只返回indel
        **kwargs: 传递给解析器的参数
        
    Returns:
        DataFrame: 变异数据
    """
    parser = VariantParser()
    df = parser.parse(file_path, file_format, **kwargs)
    
    if filter_indels and not df.empty:
        df = parser.filter_indels(df, **kwargs)
    
    return df


def read_variants_from_dataframe(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    从现有的DataFrame读取变异数据
    
    Args:
        df: 包含变异数据的DataFrame
        **kwargs: 传递给处理函数的参数
        
    Returns:
        DataFrame: 处理后的变异数据
    """
    parser = VariantParser()
    
    # 标准化DataFrame
    df_processed = parser._standardize_dataframe(df)
    
    # 验证
    parser._validate_dataframe(df_processed)
    
    # 处理indel
    df_processed = parser._add_variant_info(df_processed)
    
    # 过滤indel（如果要求）
    if kwargs.get('filter_indels', False):
        df_processed = parser.filter_indels(df_processed, **kwargs)
    
    return df_processed


def validate_variant_file(file_path: Union[str, Path]) -> Tuple[bool, List[str]]:
    """
    验证变异文件
    
    Args:
        file_path: 文件路径
        
    Returns:
        Tuple[bool, List[str]]: (是否有效, 错误信息列表)
    """
    errors = []
    
    file_path = Path(file_path)
    
    # 检查文件是否存在
    if not file_path.exists():
        errors.append(f"文件不存在: {file_path}")
        return False, errors
    
    # 检查文件是否可读
    try:
        with open(file_path, 'r') as f:
            f.readline()
    except Exception as e:
        errors.append(f"文件不可读: {e}")
        return False, errors
    
    # 尝试解析文件
    try:
        parser = VariantParser({'strict_mode': False})
        df = parser.parse(file_path, require_indel=False)
        
        if df.empty:
            errors.append("文件为空或没有有效数据")
            return False, errors
        
        # 检查是否有indel
        if 'Is_Indel' in df.columns:
            indel_count = df['Is_Indel'].sum()
            if indel_count == 0:
                errors.append("文件中没有找到indel变异")
                # 这不一定是个错误，取决于用户需求
        
        return True, errors
        
    except Exception as e:
        errors.append(f"解析文件失败: {e}")
        return False, errors
