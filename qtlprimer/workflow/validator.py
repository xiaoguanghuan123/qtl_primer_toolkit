"""
输入验证器
验证工作流输入文件的完整性和有效性。
"""
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Any
import logging

from ..core.utils import (
    validate_genome_file,
    validate_blast_database,
    FileUtils,
    SequenceUtils
)

logger = logging.getLogger(__name__)


class InputValidator:
    """输入文件验证器"""
    
    def __init__(self, config: Optional[Dict] = None):
        """
        初始化验证器
        
        Args:
            config: 配置字典
        """
        self.config = config or {}
        
        # 默认验证规则
        self.default_rules = {
            'variants_file': {
                'required_extensions': ['.vcf', '.csv', '.tsv', '.txt', '.vcf.gz'],
                'max_size_mb': 100,  # 最大文件大小
                'min_size_kb': 1,    # 最小文件大小
            },
            'reference_fasta': {
                'required_extensions': ['.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz'],
                'max_size_gb': 10,   # 最大文件大小
                'min_size_kb': 10,   # 最小文件大小
            },
            'blast_db': {
                'required_files': ['.nhr', '.nin', '.nsq'],
            }
        }
        
        logger.info("输入验证器初始化完成")
    
    def validate_workflow_inputs(
        self,
        variants_file: Optional[str] = None,
        reference_fasta: Optional[str] = None,
        blast_db: Optional[str] = None
    ) -> Tuple[bool, List[str], List[str]]:
        """
        验证工作流输入文件
        
        Args:
            variants_file: 变异文件路径
            reference_fasta: 参考基因组路径
            blast_db: BLAST数据库路径
            
        Returns:
            Tuple[bool, List[str], List[str]]: (是否有效, 错误列表, 警告列表)
        """
        errors = []
        warnings = []
        
        # 验证变异文件
        if variants_file:
            is_valid, file_errors, file_warnings = self.validate_variants_file(variants_file)
            if not is_valid:
                errors.extend(file_errors)
            warnings.extend(file_warnings)
        else:
            errors.append("未提供变异文件路径")
        
        # 验证参考基因组
        if reference_fasta:
            is_valid, file_errors, file_warnings = self.validate_reference_genome(reference_fasta)
            if not is_valid:
                errors.extend(file_errors)
            warnings.extend(file_warnings)
        else:
            errors.append("未提供参考基因组路径")
        
        # 验证BLAST数据库
        if blast_db:
            is_valid, file_errors, file_warnings = self.validate_blast_database(blast_db)
            if not is_valid:
                errors.extend(file_errors)
            warnings.extend(file_warnings)
        else:
            errors.append("未提供BLAST数据库路径")
        
        # 检查文件关联性
        if not errors:
            self._check_file_compatibility(variants_file, reference_fasta, warnings)
        
        return len(errors) == 0, errors, warnings
    
    def validate_variants_file(
        self, 
        file_path: Union[str, Path]
    ) -> Tuple[bool, List[str], List[str]]:
        """
        验证变异文件
        
        Args:
            file_path: 文件路径
            
        Returns:
            Tuple[bool, List[str], List[str]]: (是否有效, 错误列表, 警告列表)
        """
        errors = []
        warnings = []
        
        file_path = Path(file_path)
        
        # 1. 检查文件是否存在
        if not file_path.exists():
            errors.append(f"变异文件不存在: {file_path}")
            return False, errors, warnings
        
        # 2. 检查文件大小
        file_size_mb = file_path.stat().st_size / (1024 * 1024)
        max_size = self.default_rules['variants_file']['max_size_mb']
        
        if file_size_mb > max_size:
            warnings.append(f"变异文件较大 ({file_size_mb:.1f}MB > {max_size}MB)")
        
        # 3. 检查文件扩展名
        valid_extensions = self.default_rules['variants_file']['required_extensions']
        if file_path.suffix.lower() not in valid_extensions:
            warnings.append(f"变异文件扩展名不常见: {file_path.suffix}")
        
        # 4. 尝试读取文件内容
        try:
            content = FileUtils.read_file(file_path)
            lines = content.strip().split('\n')
            
            if not lines:
                errors.append("变异文件为空")
                return False, errors, warnings
            
            # 检查是否包含变异数据的关键词
            vcf_keywords = ['#CHROM', 'CHROM', 'POS', 'REF', 'ALT']
            csv_keywords = ['CHR', 'POS', 'REF', 'ALT']
            
            first_line = lines[0].upper()
            has_vcf = any(keyword in first_line for keyword in vcf_keywords)
            has_csv = any(keyword in first_line for keyword in csv_keywords)
            
            if not (has_vcf or has_csv):
                # 尝试检查其他行
                found_keywords = False
                for line in lines[:10]:
                    line_upper = line.upper()
                    if any(keyword in line_upper for keyword in vcf_keywords + csv_keywords):
                        found_keywords = True
                        break
                
                if not found_keywords:
                    warnings.append("变异文件可能不包含标准格式的变异数据")
            
            # 统计行数
            data_lines = [line for line in lines if line and not line.startswith('#')]
            if len(data_lines) == 0:
                warnings.append("变异文件没有数据行")
            else:
                logger.info(f"变异文件包含 {len(data_lines)} 个数据行")
            
            # 检查VCF格式
            if file_path.suffix.lower() in ['.vcf', '.vcf.gz']:
                vcf_errors = self._validate_vcf_format(content)
                if vcf_errors:
                    warnings.extend(vcf_errors)
        
        except Exception as e:
            errors.append(f"读取变异文件失败: {str(e)}")
            return False, errors, warnings
        
        return True, errors, warnings
    
    def validate_reference_genome(
        self, 
        file_path: Union[str, Path]
    ) -> Tuple[bool, List[str], List[str]]:
        """
        验证参考基因组文件
        
        Args:
            file_path: 文件路径
            
        Returns:
            Tuple[bool, List[str], List[str]]: (是否有效, 错误列表, 警告列表)
        """
        errors = []
        warnings = []
        
        file_path = Path(file_path)
        
        # 1. 检查文件是否存在
        if not file_path.exists():
            errors.append(f"参考基因组文件不存在: {file_path}")
            return False, errors, warnings
        
        # 2. 检查文件大小
        file_size_gb = file_path.stat().st_size / (1024 * 1024 * 1024)
        max_size = self.default_rules['reference_fasta']['max_size_gb']
        
        if file_size_gb > max_size:
            warnings.append(f"参考基因组文件较大 ({file_size_gb:.2f}GB > {max_size}GB)")
        
        # 3. 检查文件扩展名
        valid_extensions = self.default_rules['reference_fasta']['required_extensions']
        if file_path.suffix.lower() not in valid_extensions:
            warnings.append(f"参考基因组文件扩展名不常见: {file_path.suffix}")
        
        # 4. 验证FASTA格式
        is_valid, fasta_errors = validate_genome_file(file_path)
        if not is_valid:
            errors.extend(fasta_errors)
        
        if not errors:
            # 5. 尝试读取序列信息
            try:
                sequences = FileUtils.read_fasta(file_path)
                
                if not sequences:
                    errors.append("参考基因组文件为空")
                    return False, errors, warnings
                
                # 统计信息
                total_sequences = len(sequences)
                total_length = sum(len(seq) for seq in sequences.values())
                
                logger.info(f"参考基因组包含 {total_sequences} 条序列，总长度 {total_length:,} bp")
                
                # 检查序列有效性
                invalid_sequences = []
                short_sequences = []
                
                for header, seq in sequences.items():
                    if not SequenceUtils.is_valid_dna(seq):
                        invalid_sequences.append(header[:50])
                    
                    if len(seq) < 100:
                        short_sequences.append(f"{header[:50]}: {len(seq)}bp")
                
                if invalid_sequences:
                    warnings.append(f"发现 {len(invalid_sequences)} 条无效DNA序列")
                
                if short_sequences:
                    warnings.append(f"发现 {len(short_sequences)} 条短序列 (<100bp)")
                
                # 检查染色体命名
                chr_patterns = ['chr', 'chromosome', 'scaffold', 'contig']
                has_standard_names = False
                
                for header in sequences.keys():
                    header_lower = header.lower()
                    if any(pattern in header_lower for pattern in chr_patterns):
                        has_standard_names = True
                        break
                
                if not has_standard_names:
                    warnings.append("参考基因组序列名称可能不是标准染色体命名")
            
            except Exception as e:
                errors.append(f"读取参考基因组文件失败: {str(e)}")
        
        return len(errors) == 0, errors, warnings
    
    def validate_blast_database(
        self, 
        db_path: Union[str, Path]
    ) -> Tuple[bool, List[str], List[str]]:
        """
        验证BLAST数据库
        
        Args:
            db_path: 数据库路径
            
        Returns:
            Tuple[bool, List[str], List[str]]: (是否有效, 错误列表, 警告列表)
        """
        errors = []
        warnings = []
        
        db_path = Path(db_path)
        
        # 1. 检查数据库文件是否存在
        required_files = self.default_rules['blast_db']['required_files']
        missing_files = []
        
        for ext in required_files:
            db_file = f"{db_path}{ext}"
            if not Path(db_file).exists():
                missing_files.append(db_file)
        
        if missing_files:
            errors.append(f"BLAST数据库文件缺失: {', '.join(missing_files)}")
            return False, errors, warnings
        
        # 2. 验证数据库有效性
        is_valid, db_errors = validate_blast_database(db_path)
        if not is_valid:
            errors.extend(db_errors)
        
        if not errors:
            # 3. 检查数据库大小
            total_size_mb = 0
            for ext in required_files:
                db_file = f"{db_path}{ext}"
                total_size_mb += Path(db_file).stat().st_size / (1024 * 1024)
            
            logger.info(f"BLAST数据库大小: {total_size_mb:.1f} MB")
            
            if total_size_mb < 1:
                warnings.append("BLAST数据库较小，可能不完整")
        
        return len(errors) == 0, errors, warnings
    
    def _validate_vcf_format(self, content: str) -> List[str]:
        """
        验证VCF格式
        
        Args:
            content: VCF文件内容
            
        Returns:
            List[str]: 格式错误列表
        """
        warnings = []
        lines = content.strip().split('\n')
        
        # 查找标题行
        header_line = None
        for line in lines:
            if line.startswith('#CHROM'):
                header_line = line
                break
        
        if not header_line:
            warnings.append("VCF文件缺少#CHROM标题行")
            return warnings
        
        # 解析标题行
        headers = header_line[1:].split('\t')  # 移除#号
        required_headers = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        
        missing_headers = []
        for required in required_headers:
            if required not in headers:
                missing_headers.append(required)
        
        if missing_headers:
            warnings.append(f"VCF文件缺少必需列: {', '.join(missing_headers)}")
        
        return warnings
    
    def _check_file_compatibility(
        self, 
        variants_file: Union[str, Path],
        reference_fasta: Union[str, Path],
        warnings: List[str]
    ):
        """
        检查文件之间的兼容性
        
        Args:
            variants_file: 变异文件路径
            reference_fasta: 参考基因组路径
            warnings: 警告列表（会被修改）
        """
        try:
            # 读取参考基因组的染色体名称
            ref_sequences = FileUtils.read_fasta(reference_fasta)
            ref_chromosomes = set(seq[:50].upper() for seq in ref_sequences.keys())  # 取前50个字符
            
            # 读取变异文件的染色体名称
            variants_content = FileUtils.read_file(variants_file)
            variants_lines = variants_content.strip().split('\n')
            
            variant_chromosomes = set()
            for line in variants_lines:
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 1:
                        variant_chromosomes.add(parts[0].upper())
            
            # 检查染色体名称匹配
            if variant_chromosomes:
                mismatched = []
                for v_chr in variant_chromosomes:
                    found = False
                    for r_chr in ref_chromosomes:
                        if v_chr in r_chr or r_chr in v_chr:
                            found = True
                            break
                    
                    if not found:
                        # 尝试标准化染色体名称
                        v_chr_clean = v_chr.replace('CHR', '').replace('CHROMOSOME', '').strip()
                        for r_chr in ref_chromosomes:
                            r_chr_clean = r_chr.replace('CHR', '').replace('CHROMOSOME', '').strip()
                            if v_chr_clean == r_chr_clean:
                                found = True
                                break
                    
                    if not found:
                        mismatched.append(v_chr)
                
                if mismatched:
                    warnings.append(f"变异文件中的染色体在参考基因组中未找到: {', '.join(mismatched[:5])}")
        
        except Exception as e:
            warnings.append(f"检查文件兼容性时出错: {str(e)}")
    
    def validate_output_directory(
        self,
        output_dir: Union[str, Path],
        overwrite: bool = False
    ) -> Tuple[bool, List[str], List[str]]:
        """
        验证输出目录
        
        Args:
            output_dir: 输出目录路径
            overwrite: 是否允许覆盖
            
        Returns:
            Tuple[bool, List[str], List[str]]: (是否有效, 错误列表, 警告列表)
        """
        errors = []
        warnings = []
        
        output_dir = Path(output_dir)
        
        # 1. 检查目录是否存在
        if output_dir.exists():
            if not output_dir.is_dir():
                errors.append(f"输出路径不是目录: {output_dir}")
                return False, errors, warnings
            
            # 2. 检查目录是否为空
            try:
                contents = list(output_dir.iterdir())
                if contents and not overwrite:
                    warnings.append(f"输出目录不为空: {output_dir}")
                    
                    # 检查是否有重要文件
                    important_files = ['.csv', '.xlsx', '.html', '.png']
                    has_important = any(
                        any(file.suffix.lower() == ext for ext in important_files)
                        for file in contents
                    )
                    
                    if has_important:
                        errors.append(f"输出目录包含重要文件，请使用不同的目录或启用覆盖模式")
                        return False, errors, warnings
            except Exception as e:
                warnings.append(f"检查输出目录内容时出错: {e}")
        
        # 3. 检查写入权限
        try:
            test_file = output_dir / ".write_test"
            test_file.touch()
            test_file.unlink()
        except Exception as e:
            errors.append(f"输出目录不可写: {e}")
        
        return len(errors) == 0, errors, warnings
    
    def validate_primers_dataframe(self, primers_df) -> Tuple[bool, List[str], List[str]]:
        """
        验证引物DataFrame
        
        Args:
            primers_df: 引物DataFrame
            
        Returns:
            Tuple[bool, List[str], List[str]]: (是否有效, 错误列表, 警告列表)
        """
        errors = []
        warnings = []
        
        if primers_df is None or primers_df.empty:
            warnings.append("引物DataFrame为空")
            return len(errors) == 0, errors, warnings
        
        # 检查必需列
        required_columns = ['ID', 'LEFT_PRIMER', 'RIGHT_PRIMER', 'PRODUCT_SIZE']
        missing_columns = []
        
        for col in required_columns:
            if col not in primers_df.columns:
                missing_columns.append(col)
        
        if missing_columns:
            errors.append(f"引物DataFrame缺少必需列: {', '.join(missing_columns)}")
        
        # 检查引物序列
        invalid_primers = []
        for idx, row in primers_df.iterrows():
            if 'LEFT_PRIMER' in primers_df.columns and 'LEFT_PRIMER' in row:
                left_primer = str(row['LEFT_PRIMER'])
                if left_primer and left_primer not in ['FAILED', '']:
                    if not SequenceUtils.is_valid_dna(left_primer):
                        invalid_primers.append(f"行{idx} 左引物: {left_primer[:20]}...")
            
            if 'RIGHT_PRIMER' in primers_df.columns and 'RIGHT_PRIMER' in row:
                right_primer = str(row['RIGHT_PRIMER'])
                if right_primer and right_primer not in ['FAILED', '']:
                    if not SequenceUtils.is_valid_dna(right_primer):
                        invalid_primers.append(f"行{idx} 右引物: {right_primer[:20]}...")
        
        if invalid_primers:
            warnings.append(f"发现 {len(invalid_primers)} 个无效引物序列")
        
        return len(errors) == 0, errors, warnings


@dataclass
class ValidationResult:
    """验证结果"""
    is_valid: bool
    errors: List[str]
    warnings: List[str]
    details: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self):
        """转换为字典"""
        return {
            'is_valid': self.is_valid,
            'errors': self.errors,
            'warnings': self.warnings,
            'details': self.details
        }
    
    def to_string(self, include_warnings: bool = True) -> str:
        """转换为字符串"""
        lines = []
        lines.append("=" * 60)
        lines.append("输入文件验证结果")
        lines.append("=" * 60)
        
        if self.is_valid:
            lines.append("✅ 验证通过")
        else:
            lines.append("❌ 验证失败")
        
        if self.errors:
            lines.append("")
            lines.append("错误:")
            for error in self.errors:
                lines.append(f"  - {error}")
        
        if include_warnings and self.warnings:
            lines.append("")
            lines.append("警告:")
            for warning in self.warnings:
                lines.append(f"  - {warning}")
        
        if self.details:
            lines.append("")
            lines.append("详细信息:")
            for key, value in self.details.items():
                lines.append(f"  {key}: {value}")
        
        lines.append("=" * 60)
        return "\n".join(lines)


# 便捷函数
def validate_workflow_inputs(
    variants_file: Optional[str] = None,
    reference_fasta: Optional[str] = None,
    blast_db: Optional[str] = None
) -> Tuple[bool, List[str], List[str]]:
    """
    验证工作流输入文件（便捷函数）
    
    Args:
        variants_file: 变异文件路径
        reference_fasta: 参考基因组路径
        blast_db: BLAST数据库路径
        
    Returns:
        Tuple[bool, List[str], List[str]]: (是否有效, 错误列表, 警告列表)
    """
    validator = InputValidator()
    return validator.validate_workflow_inputs(variants_file, reference_fasta, blast_db)


def validate_all_inputs(
    variants_file: str,
    reference_fasta: str,
    blast_db: str,
    output_dir: Optional[str] = None
) -> ValidationResult:
    """
    验证所有输入文件和输出目录
    
    Args:
        variants_file: 变异文件路径
        reference_fasta: 参考基因组路径
        blast_db: BLAST数据库路径
        output_dir: 输出目录路径（可选）
        
    Returns:
        ValidationResult: 验证结果
    """
    validator = InputValidator()
    
    # 验证输入文件
    is_valid, errors, warnings = validator.validate_workflow_inputs(
        variants_file, reference_fasta, blast_db
    )
    
    details = {}
    
    # 验证输出目录
    if output_dir:
        output_is_valid, output_errors, output_warnings = validator.validate_output_directory(
            output_dir, overwrite=False
        )
        
        if not output_is_valid:
            is_valid = False
            errors.extend(output_errors)
        
        warnings.extend(output_warnings)
    
    # 收集详细信息
    try:
        # 文件大小信息
        if Path(variants_file).exists():
            size_mb = Path(variants_file).stat().st_size / (1024 * 1024)
            details['variants_file_size'] = f"{size_mb:.1f} MB"
        
        if Path(reference_fasta).exists():
            size_gb = Path(reference_fasta).stat().st_size / (1024 * 1024 * 1024)
            details['reference_genome_size'] = f"{size_gb:.2f} GB"
    except Exception as e:
        warnings.append(f"收集文件信息时出错: {e}")
    
    return ValidationResult(
        is_valid=is_valid,
        errors=errors,
        warnings=warnings,
        details=details
    )
