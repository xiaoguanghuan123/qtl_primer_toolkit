"""
QTL Primer Toolkit Core Module
核心算法模块：变异解析、引物设计、特异性验证
"""

from .variant_parser import VariantParser, parse_variants, read_variants_from_dataframe
from .primer_designer import PrimerDesigner, design_primers_from_dataframe
from .specificity_checker import SpecificityChecker, check_primers_specificity
from .utils import (
    validate_genome_file,
    validate_blast_database,
    calculate_tm,
    calculate_gc_content,
    reverse_complement
)

__all__ = [
    # 变异解析
    'VariantParser',
    'parse_variants',
    'read_variants_from_dataframe',
    
    # 引物设计
    'PrimerDesigner',
    'design_primers_from_dataframe',
    
    # 特异性验证
    'SpecificityChecker',
    'check_primers_specificity',
    
    # 工具函数
    'validate_genome_file',
    'validate_blast_database',
    'calculate_tm',
    'calculate_gc_content',
    'reverse_complement',
]

# 版本信息
__version__ = "1.0.0"
__author__ = "QTL Primer Toolkit Development Team"
__email__ = "guagnchaosun@sicau.edu.cn"

# 模块初始化日志
import logging
logger = logging.getLogger(__name__)
logger.info(f"QTL Primer Toolkit Core Module v{__version__} initialized")
