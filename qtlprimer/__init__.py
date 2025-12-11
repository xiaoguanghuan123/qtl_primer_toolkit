"""QTL Primer Toolkit - A comprehensive toolkit for designing PCR primers."""

__version__ = "1.0.0"
__author__ = "Guangchao Sun & Changchuang Liu"
__email__ = "guangchaosun@sicau.edu.cn"
__license__ = "MIT"

import logging
from .utils.logging_utils import setup_logging

# 设置默认日志
setup_logging()

# 主模块导入
from .core.variant_parser import VariantParser
from .core.primer_designer import PrimerDesigner
from .core.specificity_checker import SpecificityChecker
from .visualization.gel_electrophoresis import GelElectrophoresisVisualizer
from .workflow.manager import QTLFineMappingWorkflow

# 便捷函数
def run_workflow(variants_file, reference_fasta, blast_db, **kwargs):
    """一键运行完整工作流"""
    from .workflow.manager import run_qtl_fine_mapping_pipeline
    return run_qtl_fine_mapping_pipeline(
        variants_file=variants_file,
        reference_fasta=reference_fasta,
        blast_db=blast_db,
        **kwargs
    )

def design_primers(variants_file, reference_fasta, **kwargs):
    """仅设计引物"""
    from .core.primer_designer import design_primers_from_dataframe
    from .core.variant_parser import parse_variants
    
    df = parse_variants(variants_file)
    return design_primers_from_dataframe(df, reference_fasta, **kwargs)

__all__ = [
    'VariantParser',
    'PrimerDesigner',
    'SpecificityChecker',
    'GelElectrophoresisVisualizer',
    'QTLFineMappingWorkflow',
    'run_workflow',
    'design_primers',
]
