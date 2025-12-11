"""
QTL Primer Toolkit 工作流模块
提供完整的工作流管理、报告生成和输入验证功能。
"""

from .manager import (
    QTLFineMappingWorkflow,
    run_qtl_fine_mapping_pipeline,
    WorkflowConfig,
    WorkflowResult
)

from .reporter import (
    ReportGenerator,
    generate_excel_report,
    generate_html_report,
    generate_markdown_report
)

from .validator import (
    InputValidator,
    validate_workflow_inputs,
    ValidationResult
)

__version__ = "1.0.0"
__author__ = "QTL Primer Toolkit Development Team"

# 模块初始化日志
import logging
logger = logging.getLogger(__name__)
logger.info(f"QTL Primer Toolkit Workflow Module v{__version__} initialized")

# 导出主要类和函数
__all__ = [
    # 工作流管理
    'QTLFineMappingWorkflow',
    'run_qtl_fine_mapping_pipeline',
    'WorkflowConfig',
    'WorkflowResult',
    
    # 报告生成
    'ReportGenerator',
    'generate_excel_report',
    'generate_html_report',
    'generate_markdown_report',
    
    # 输入验证
    'InputValidator',
    'validate_workflow_inputs',
    'ValidationResult',
]

# 便捷函数
def run_standard_workflow(variants_file, reference_fasta, blast_db, output_dir=None):
    """
    运行标准工作流（简化参数）
    
    Args:
        variants_file: 变异文件路径
        reference_fasta: 参考基因组路径
        blast_db: BLAST数据库路径
        output_dir: 输出目录（可选）
        
    Returns:
        WorkflowResult: 工作流结果
    """
    from .manager import run_qtl_fine_mapping_pipeline
    
    if output_dir is None:
        import os
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = f"./qtl_results_{timestamp}"
    
    return run_qtl_fine_mapping_pipeline(
        variants_file=variants_file,
        reference_fasta=reference_fasta,
        blast_db=blast_db,
        output_dir=output_dir
    )

def validate_and_run(variants_file, reference_fasta, blast_db, output_dir=None):
    """
    验证输入文件后运行工作流
    
    Args:
        variants_file: 变异文件路径
        reference_fasta: 参考基因组路径
        blast_db: BLAST数据库路径
        output_dir: 输出目录（可选）
        
    Returns:
        tuple: (是否成功, 结果或错误信息)
    """
    # 先验证输入
    from .validator import validate_workflow_inputs
    
    is_valid, errors, warnings = validate_workflow_inputs(
        variants_file=variants_file,
        reference_fasta=reference_fasta,
        blast_db=blast_db
    )
    
    if not is_valid:
        logger.error(f"输入验证失败: {errors}")
        return False, errors
    
    if warnings:
        logger.warning(f"输入验证警告: {warnings}")
    
    # 运行工作流
    try:
        result = run_standard_workflow(variants_file, reference_fasta, blast_db, output_dir)
        return True, result
    except Exception as e:
        logger.error(f"工作流执行失败: {e}")
        return False, str(e)

def generate_quick_report(workflow_result, output_dir=None, formats=None):
    """
    为工作流结果快速生成报告
    
    Args:
        workflow_result: WorkflowResult对象
        output_dir: 输出目录（可选）
        formats: 报告格式列表（可选）
        
    Returns:
        dict: 报告文件路径字典
    """
    from .reporter import ReportGenerator
    
    if output_dir is None:
        output_dir = workflow_result.output_dir / "reports"
    
    if formats is None:
        formats = ['excel', 'html', 'markdown']
    
    reporter = ReportGenerator()
    reports = reporter.generate_all_reports(workflow_result, output_dir, formats)
    
    return reports
