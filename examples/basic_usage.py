#!/usr/bin/env python3
"""
基本使用示例
"""
import sys
from pathlib import Path

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from qtlprimer import run_workflow

def main():
    """基本使用示例"""
    
    # 示例数据路径
    variants_file = "examples/data/test_indels.vcf"
    reference_fasta = "examples/data/test_genome.fasta"
    blast_db = "examples/data/test_blastdb"
    
    # 检查文件是否存在
    for file_path in [variants_file, reference_fasta]:
        if not Path(file_path).exists():
            print(f"警告: 文件不存在: {file_path}")
            print("请先创建示例数据文件")
            return
    
    try:
        # 运行工作流
        print("🚀 运行QTL精细定位引物设计工作流...")
        
        results = run_workflow(
            variants_file=variants_file,
            reference_fasta=reference_fasta,
            blast_db=blast_db,
            output_dir="./example_results"
        )
        
        if results:
            print("✅ 工作流完成!")
            print(f"输出目录: {results.get('output_dir')}")
            print(f"找到Indel数量: {results.get('indel_count', 0)}")
            print(f"设计引物对数: {results.get('primer_count', 0)}")
        else:
            print("❌ 工作流失败")
    
    except Exception as e:
        print(f"❌ 错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
