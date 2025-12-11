#!/usr/bin/env python3
"""
QTL精细定位引物设计工具包 - 主入口点
"""
import argparse
import sys
import logging
from pathlib import Path

from primer_toolkit.workflow.workflow_manager import run_qtl_fine_mapping_pipeline

def setup_logging(verbose: bool):
    """设置日志"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
        ]
    )

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='QTL精细定位引物设计工具包',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  %(prog)s design --variants my_indels.vcf --reference genome.fa --blast-db genome_db
  %(prog)s design --variants variants.csv --reference genome.fa --output ./my_results
  
输出文件:
  - indels_filtered.csv: 过滤后的候选Indel
  - primers_designed.csv: 设计的引物对
  - primers_with_specificity.csv: 特异性验证和排名
  - gel_visualizations/: 电泳结果模拟图
  - qtl_primers_full_report.xlsx: Excel完整报告
  - qtl_primers_report.html: HTML交互报告
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='子命令')
    
    # design 子命令
    parser_design = subparsers.add_parser('design', help='运行完整工作流')
    parser_design.add_argument('--variants', '-v', required=True,
                              help='变异文件路径 (VCF/CSV)')
    parser_design.add_argument('--reference', '-r', required=True,
                              help='参考基因组FASTA文件路径')
    parser_design.add_argument('--blast-db', '-b', required=True,
                              help='BLAST数据库路径')
    parser_design.add_argument('--output', '-o', default='./qtl_results',
                              help='输出目录 (默认: ./qtl_results)')
    parser_design.add_argument('--sample', '-s', default='QTL',
                              help='样本前缀 (默认: QTL)')
    parser_design.add_argument('--verbose', '-V', action='store_true',
                              help='显示详细日志')
    
    # visualize 子命令
    parser_viz = subparsers.add_parser('visualize', help='仅生成电泳可视化')
    parser_viz.add_argument('--primer-pair', '-p', required=True,
                           help='引物对ID')
    parser_viz.add_argument('--size-a', type=int, default=350,
                           help='标记A产物大小 (默认: 350bp)')
    parser_viz.add_argument('--size-b', type=int, default=500,
                           help='标记B产物大小 (默认: 500bp)')
    parser_viz.add_argument('--output', '-o', default='./gel_images',
                           help='输出目录 (默认: ./gel_images)')
    
    args = parser.parse_args()
    
    if args.command == 'design':
        setup_logging(args.verbose)
        
        try:
            print(f"""
            ============================================
            QTL精细定位引物设计工具包
            ============================================
            输入文件:
              - 变异文件: {args.variants}
              - 参考基因组: {args.reference}
              - BLAST数据库: {args.blast_db}
            输出目录: {args.output}
            ============================================
            """)
            
            # 运行完整工作流
            results = run_qtl_fine_mapping_pipeline(
                variants_file=args.variants,
                reference_fasta=args.reference,
                blast_db=args.blast_db,
                output_dir=args.output
            )
            
            if results and results.get('output_dir'):
                print(f"""
            ============================================
            工作流完成！
            ============================================
            主要结果文件:
              - Indel筛选结果: {results['output_dir']}/01_indels_filtered.csv
              - 引物设计结果: {results['output_dir']}/02_primers_designed.csv
              - 特异性排名: {results['output_dir']}/03_primers_with_specificity.csv
              - 电泳模拟图: {results['output_dir']}/04_gel_visualizations/
              - 完整报告: {results['output_dir']}/qtl_primers_full_report.xlsx
              
            请查看 {results['output_dir']}/README.md 获取详细说明
            ============================================
                """)
            
        except Exception as e:
            print(f"错误: {str(e)}", file=sys.stderr)
            sys.exit(1)
    
    elif args.command == 'visualize':
        from primer_toolkit.visualization import visualize_recombination_scenario
        
        output_dir = Path(args.output)
        output_dir.mkdir(exist_ok=True)
        
        # 生成三种场景的可视化
        scenarios = [
            ('AA', 'BB', 'non_recombinant'),
            ('AA', 'Bb', 'recombinant'),
            ('Aa', 'Bb', 'double_heterozygous')
        ]
        
        for geno_a, geno_b, scenario in scenarios:
            output_file = output_dir / f"{args.primer_pair}_{scenario}.png"
            
            visualize_recombination_scenario(
                marker_a_genotype=geno_a,
                marker_b_genotype=geno_b,
                size_a=args.size_a,
                size_b=args.size_b,
                output_file=str(output_file)
            )
            
            print(f"已生成: {output_file}")
        
        print(f"\n所有可视化文件已保存至: {output_dir}")
    
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
