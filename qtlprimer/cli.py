#!/usr/bin/env python3
"""
QTL Primer Toolkit - 命令行接口
"""
import argparse
import sys
import logging
from pathlib import Path

from qtlprimer import __version__, __author__
from qtlprimer.utils.logging_utils import setup_logging
from qtlprimer.workflow.manager import run_qtl_fine_mapping_pipeline
from qtlprimer.visualization.gel_electrophoresis import (
    GelElectrophoresisVisualizer, 
    visualize_recombination_scenario
)

def main():
    """主命令行入口"""
    parser = argparse.ArgumentParser(
        prog='qtlprimer',
        description='QTL精细定位引物设计工具包 v{}'.format(__version__),
        epilog='作者: {} | 更多信息: https://github.com/yourusername/qtl-primer-toolkit'.format(__author__),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-V', '--version',
        action='version',
        version='%(prog)s v{}'.format(__version__)
    )
    
    # 创建子命令
    subparsers = parser.add_subparsers(
        dest='command',
        title='可用命令',
        metavar='<command>'
    )
    
    # ==================== design 命令 ====================
    parser_design = subparsers.add_parser(
        'design',
        help='运行完整引物设计工作流'
    )
    
    # 必需参数
    parser_design.add_argument(
        '-v', '--variants',
        required=True,
        help='变异文件路径 (VCF/CSV格式)'
    )
    parser_design.add_argument(
        '-r', '--reference',
        required=True,
        help='参考基因组FASTA文件路径'
    )
    parser_design.add_argument(
        '-b', '--blast-db',
        required=True,
        help='BLAST数据库路径'
    )
    
    # 可选参数
    parser_design.add_argument(
        '-o', '--output',
        default='./qtl_results',
        help='输出目录 (默认: ./qtl_results)'
    )
    parser_design.add_argument(
        '-s', '--sample',
        default='QTL',
        help='样本前缀 (默认: QTL)'
    )
    parser_design.add_argument(
        '-c', '--config',
        help='配置文件路径 (YAML/JSON)'
    )
    
    # 设计参数
    design_group = parser_design.add_argument_group('引物设计参数')
    design_group.add_argument(
        '--min-size',
        type=int,
        default=18,
        help='引物最小长度 (默认: 18)'
    )
    design_group.add_argument(
        '--max-size',
        type=int,
        default=24,
        help='引物最大长度 (默认: 24)'
    )
    design_group.add_argument(
        '--tm-range',
        type=float,
        nargs=2,
        default=[50.0, 60.0],
        help='熔解温度范围 (默认: 50.0 60.0)'
    )
    
    # 输出参数
    output_group = parser_design.add_argument_group('输出参数')
    output_group.add_argument(
        '--formats',
        nargs='+',
        choices=['csv', 'excel', 'html', 'all'],
        default=['csv', 'excel'],
        help='输出格式 (默认: csv excel)'
    )
    output_group.add_argument(
        '--top-n',
        type=int,
        default=10,
        help='显示前N个最佳引物 (默认: 10)'
    )
    
    # 其他参数
    parser_design.add_argument(
        '--verbose',
        action='store_true',
        help='显示详细日志'
    )
    parser_design.add_argument(
        '--dry-run',
        action='store_true',
        help='模拟运行，不生成实际文件'
    )
    
    # ==================== visualize 命令 ====================
    parser_viz = subparsers.add_parser(
        'visualize',
        help='生成电泳可视化'
    )
    
    parser_viz.add_argument(
        '-p', '--primer-pair',
        required=True,
        help='引物对ID或名称'
    )
    parser_viz.add_argument(
        '-s', '--sizes',
        type=int,
        nargs=2,
        required=True,
        metavar=('SIZE_A', 'SIZE_B'),
        help='两个标记的产物大小 (bp)'
    )
    parser_viz.add_argument(
        '-o', '--output',
        default='./gel_images',
        help='输出目录 (默认: ./gel_images)'
    )
    
    # 基因型选项
    genotype_group = parser_viz.add_argument_group('基因型选项')
    genotype_group.add_argument(
        '--scenario',
        choices=['all', 'non_recombinant', 'recombinant', 'double_het'],
        default='all',
        help='要生成的可视化场景 (默认: all)'
    )
    genotype_group.add_argument(
        '--genotype-a',
        default='AA',
        help='标记A基因型 (默认: AA)'
    )
    genotype_group.add_argument(
        '--genotype-b',
        default='BB',
        help='标记B基因型 (默认: BB)'
    )
    
    # ==================== check 命令 ====================
    parser_check = subparsers.add_parser(
        'check',
        help='检查输入文件有效性'
    )
    
    parser_check.add_argument(
        '-v', '--variants',
        help='变异文件路径'
    )
    parser_check.add_argument(
        '-r', '--reference',
        help='参考基因组文件路径'
    )
    parser_check.add_argument(
        '-b', '--blast-db',
        help='BLAST数据库路径'
    )
    parser_check.add_argument(
        '--all',
        action='store_true',
        help='检查所有必需文件'
    )
    
    # ==================== config 命令 ====================
    parser_config = subparsers.add_parser(
        'config',
        help='生成或验证配置文件'
    )
    
    parser_config.add_argument(
        'action',
        choices=['generate', 'validate', 'show'],
        help='操作: generate=生成, validate=验证, show=显示'
    )
    parser_config.add_argument(
        '-o', '--output',
        help='输出配置文件路径'
    )
    parser_config.add_argument(
        '-f', '--format',
        choices=['yaml', 'json', 'python'],
        default='yaml',
        help='配置文件格式 (默认: yaml)'
    )
    
    # ==================== 解析参数 ====================
    args = parser.parse_args()
    
    # 设置日志
    log_level = logging.DEBUG if getattr(args, 'verbose', False) else logging.INFO
    setup_logging(level=log_level)
    
    logger = logging.getLogger(__name__)
    
    # 处理命令
    if args.command == 'design':
        logger.info("=" * 60)
        logger.info("QTL精细定位引物设计工作流")
        logger.info("=" * 60)
        
        try:
            # 准备配置
            config = {}
            if args.config:
                import yaml
                with open(args.config, 'r') as f:
                    config = yaml.safe_load(f)
            
            # 更新命令行参数
            if args.min_size or args.max_size or args.tm_range:
                config.setdefault('primer3', {})
                if args.min_size:
                    config['primer3']['PRIMER_MIN_SIZE'] = args.min_size
                if args.max_size:
                    config['primer3']['PRIMER_MAX_SIZE'] = args.max_size
                if args.tm_range:
                    config['primer3']['PRIMER_MIN_TM'] = args.tm_range[0]
                    config['primer3']['PRIMER_MAX_TM'] = args.tm_range[1]
            
            if args.formats:
                config.setdefault('tool', {})
                config['tool']['OUTPUT_FORMATS'] = args.formats
            
            if args.top_n:
                config.setdefault('tool', {})
                config['tool']['TOP_N_RESULTS'] = args.top_n
            
            # 运行工作流
            if args.dry_run:
                logger.info("模拟运行 - 检查输入文件...")
                
                from qtlprimer.workflow.validator import validate_inputs
                is_valid, errors = validate_inputs(
                    variants_file=args.variants,
                    reference_fasta=args.reference,
                    blast_db=args.blast_db
                )
                
                if is_valid:
                    logger.info("✅ 输入文件验证通过")
                    logger.info(f"变异文件: {args.variants}")
                    logger.info(f"参考基因组: {args.reference}")
                    logger.info(f"BLAST数据库: {args.blast_db}")
                    logger.info(f"输出目录: {args.output}")
                    logger.info("模拟运行完成 - 可以正式运行")
                else:
                    logger.error("❌ 输入文件验证失败:")
                    for error in errors:
                        logger.error(f"  - {error}")
                    sys.exit(1)
            else:
                results = run_qtl_fine_mapping_pipeline(
                    variants_file=args.variants,
                    reference_fasta=args.reference,
                    blast_db=args.blast_db,
                    output_dir=args.output,
                    config=config
                )
                
                if results and results.get('output_dir'):
                    logger.info("=" * 60)
                    logger.info("🎉 工作流完成!")
                    logger.info("=" * 60)
                    logger.info(f"输出目录: {results['output_dir']}")
                    logger.info(f"找到Indel数量: {results.get('indel_count', 0)}")
                    logger.info(f"设计引物对数: {results.get('primer_count', 0)}")
                    logger.info("")
                    logger.info("📁 生成的文件:")
                    logger.info(f"  - 候选Indel: {results['output_dir']}/01_indels_filtered.csv")
                    logger.info(f"  - 引物设计: {results['output_dir']}/02_primers_designed.csv")
                    logger.info(f"  - 特异性验证: {results['output_dir']}/03_primers_with_specificity.csv")
                    logger.info(f"  - 电泳模拟图: {results['output_dir']}/04_gel_visualizations/")
                    logger.info(f"  - 完整报告: {results['output_dir']}/qtl_primers_full_report.xlsx")
                    logger.info(f"  - HTML报告: {results['output_dir']}/qtl_primers_report.html")
                    logger.info("")
                    logger.info("📖 详细说明请查看 README.md")
                    logger.info("=" * 60)
        
        except Exception as e:
            logger.error(f"工作流执行失败: {str(e)}")
            if getattr(args, 'verbose', False):
                import traceback
                traceback.print_exc()
            sys.exit(1)
    
    elif args.command == 'visualize':
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        visualizer = GelElectrophoresisVisualizer()
        
        scenarios = []
        if args.scenario == 'all':
            scenarios = [
                ('AA', 'BB', 'non_recombinant'),
                ('AA', 'Bb', 'recombinant'),
                ('Aa', 'Bb', 'double_heterozygous')
            ]
        elif args.scenario == 'non_recombinant':
            scenarios = [(args.genotype_a, args.genotype_b, 'non_recombinant')]
        elif args.scenario == 'recombinant':
            scenarios = [(args.genotype_a, args.genotype_b, 'recombinant')]
        elif args.scenario == 'double_het':
            scenarios = [(args.genotype_a, args.genotype_b, 'double_heterozygous')]
        
        for geno_a, geno_b, scenario in scenarios:
            output_file = output_dir / f"{args.primer_pair}_{scenario}.png"
            
            visualize_recombination_scenario(
                marker_a_genotype=geno_a,
                marker_b_genotype=geno_b,
                size_a=args.sizes[0],
                size_b=args.sizes[1],
                output_file=str(output_file),
                sample_name=args.primer_pair
            )
            
            print(f"✅ 已生成: {output_file}")
        
        print(f"\n🎨 所有可视化文件已保存至: {output_dir}")
    
    elif args.command == 'check':
        from qtlprimer.workflow.validator import validate_inputs
        
        files_to_check = {}
        if args.variants:
            files_to_check['variants_file'] = args.variants
        if args.reference:
            files_to_check['reference_fasta'] = args.reference
        if args.blast_db:
            files_to_check['blast_db'] = args.blast_db
        
        if args.all or not files_to_check:
            # 检查所有必需文件
            files_to_check = {
                'variants_file': args.variants,
                'reference_fasta': args.reference,
                'blast_db': args.blast_db
            }
        
        is_valid, errors = validate_inputs(**files_to_check)
        
        if is_valid:
            print("✅ 所有输入文件验证通过")
        else:
            print("❌ 输入文件验证失败:")
            for error in errors:
                print(f"  - {error}")
            sys.exit(1)
    
    elif args.command == 'config':
        from qtlprimer.config.defaults import get_config, validate_config
        
        if args.action == 'generate':
            config = get_config()
            
            if args.format == 'yaml':
                import yaml
                output = yaml.dump(config, default_flow_style=False, sort_keys=False)
            elif args.format == 'json':
                import json
                output = json.dumps(config, indent=2, ensure_ascii=False)
            elif args.format == 'python':
                import pprint
                output = pprint.pformat(config, indent=2)
            
            if args.output:
                with open(args.output, 'w') as f:
                    f.write(output)
                print(f"✅ 配置文件已生成: {args.output}")
            else:
                print(output)
        
        elif args.action == 'validate':
            if args.output:
                import yaml
                with open(args.output, 'r') as f:
                    config = yaml.safe_load(f)
                
                is_valid, errors = validate_config(config)
                
                if is_valid:
                    print("✅ 配置文件验证通过")
                else:
                    print("❌ 配置文件验证失败:")
                    for error in errors:
                        print(f"  - {error}")
                    sys.exit(1)
            else:
                print("❌ 请使用 -o 指定要验证的配置文件")
                sys.exit(1)
        
        elif args.action == 'show':
            config = get_config()
            import yaml
            print(yaml.dump(config, default_flow_style=False, sort_keys=False))
    
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == '__main__':
    main()
