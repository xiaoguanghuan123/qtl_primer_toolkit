"""
QTL Primer Toolkit 可视化模块
用于生成F2:3群体重组检测的电泳模拟图和遗传学可视化。
"""

from .gel_electrophoresis import (
    GelElectrophoresisVisualizer,
    visualize_recombination_scenario,
    create_tutorial_plot,
)
from .recombination_plots import (
    RecombinationPlotter,
    plot_recombination_interval,
    create_genetic_map,
)

__version__ = "1.0.0"
__author__ = "QTL Primer Toolkit Development Team"

# 模块初始化日志
import logging
logger = logging.getLogger(__name__)
logger.info(f"QTL Primer Toolkit Visualization Module v{__version__} initialized")

# 导出主要类和函数
__all__ = [
    # 凝胶电泳可视化
    'GelElectrophoresisVisualizer',
    'visualize_recombination_scenario',
    'create_tutorial_plot',
    
    # 重组图可视化
    'RecombinationPlotter',
    'plot_recombination_interval',
    'create_genetic_map',
]

# 便捷函数
def create_all_visualizations(primer_pairs, output_dir, config=None):
    """
    为所有引物对创建完整的可视化
    
    Args:
        primer_pairs: 引物对数据列表
        output_dir: 输出目录
        config: 配置字典
        
    Returns:
        dict: 可视化结果字典
    """
    import os
    from pathlib import Path
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # 创建凝胶电泳可视化
    gel_dir = output_dir / "gel_electrophoresis"
    gel_dir.mkdir(exist_ok=True)
    
    visualizer = GelElectrophoresisVisualizer(config=config)
    
    for i, primer_pair in enumerate(primer_pairs):
        primer_id = primer_pair.get('id', f"primer_{i}")
        
        # 创建电泳图
        gel_results = visualizer.create_primer_visualizations(
            primer_id=primer_id,
            product_sizes=primer_pair.get('product_sizes', [300, 500]),
            output_dir=gel_dir
        )
        
        results[f"{primer_id}_gel"] = gel_results
    
    # 创建重组图（如果提供了位置信息）
    plotter = RecombinationPlotter(config=config)
    
    for i, primer_pair in enumerate(primer_pairs):
        if 'positions' in primer_pair and 'chromosome' in primer_pair:
            primer_id = primer_pair.get('id', f"primer_{i}")
            
            # 创建重组区间图
            plot_file = output_dir / f"{primer_id}_recombination_plot.png"
            plotter.plot_recombination_interval(
                chromosome=primer_pair['chromosome'],
                positions=primer_pair['positions'],
                output_file=str(plot_file)
            )
            
            results[f"{primer_id}_plot"] = str(plot_file)
    
    logger.info(f"创建了 {len(results)} 个可视化文件")
    return results
