"""
重组图可视化模块
用于展示QTL区间、标记位置和重组事件的遗传学图表。
"""
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path
import logging

from ..config.defaults import VISUALIZATION_DEFAULTS

logger = logging.getLogger(__name__)


class RecombinationPlotter:
    """重组事件可视化器"""
    
    def __init__(self, config: Optional[Dict] = None):
        """
        初始化可视化器
        
        Args:
            config: 配置字典，覆盖默认值
        """
        self.config = VISUALIZATION_DEFAULTS.copy()
        if config:
            self.config.update(config)
        
        # 图形设置
        self.figsize = self.config.get('FIG_SIZE', (10, 6))
        self.dpi = self.config.get('DPI', 300)
        
        # 颜色方案
        self.colors = self.config.get('COLORS', {
            'chromosome': '#6c757d',
            'qtl_region': '#ffc107',
            'marker': '#007bff',
            'recombinant': '#dc3545',
            'non_recombinant': '#28a745',
            'heterozygous': '#6f42c1',
            'text': '#212529',
        })
        
        logger.info("重组图可视化器初始化完成")
    
    def plot_recombination_interval(
        self,
        chromosome: str,
        positions: List[int],
        marker_names: Optional[List[str]] = None,
        genotypes: Optional[List[str]] = None,
        qtl_position: Optional[int] = None,
        output_file: Optional[Union[str, Path]] = None
    ) -> plt.Figure:
        """
        绘制重组区间图
        
        Args:
            chromosome: 染色体名称
            positions: 标记位置列表
            marker_names: 标记名称列表
            genotypes: 基因型列表
            qtl_position: QTL推测位置
            output_file: 输出文件路径
            
        Returns:
            matplotlib Figure对象
        """
        if len(positions) < 2:
            raise ValueError("至少需要两个标记位置")
        
        # 设置默认值
        if marker_names is None:
            marker_names = [f"M{i+1}" for i in range(len(positions))]
        if genotypes is None:
            genotypes = ["Homozygous"] * len(positions)
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=self.figsize, dpi=self.dpi,
                                       gridspec_kw={'height_ratios': [1, 2]})
        
        # 1. 绘制染色体和标记位置
        self._draw_chromosome_map(ax1, chromosome, positions, marker_names, qtl_position)
        
        # 2. 绘制重组事件示意图
        self._draw_recombination_diagram(ax2, positions, marker_names, genotypes)
        
        # 设置总标题
        fig.suptitle(f"QTL精细定位重组分析 - 染色体 {chromosome}", 
                    fontsize=14, fontweight='bold', y=0.95)
        
        plt.tight_layout()
        
        # 保存或显示
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(output_file), dpi=self.dpi, bbox_inches='tight')
            logger.info(f"重组图已保存至: {output_file}")
        
        return fig
    
    def _draw_chromosome_map(self, ax, chromosome, positions, marker_names, qtl_position):
        """绘制染色体和标记位置"""
        # 染色体表示
        chr_length = max(positions) * 1.2  # 染色体长度
        chr_y = 0.5
        
        # 绘制染色体
        ax.add_patch(plt.Rectangle((0, chr_y - 0.05), chr_length, 0.1,
                                  facecolor=self.colors['chromosome'], alpha=0.3,
                                  edgecolor=self.colors['chromosome'], linewidth=2))
        
        # 绘制着丝粒示意
        centromere_pos = chr_length / 2
        ax.add_patch(plt.Rectangle((centromere_pos - chr_length*0.02, chr_y - 0.08),
                                  chr_length*0.04, 0.16,
                                  facecolor='gray', alpha=0.5))
        
        # 绘制标记位置
        for pos, name in zip(positions, marker_names):
            ax.plot([pos, pos], [chr_y - 0.1, chr_y - 0.2], 
                   color=self.colors['marker'], linewidth=1.5)
            ax.scatter(pos, chr_y - 0.2, color=self.colors['marker'], 
                      s=100, zorder=5)
            ax.text(pos, chr_y - 0.25, f"{name}\n{pos:,} bp",
                   ha='center', va='top', fontsize=9)
        
        # 绘制QTL位置（如果提供）
        if qtl_position:
            ax.axvline(x=qtl_position, ymin=chr_y-0.05, ymax=chr_y+0.05,
                      color=self.colors['qtl_region'], linewidth=3, linestyle='--', alpha=0.7)
            ax.text(qtl_position, chr_y + 0.15, "QTL区域",
                   ha='center', va='bottom', fontsize=10,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor=self.colors['qtl_region'], alpha=0.3))
        
        # 设置坐标轴
        ax.set_xlim(0, chr_length)
        ax.set_ylim(0, 1)
        ax.set_xlabel("染色体位置 (bp)", fontsize=10)
        ax.set_ylabel("染色体示意图", fontsize=10)
        ax.set_title(f"染色体 {chromosome} - 标记位置", fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # 隐藏y轴刻度
        ax.set_yticks([])
    
    def _draw_recombination_diagram(self, ax, positions, marker_names, genotypes):
        """绘制重组事件示意图"""
        num_markers = len(positions)
        
        # 创建基因型颜色映射
        genotype_colors = {
            'AA': self.colors['non_recombinant'],
            'BB': self.colors['non_recombinant'],
            'aa': self.colors['non_recombinant'],
            'bb': self.colors['non_recombinant'],
            'Aa': self.colors['heterozygous'],
            'Bb': self.colors['heterozygous'],
            'AB': self.colors['recombinant'],
            'ab': self.colors['recombinant'],
            'Homozygous': self.colors['non_recombinant'],
            'Heterozygous': self.colors['heterozygous'],
            'Recombinant': self.colors['recombinant'],
        }
        
        # 绘制标记和基因型
        for i, (pos, name, genotype) in enumerate(zip(positions, marker_names, genotypes)):
            y_pos = num_markers - i  # 从上到下排列
            
            # 绘制标记点
            color = genotype_colors.get(genotype, self.colors['marker'])
            ax.scatter(pos, y_pos, color=color, s=200, zorder=5, edgecolor='black', linewidth=1)
            
            # 标注标记名称和基因型
            ax.text(pos, y_pos + 0.3, f"{name}\n{genotype}",
                   ha='center', va='bottom', fontsize=9,
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
        
        # 绘制连接线表示可能的重组事件
        if len(positions) >= 2:
            for i in range(len(positions) - 1):
                y_pos1 = num_markers - i
                y_pos2 = num_markers - (i + 1)
                
                # 检查是否可能发生重组
                geno1 = genotypes[i]
                geno2 = genotypes[i + 1]
                
                is_recombinant = self._is_potential_recombinant(geno1, geno2)
                
                if is_recombinant:
                    # 绘制重组线（红色虚线）
                    ax.plot([positions[i], positions[i+1]], [y_pos1, y_pos2],
                           color=self.colors['recombinant'], linestyle='--', 
                           linewidth=2, alpha=0.7)
                    
                    # 标注重组事件
                    mid_x = (positions[i] + positions[i+1]) / 2
                    mid_y = (y_pos1 + y_pos2) / 2
                    
                    ax.annotate('重组事件', xy=(mid_x, mid_y),
                               xytext=(mid_x, mid_y - 0.5),
                               arrowprops=dict(arrowstyle='->', color=self.colors['recombinant']),
                               ha='center', va='top', fontsize=9,
                               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=self.colors['recombinant']))
                else:
                    # 绘制非重组线（绿色实线）
                    ax.plot([positions[i], positions[i+1]], [y_pos1, y_pos2],
                           color=self.colors['non_recombinant'], linestyle='-',
                           linewidth=2, alpha=0.5)
        
        # 设置坐标轴
        ax.set_xlim(min(positions) * 0.9, max(positions) * 1.1)
        ax.set_ylim(0, num_markers + 1)
        ax.set_xlabel("染色体位置 (bp)", fontsize=10)
        ax.set_ylabel("标记", fontsize=10)
        ax.set_title("重组事件分析", fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # 设置y轴刻度为标记名称
        ax.set_yticks(range(1, num_markers + 1))
        ax.set_yticklabels(reversed(marker_names))
        
        # 添加图例
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=self.colors['non_recombinant'], label='纯合/非重组'),
            Patch(facecolor=self.colors['heterozygous'], label='杂合'),
            Patch(facecolor=self.colors['recombinant'], label='重组事件'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    def _is_potential_recombinant(self, geno1: str, geno2: str) -> bool:
        """判断两个标记之间是否可能发生重组"""
        # 简化判断：一个纯合一个杂合可能表示重组
        def is_homozygous(genotype):
            genotype = str(genotype).upper()
            return genotype in ['AA', 'BB', 'AA/AA', 'BB/BB', 'HOMOZYGOUS']
        
        def is_heterozygous(genotype):
            genotype = str(genotype).upper()
            return genotype in ['AB', 'AB/AB', 'Aa', 'Bb', 'Aa/Aa', 'Bb/Bb', 'HETEROZYGOUS']
        
        return (is_homozygous(geno1) and is_heterozygous(geno2)) or \
               (is_heterozygous(geno1) and is_homozygous(geno2))
    
    def create_genetic_map(
        self,
        chromosomes: List[str],
        positions: List[List[int]],
        lod_scores: Optional[List[List[float]]] = None,
        qtl_regions: Optional[List[Dict]] = None,
        output_file: Optional[Union[str, Path]] = None
    ) -> plt.Figure:
        """
        创建遗传连锁图
        
        Args:
            chromosomes: 染色体列表
            positions: 每个染色体的标记位置列表
            lod_scores: LOD分数列表（可选）
            qtl_regions: QTL区域信息列表（可选）
            output_file: 输出文件路径
            
        Returns:
            matplotlib Figure对象
        """
        fig, axes = plt.subplots(len(chromosomes), 1, figsize=(10, 3*len(chromosomes)), 
                                dpi=self.dpi, squeeze=False)
        axes = axes.flatten()
        
        for idx, (chrom, chrom_positions) in enumerate(zip(chromosomes, positions)):
            ax = axes[idx]
            
            # 绘制染色体标记
            ax.plot(chrom_positions, np.ones_like(chrom_positions) * 0.5, 
                   'o', color=self.colors['marker'], markersize=8, alpha=0.7)
            
            # 绘制LOD分数（如果提供）
            if lod_scores and idx < len(lod_scores):
                chrom_lod = lod_scores[idx]
                if len(chrom_lod) == len(chrom_positions):
                    # 创建LOD分数柱状图
                    bars = ax.bar(chrom_positions, chrom_lod, 
                                 width=np.diff(sorted(chrom_positions)).min() * 0.3 if len(chrom_positions) > 1 else 10,
                                 alpha=0.5, color='orange')
                    
                    # 标注显著LOD分数
                    for pos, lod in zip(chrom_positions, chrom_lod):
                        if lod > 3.0:  # 常用LOD阈值
                            ax.text(pos, lod + 0.1, f"{lod:.1f}", 
                                   ha='center', va='bottom', fontsize=8)
            
            # 绘制QTL区域（如果提供）
            if qtl_regions:
                for qtl in qtl_regions:
                    if qtl.get('chromosome') == chrom:
                        start = qtl.get('start', min(chrom_positions))
                        end = qtl.get('end', max(chrom_positions))
                        
                        ax.axvspan(start, end, alpha=0.2, color=self.colors['qtl_region'])
                        
                        # 标注QTL
                        qtl_pos = (start + end) / 2
                        ax.text(qtl_pos, ax.get_ylim()[1] * 0.9, "QTL",
                               ha='center', va='top', fontsize=10,
                               bbox=dict(boxstyle='round,pad=0.3', 
                                        facecolor=self.colors['qtl_region'], alpha=0.5))
            
            # 设置坐标轴
            ax.set_xlim(min(chrom_positions) * 0.9, max(chrom_positions) * 1.1)
            ax.set_ylim(0, max(ax.get_ylim()[1], 10))
            ax.set_xlabel(f"染色体 {chrom} 位置 (cM)", fontsize=9)
            ax.set_ylabel("LOD分数", fontsize=9)
            ax.set_title(f"染色体 {chrom} 遗传连锁图", fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
        
        plt.suptitle("遗传连锁图与QTL定位", fontsize=14, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # 保存或显示
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(output_file), dpi=self.dpi, bbox_inches='tight')
            logger.info(f"遗传连锁图已保存至: {output_file}")
        
        return fig


# 便捷函数
def plot_recombination_interval(
    chromosome: str,
    positions: List[int],
    marker_names: Optional[List[str]] = None,
    genotypes: Optional[List[str]] = None,
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict] = None
) -> plt.Figure:
    """
    绘制重组区间图（便捷函数）
    
    Args:
        chromosome: 染色体名称
        positions: 标记位置列表
        marker_names: 标记名称列表
        genotypes: 基因型列表
        output_file: 输出文件路径
        config: 配置字典
        
    Returns:
        matplotlib Figure对象
    """
    plotter = RecombinationPlotter(config)
    
    fig = plotter.plot_recombination_interval(
        chromosome=chromosome,
        positions=positions,
        marker_names=marker_names,
        genotypes=genotypes,
        output_file=output_file
    )
    
    return fig


def create_genetic_map(
    chromosomes: List[str],
    positions: List[List[int]],
    lod_scores: Optional[List[List[float]]] = None,
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict] = None
) -> plt.Figure:
    """
    创建遗传连锁图（便捷函数）
    
    Args:
        chromosomes: 染色体列表
        positions: 每个染色体的标记位置列表
        lod_scores: LOD分数列表
        output_file: 输出文件路径
        config: 配置字典
        
    Returns:
        matplotlib Figure对象
    """
    plotter = RecombinationPlotter(config)
    
    fig = plotter.create_genetic_map(
        chromosomes=chromosomes,
        positions=positions,
        lod_scores=lod_scores,
        output_file=output_file
    )
    
    return fig
