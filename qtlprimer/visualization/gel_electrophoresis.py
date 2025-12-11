"""
凝胶电泳可视化模块
模拟F2:3群体中Indel标记的PCR电泳结果。
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path
import logging

from ..config.defaults import VISUALIZATION_DEFAULTS

logger = logging.getLogger(__name__)


class GelElectrophoresisVisualizer:
    """琼脂糖凝胶电泳结果模拟器"""
    
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
        self.figsize = self.config.get('FIG_SIZE', (12, 8))
        self.dpi = self.config.get('DPI', 300)
        self.font_size = self.config.get('FONT_SIZE', 10)
        
        # 颜色方案
        self.colors = self.config.get('COLORS', {
            'band': '#2E86AB',
            'hetero_band': '#A23B72',
            'gel_background': '#F8F9FA',
            'well': '#6C757D',
            'marker': '#495057',
            'text': '#212529',
            'good': '#28a745',
            'warning': '#ffc107',
            'error': '#dc3545',
        })
        
        # 电泳参数
        self.gel_params = self.config.get('GEL_PARAMS', {
            'band_height': 0.6,
            'band_width_scale': 0.8,
            'hetero_spacing': 0.15,
            'blur_sigma': 0.7,
            'band_intensity': 0.9,
        })
        
        # 基因型示例
        self.genotype_examples = self.config.get('GENOTYPE_EXAMPLES', {
            'non_recombinant': ('AA', 'BB'),
            'recombinant': ('AA', 'Bb'),
            'double_heterozygous': ('Aa', 'Bb'),
        })
        
        logger.info("凝胶电泳可视化器初始化完成")
    
    def visualize_recombination_test(
        self, 
        marker_a_genotype: str, 
        marker_b_genotype: str,
        marker_a_product_size: int,
        marker_b_product_size: int,
        sample_name: str = "F2 Plant",
        output_path: Optional[Union[str, Path]] = None
    ) -> plt.Figure:
        """
        可视化重组检测的电泳结果
        
        Args:
            marker_a_genotype: 标记A基因型 ('AA', 'Aa', 'aa')
            marker_b_genotype: 标记B基因型 ('BB', 'Bb', 'bb')
            marker_a_product_size: 标记A PCR产物大小 (bp)
            marker_b_product_size: 标记B PCR产物大小 (bp)
            sample_name: 样本名称
            output_path: 输出文件路径，如果为None则显示图像
            
        Returns:
            matplotlib Figure对象
        """
        # 1. 创建图形
        fig, axes = self._create_gel_layout()
        
        # 2. 确定重组状态
        recombination_status = self._determine_recombination_status(
            marker_a_genotype, marker_b_genotype
        )
        
        # 3. 绘制Marker泳道
        self._draw_marker_lane(axes['gel'])
        
        # 4. 绘制样本泳道
        self._draw_sample_lane(
            axes['gel'], 
            marker_a_genotype, 
            marker_b_genotype,
            marker_a_product_size, 
            marker_b_product_size,
            sample_name
        )
        
        # 5. 添加解释文本
        self._add_explanation_text(
            axes['text'], 
            marker_a_genotype, 
            marker_b_genotype,
            recombination_status,
            marker_a_product_size,
            marker_b_product_size
        )
        
        # 6. 设置标题
        title_text = self._generate_title(recombination_status, sample_name)
        fig.suptitle(title_text, fontsize=14, fontweight='bold', y=0.98)
        
        plt.tight_layout()
        
        # 保存或显示
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(output_path), dpi=self.dpi, bbox_inches='tight')
            logger.info(f"电泳图已保存至: {output_path}")
        
        return fig
    
    def _create_gel_layout(self) -> Dict[str, plt.Axes]:
        """创建凝胶图的布局"""
        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        
        # 创建网格布局：凝胶图占70%，解释文本占30%
        gs = fig.add_gridspec(1, 2, width_ratios=[0.7, 0.3], wspace=0.05)
        
        axes = {
            'gel': fig.add_subplot(gs[0]),
            'text': fig.add_subplot(gs[1])
        }
        
        # 设置凝胶图背景和样式
        axes['gel'].set_facecolor(self.colors['gel_background'])
        axes['gel'].set_xlim(-1, 3)
        axes['gel'].set_ylim(0, 100)
        axes['gel'].invert_yaxis()  # 电泳方向从上到下
        axes['gel'].set_xlabel('泳道', fontsize=self.font_size)
        axes['gel'].set_ylabel('片段大小 (bp)', fontsize=self.font_size)
        axes['gel'].grid(True, alpha=0.3, linestyle='--')
        axes['gel'].set_axisbelow(True)
        
        # 隐藏文本区域的坐标轴
        axes['text'].axis('off')
        
        return {'fig': fig, **axes}
    
    def _determine_recombination_status(self, geno_a: str, geno_b: str) -> Dict:
        """根据两个标记的基因型确定重组状态"""
        # 简化判断：大写字母表示纯合，小写字母表示杂合中的另一个等位基因
        def is_homozygous(genotype):
            if not genotype or len(genotype) < 2:
                return False
            return genotype[0] == genotype[1] or genotype.isupper()
        
        is_homozygous_a = is_homozygous(geno_a)
        is_homozygous_b = is_homozygous(geno_b)
        
        # 判断逻辑
        if is_homozygous_a and is_homozygous_b:
            status = "no_recombination"
            explanation = "两个标记均为纯合，未检测到重组"
        elif (is_homozygous_a and not is_homozygous_b) or (not is_homozygous_a and is_homozygous_b):
            status = "recombination_detected"
            explanation = "一个标记纯合，一个标记杂合，检测到重组事件"
        else:
            status = "double_heterozygous"
            explanation = "两个标记均为杂合，需要进一步分析F3家系"
        
        return {
            'status': status,
            'explanation': explanation,
            'is_homozygous_a': is_homozygous_a,
            'is_homozygous_b': is_homozygous_b
        }
    
    def _draw_marker_lane(self, ax: plt.Axes):
        """绘制DNA Marker泳道"""
        # 常见DNA Marker的片段大小
        marker_sizes = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        
        # 绘制上样孔
        well_x = 0.5
        ax.add_patch(patches.Rectangle(
            (well_x - 0.15, -2), 0.3, 1.5,
            facecolor=self.colors['well'], alpha=0.8
        ))
        
        # 绘制Marker条带
        for size in marker_sizes:
            if size % 100 == 0:  # 每100bp显示一个主要条带
                y_pos = size / 10  # 缩放以适配坐标轴
                
                # 条带宽度随片段大小变化
                width = self.gel_params['band_width_scale'] * (1 - size/2000)
                
                ax.add_patch(patches.Rectangle(
                    (well_x - width/2, y_pos - self.gel_params['band_height']/2),
                    width, self.gel_params['band_height'],
                    facecolor=self.colors['marker'], alpha=0.9,
                    edgecolor='black', linewidth=0.5
                ))
                
                # 标注主要条带
                if size in [100, 500, 1000]:
                    ax.text(well_x + 0.2, y_pos, f'{size} bp',
                           fontsize=8, va='center')
        
        ax.text(well_x, -4, 'DNA\nMarker', fontsize=9, 
                ha='center', va='top', fontweight='bold')
    
    def _draw_sample_lane(
        self, 
        ax: plt.Axes,
        geno_a: str,
        geno_b: str,
        size_a: int,
        size_b: int,
        sample_name: str
    ):
        """绘制样本泳道"""
        well_x = 1.5
        
        # 绘制上样孔
        ax.add_patch(patches.Rectangle(
            (well_x - 0.15, -2), 0.3, 1.5,
            facecolor=self.colors['well'], alpha=0.8
        ))
        
        # 绘制标记A的条带
        self._draw_bands_for_marker(
            ax, well_x, geno_a, size_a, 'A', offset=-0.2
        )
        
        # 绘制标记B的条带
        self._draw_bands_for_marker(
            ax, well_x, geno_b, size_b, 'B', offset=0.2
        )
        
        # 标注泳道
        ax.text(well_x, -4, sample_name, fontsize=9,
                ha='center', va='top', fontweight='bold')
    
    def _draw_bands_for_marker(
        self, 
        ax: plt.Axes,
        x_center: float,
        genotype: str,
        product_size: int,
        marker_label: str,
        offset: float
    ):
        """为单个标记绘制条带"""
        y_pos = product_size / 10  # 缩放
        
        # 判断基因型
        def is_homozygous(genotype):
            if not genotype or len(genotype) < 2:
                return True  # 默认视为纯合
            return genotype[0] == genotype[1] or genotype.isupper() or genotype.lower() in ['aa', 'bb']
        
        if is_homozygous(genotype):
            # 纯合：单一条带
            self._draw_single_band(
                ax, x_center + offset, y_pos, 
                f"Marker {marker_label}\n{genotype}\n{product_size}bp",
                is_heterozygous=False
            )
        else:
            # 杂合：双带（模拟Indel杂合子的特征）
            self._draw_heterozygous_bands(
                ax, x_center + offset, y_pos, product_size,
                f"Marker {marker_label}\n{genotype}\n~{product_size}bp"
            )
    
    def _draw_single_band(
        self, 
        ax: plt.Axes,
        x: float, 
        y: float,
        label: str,
        is_heterozygous: bool = False
    ):
        """绘制单一条带"""
        color = self.colors['hetero_band'] if is_heterozygous else self.colors['band']
        
        # 主条带
        band = patches.Rectangle(
            (x - self.gel_params['band_width_scale']/2, 
             y - self.gel_params['band_height']/2),
            self.gel_params['band_width_scale'], 
            self.gel_params['band_height'],
            facecolor=color, alpha=self.gel_params['band_intensity'],
            edgecolor='black', linewidth=1,
            zorder=5
        )
        ax.add_patch(band)
        
        # 添加模糊效果（模拟真实电泳）
        for i in range(3):
            blur_alpha = 0.3 - i * 0.08
            blur_height = self.gel_params['band_height'] + 0.1 * (i + 1)
            blur_band = patches.Rectangle(
                (x - self.gel_params['band_width_scale']/2 * (1 - i*0.05),
                 y - blur_height/2),
                self.gel_params['band_width_scale'] * (1 - i*0.05),
                blur_height,
                facecolor=color, alpha=blur_alpha,
                edgecolor='none', zorder=4-i
            )
            ax.add_patch(blur_band)
        
        # 添加标签
        ax.text(x + self.gel_params['band_width_scale']/2 + 0.1, y, 
               label, fontsize=8, va='center')
    
    def _draw_heterozygous_bands(
        self, 
        ax: plt.Axes,
        x: float, 
        y: float,
        product_size: int,
        label: str
    ):
        """绘制杂合基因型的双带"""
        spacing = self.gel_params['hetero_spacing']
        
        # 第一条带（略小片段）
        self._draw_single_band(
            ax, x - spacing/2, y * 0.98, "", is_heterozygous=True
        )
        
        # 第二条带（略大片段）
        self._draw_single_band(
            ax, x + spacing/2, y * 1.02, "", is_heterozygous=True
        )
        
        # 添加总标签
        ax.text(x + spacing/2 + 0.15, y, label, fontsize=8, va='center')
    
    def _add_explanation_text(
        self,
        ax: plt.Axes,
        geno_a: str,
        geno_b: str,
        recombination_status: Dict,
        size_a: int,
        size_b: int
    ):
        """添加解释文本"""
        y_pos = 0.95
        line_height = 0.08
        
        # 标题
        ax.text(0.05, y_pos, '结果解释', fontsize=12, 
                fontweight='bold', transform=ax.transAxes)
        y_pos -= line_height * 2
        
        # 基因型信息
        ax.text(0.05, y_pos, f'标记A基因型: {geno_a}', fontsize=10,
                transform=ax.transAxes)
        y_pos -= line_height
        
        ax.text(0.05, y_pos, f'标记B基因型: {geno_b}', fontsize=10,
                transform=ax.transAxes)
        y_pos -= line_height
        
        ax.text(0.05, y_pos, f'产物大小: {size_a}bp / {size_b}bp', 
                fontsize=10, transform=ax.transAxes)
        y_pos -= line_height * 1.5
        
        # 重组状态
        status_display = {
            'no_recombination': '❌ 未检测到重组',
            'recombination_detected': '✅ 检测到重组',
            'double_heterozygous': '⚠️ 双杂合需进一步分析'
        }
        
        status_text = status_display.get(
            recombination_status['status'], 
            recombination_status['explanation']
        )
        
        ax.text(0.05, y_pos, f'重组状态:', fontsize=10,
                fontweight='bold', transform=ax.transAxes)
        y_pos -= line_height
        
        ax.text(0.08, y_pos, status_text, fontsize=10,
                style='italic', transform=ax.transAxes)
        y_pos -= line_height * 1.5
        
        # 详细解释
        ax.text(0.05, y_pos, '遗传学解释:', fontsize=10,
                fontweight='bold', transform=ax.transAxes)
        y_pos -= line_height
        
        explanation_lines = self._wrap_text(
            recombination_status['explanation'], 40
        )
        
        for line in explanation_lines:
            ax.text(0.08, y_pos, line, fontsize=9,
                    transform=ax.transAxes)
            y_pos -= line_height
        
        # F2:3群体特异性说明
        y_pos -= line_height
        ax.text(0.05, y_pos, 'F2:3群体提示:', fontsize=10,
                fontweight='bold', transform=ax.transAxes)
        y_pos -= line_height
        
        f2_notes = [
            "• F2植株的PCR结果需通过F3家系验证",
            "• 重组单株的后代将呈现分离",
            "• 用于精细定位QTL位置"
        ]
        
        for note in f2_notes:
            ax.text(0.08, y_pos, note, fontsize=8,
                    transform=ax.transAxes)
            y_pos -= line_height * 0.8
    
    def _generate_title(self, recombination_status: Dict, sample_name: str) -> str:
        """生成图表标题"""
        base_title = f"F2:3群体Indel标记电泳模拟 - {sample_name}"
        
        if recombination_status['status'] == 'recombination_detected':
            return f"{base_title} - ✅ 重组单株检测"
        elif recombination_status['status'] == 'no_recombination':
            return f"{base_title} - ❌ 非重组单株"
        else:
            return f"{base_title} - ⚠️ 需进一步分析"
    
    def _wrap_text(self, text: str, max_line_length: int) -> List[str]:
        """文本换行"""
        words = text.split()
        lines = []
        current_line = []
        
        for word in words:
            if len(' '.join(current_line + [word])) <= max_line_length:
                current_line.append(word)
            else:
                lines.append(' '.join(current_line))
                current_line = [word]
        
        if current_line:
            lines.append(' '.join(current_line))
        
        return lines
    
    def create_tutorial_figure(self, output_path: Optional[Union[str, Path]] = None) -> plt.Figure:
        """
        创建教学图，展示所有可能的基因型组合
        
        Args:
            output_path: 输出文件路径
            
        Returns:
            matplotlib Figure对象
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=self.dpi)
        axes = axes.flatten()
        
        # 定义测试用例
        test_cases = [
            ("AA", "BB", "无重组 - 双纯合亲本型1"),
            ("aa", "bb", "无重组 - 双纯合亲本型2"),
            ("AA", "Bb", "检测到重组 - 标记B杂合"),
            ("aa", "Bb", "检测到重组 - 标记B杂合"),
            ("Aa", "BB", "检测到重组 - 标记A杂合"),
            ("Aa", "Bb", "双杂合 - 需F3验证"),
        ]
        
        for idx, (geno_a, geno_b, title) in enumerate(test_cases):
            if idx < len(axes):
                ax_gel = axes[idx]
                ax_text = None
                
                # 简化绘制，只画凝胶部分
                self._draw_simple_gel(
                    ax_gel, geno_a, geno_b, 
                    product_size_a=350, product_size_b=550
                )
                ax_gel.set_title(title, fontsize=10, fontweight='bold')
        
        plt.suptitle('F2:3群体Indel标记检测 - 所有基因型组合示例', 
                    fontsize=14, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(output_path), dpi=self.dpi, bbox_inches='tight')
            logger.info(f"教学图已保存至: {output_path}")
        
        return fig
    
    def _draw_simple_gel(self, ax, geno_a, geno_b, product_size_a, product_size_b):
        """简化的凝胶绘制（用于教学图）"""
        ax.set_facecolor(self.colors['gel_background'])
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(0, 100)
        ax.invert_yaxis()
        ax.set_axis_off()
        
        # 绘制两个标记的条带
        self._draw_single_band(ax, 0.3, product_size_a/10, f"A: {geno_a}")
        self._draw_single_band(ax, 0.7, product_size_b/10, f"B: {geno_b}")
        
        # 如果是杂合基因型，绘制双带
        def is_heterozygous(genotype):
            if not genotype or len(genotype) < 2:
                return False
            return genotype[0] != genotype[1] or genotype.islower()
        
        if is_heterozygous(geno_a):
            self._draw_heterozygous_bands(ax, 0.3, product_size_a/10, product_size_a, "")
        if is_heterozygous(geno_b):
            self._draw_heterozygous_bands(ax, 0.7, product_size_b/10, product_size_b, "")
    
    def create_primer_visualizations(self, primer_id: str, product_sizes: List[int], 
                                   output_dir: Union[str, Path]) -> Dict[str, str]:
        """
        为引物对创建所有场景的可视化
        
        Args:
            primer_id: 引物对ID
            product_sizes: 产物大小列表 [size_a, size_b]
            output_dir: 输出目录
            
        Returns:
            Dict[str, str]: 可视化文件路径字典
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        # 为每种场景创建图像
        for scenario, (geno_a, geno_b) in self.genotype_examples.items():
            output_file = output_dir / f"{primer_id}_{scenario}.png"
            
            fig = self.visualize_recombination_test(
                marker_a_genotype=geno_a,
                marker_b_genotype=geno_b,
                marker_a_product_size=product_sizes[0],
                marker_b_product_size=product_sizes[1],
                sample_name=primer_id,
                output_path=str(output_file)
            )
            
            plt.close(fig)  # 关闭图形以释放内存
            results[scenario] = str(output_file)
        
        # 创建教学图
        tutorial_file = output_dir / f"{primer_id}_tutorial.png"
        tutorial_fig = self.create_tutorial_figure(str(tutorial_file))
        plt.close(tutorial_fig)
        results['tutorial'] = str(tutorial_file)
        
        logger.info(f"为引物对 {primer_id} 创建了 {len(results)} 个可视化文件")
        return results


# 便捷函数
def visualize_recombination_scenario(
    marker_a_genotype: str,
    marker_b_genotype: str,
    size_a: int = 350,
    size_b: int = 550,
    output_file: Optional[Union[str, Path]] = None,
    sample_name: str = "F2 Plant",
    config: Optional[Dict] = None
) -> plt.Figure:
    """
    快速可视化重组场景的便捷函数
    
    Args:
        marker_a_genotype: 标记A基因型
        marker_b_genotype: 标记B基因型
        size_a: 标记A产物大小
        size_b: 标记B产物大小
        output_file: 输出文件路径
        sample_name: 样本名称
        config: 配置字典
        
    Returns:
        matplotlib Figure对象
    """
    visualizer = GelElectrophoresisVisualizer(config)
    
    fig = visualizer.visualize_recombination_test(
        marker_a_genotype=marker_a_genotype,
        marker_b_genotype=marker_b_genotype,
        marker_a_product_size=size_a,
        marker_b_product_size=size_b,
        sample_name=sample_name,
        output_path=output_file
    )
    
    return fig


def create_tutorial_plot(output_file: Optional[Union[str, Path]] = None, 
                        config: Optional[Dict] = None) -> plt.Figure:
    """
    创建包含所有可能情况的教学图
    
    Args:
        output_file: 输出文件路径
        config: 配置字典
        
    Returns:
        matplotlib Figure对象
    """
    visualizer = GelElectrophoresisVisualizer(config)
    
    fig = visualizer.create_tutorial_figure(output_file)
    return fig
