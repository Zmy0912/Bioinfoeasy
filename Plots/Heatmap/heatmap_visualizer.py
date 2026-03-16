"""
专业热图数据可视化工具
支持CSV/Excel数据读取、聚类分析和多种配色方案
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from pathlib import Path
import sys


class HeatmapVisualizer:
    """热图可视化类"""
    
    def __init__(self, data_file=None, title="热图分析", 
                 xlabel="列", ylabel="行", cmap="viridis",
                 dpi=300, cluster_rows=True, cluster_cols=True,
                 vmin=None, vmax=None,
                 title_size=16, label_size=12, tick_size=8,
                 cbar_width_ratio=0.15, cbar_height_ratio=None):
        """
        初始化热图可视化器
        
        参数:
            data_file: 数据文件路径 (CSV或Excel)
            title: 图表标题
            xlabel: X轴标签
            ylabel: Y轴标签
            cmap: 颜色方案 ('viridis', 'coolwarm', 'plasma', 'inferno', 
                  'RdYlBu', 'RdBu', 'PiYG', 'PRGn' 等)
            dpi: 输出图片分辨率
            cluster_rows: 是否对行进行聚类
            cluster_cols: 是否对列进行聚类
            vmin: 颜色标尺最小值（None表示自动）
            vmax: 颜色标尺最大值（None表示自动）
            title_size: 标题字体大小
            label_size: 轴标签字体大小
            tick_size: 刻度标签字体大小
            cbar_width_ratio: 颜色标尺宽度占主图的比例 (0.05-0.3)
            cbar_height_ratio: 颜色标尺高度占主图的比例 (None表示自适应)
        """
        self.data_file = data_file
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cmap = cmap
        self.dpi = dpi
        self.cluster_rows = cluster_rows
        self.cluster_cols = cluster_cols
        self.vmin = vmin
        self.vmax = vmax
        self.title_size = title_size
        self.label_size = label_size
        self.tick_size = tick_size
        self.cbar_width_ratio = cbar_width_ratio
        self.cbar_height_ratio = cbar_height_ratio
        self.data = None
        self.fig = None
        self.ax = None
        
        # 设置中文字体支持
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
        plt.rcParams['axes.unicode_minus'] = False
    
    def load_data(self, data_file=None):
        """
        加载数据文件
        
        参数:
            data_file: 数据文件路径 (CSV或Excel)
        """
        if data_file:
            self.data_file = data_file
        
        if not self.data_file:
            raise ValueError("请提供数据文件路径")
        
        file_path = Path(self.data_file)
        
        if not file_path.exists():
            raise FileNotFoundError(f"文件不存在: {self.data_file}")
        
        # 根据文件扩展名选择读取方式
        if file_path.suffix.lower() == '.csv':
            # 尝试多种编码格式
            encodings = ['utf-8', 'gbk', 'gb18030', 'gb2312', 'utf-16', 'latin1']
            
            for encoding in encodings:
                try:
                    # 先尝试读取，忽略可能的格式问题
                    try:
                        self.data = pd.read_csv(self.data_file, index_col=0, encoding=encoding)
                        print(f"使用编码 {encoding} 成功读取文件")
                        break
                    except pd.errors.ParserError as e:
                        # 尝试其他分隔符
                        print(f"编码 {encoding} 解析失败，尝试其他分隔符...")
                        for sep in [',', ';', '\t', '|']:
                            try:
                                self.data = pd.read_csv(
                                    self.data_file, 
                                    index_col=0, 
                                    encoding=encoding,
                                    sep=sep,
                                    on_bad_lines='warn'  # 警告但不跳过坏行
                                )
                                print(f"使用编码 {encoding} 和分隔符 '{sep}' 成功读取文件")
                                break
                            except (pd.errors.ParserError, UnicodeDecodeError, UnicodeError):
                                continue
                        else:
                            raise ValueError(
                                f"CSV 文件格式错误：\n"
                                f"- 文件可能包含不一致的列数\n"
                                f"- 检查第 4 行附近的数据格式\n"
                                f"- 确保所有行的列数相同\n"
                                f"- 建议使用 Excel 打开并另存为 .xlsx 格式"
                            )
                except (UnicodeDecodeError, UnicodeError):
                    continue
            else:
                raise ValueError(f"无法识别文件编码，请确保文件使用 UTF-8、GBK 或 GB18030 编码")
        elif file_path.suffix.lower() in ['.xlsx', '.xls']:
            self.data = pd.read_excel(self.data_file, index_col=0)
        else:
            raise ValueError("不支持的文件格式，请使用CSV或Excel文件")
        
        # 确保数据为数值类型
        self.data = self.data.astype(float)
        
        print(f"数据加载成功: {self.data.shape[0]} 行 × {self.data.shape[1]} 列")
        return self.data
    
    def generate_sample_data(self, n_rows=15, n_cols=12, seed=42):
        """
        生成示例数据（用于测试）
        
        参数:
            n_rows: 行数
            n_cols: 列数
            seed: 随机种子
        """
        np.random.seed(seed)
        
        # 生成有模式的示例数据
        base_pattern = np.linspace(0, 10, n_cols)
        self.data = pd.DataFrame(
            np.random.randn(n_rows, n_cols) * 0.5 + np.tile(base_pattern, (n_rows, 1)),
            index=[f'样本_{i+1}' for i in range(n_rows)],
            columns=[f'基因_{j+1}' for j in range(n_cols)]
        )
        
        print(f"示例数据生成成功: {self.data.shape[0]} 行 × {self.data.shape[1]} 列")
        return self.data
    
    def _cluster_data(self):
        """
        执行聚类分析并返回排序后的数据
        """
        clustered_data = self.data.copy()
        
        # 行聚类
        if self.cluster_rows:
            row_linkage = linkage(clustered_data.values, method='average')
            row_order = dendrogram(row_linkage, no_plot=True)['leaves']
            clustered_data = clustered_data.iloc[row_order, :]
        
        # 列聚类
        if self.cluster_cols:
            col_linkage = linkage(clustered_data.values.T, method='average')
            col_order = dendrogram(col_linkage, no_plot=True)['leaves']
            clustered_data = clustered_data.iloc[:, col_order]
        
        return clustered_data
    
    def plot_heatmap(self, save_path=None):
        """
        绘制热图
        
        参数:
            save_path: 保存路径（PNG格式），如果不指定则不保存
        """
        if self.data is None:
            raise ValueError("请先加载数据或生成示例数据")
        
        # 执行聚类
        clustered_data = self._cluster_data()
        
        # 计算图形尺寸
        n_rows, n_cols = clustered_data.shape
        fig_width = max(8, n_cols * 0.6 + 2)
        fig_height = max(6, n_rows * 0.6 + 1.5)
        
        # 创建图形
        self.fig = plt.figure(figsize=(fig_width, fig_height))
        
        # 计算网格布局（颜色标尺单独放在右侧空白处）
        # 根据用户设置的颜色标尺宽度比例计算布局
        main_area_ratio = 1 - self.cbar_width_ratio
        main_area_ratio = max(0.7, min(0.95, main_area_ratio))  # 限制在 0.7-0.95 范围内
        
        if self.cluster_rows and self.cluster_cols:
            # 主图形区域，颜色标尺占 cbar_width_ratio
            gs = self.fig.add_gridspec(
                nrows=2, ncols=2, 
                width_ratios=[main_area_ratio, self.cbar_width_ratio],
                height_ratios=[0.15, 0.85],
                wspace=0.02, hspace=0.02,
                left=0.02, right=0.95, top=0.95, bottom=0.02
            )
            ax_dendro_col = self.fig.add_subplot(gs[0, 0])
            ax_dendro_row = self.fig.add_subplot(gs[1, 0])
            ax_heatmap = self.fig.add_subplot(gs[1, 0])
            cax = self.fig.add_subplot(gs[0:2, 1])
        elif self.cluster_rows:
            gs = self.fig.add_gridspec(
                nrows=1, ncols=2,
                width_ratios=[main_area_ratio, self.cbar_width_ratio],
                wspace=0.02,
                left=0.02, right=0.95, top=0.95, bottom=0.05
            )
            ax_dendro_row = self.fig.add_subplot(gs[0, 0])
            ax_heatmap = self.fig.add_subplot(gs[0, 0])
            cax = self.fig.add_subplot(gs[0, 1])
            ax_dendro_col = None
        elif self.cluster_cols:
            gs = self.fig.add_gridspec(
                nrows=2, ncols=2,
                height_ratios=[0.15, 0.85],
                width_ratios=[main_area_ratio, self.cbar_width_ratio],
                hspace=0.02, wspace=0.02,
                left=0.02, right=0.95, top=0.95, bottom=0.02
            )
            ax_dendro_col = self.fig.add_subplot(gs[0, 0])
            ax_heatmap = self.fig.add_subplot(gs[1, 0])
            cax = self.fig.add_subplot(gs[0:2, 1])
            ax_dendro_row = None
        else:
            # 无聚类时，使用更简单的布局
            gs = self.fig.add_gridspec(
                nrows=1, ncols=2,
                width_ratios=[main_area_ratio, self.cbar_width_ratio],
                wspace=0.02,
                left=0.05, right=0.95, top=0.90, bottom=0.1
            )
            ax_heatmap = self.fig.add_subplot(gs[0, 0])
            cax = self.fig.add_subplot(gs[0, 1])
            ax_dendro_row = None
            ax_dendro_col = None
        
        # 绘制列树状图
        if self.cluster_cols and ax_dendro_col is not None:
            col_linkage = linkage(clustered_data.values.T, method='average')
            dendrogram(col_linkage, ax=ax_dendro_col, orientation='top')
            ax_dendro_col.set_xticks([])
            ax_dendro_col.set_yticks([])
            ax_dendro_col.spines['top'].set_visible(False)
            ax_dendro_col.spines['right'].set_visible(False)
            ax_dendro_col.spines['bottom'].set_visible(False)
            ax_dendro_col.spines['left'].set_visible(False)
        else:
            ax_dendro_col = None
        
        # 绘制行树状图
        if self.cluster_rows and ax_dendro_row is not None:
            row_linkage = linkage(clustered_data.values, method='average')
            dendrogram(row_linkage, ax=ax_dendro_row, orientation='left')
            ax_dendro_row.set_xticks([])
            ax_dendro_row.set_yticks([])
            ax_dendro_row.spines['top'].set_visible(False)
            ax_dendro_row.spines['right'].set_visible(False)
            ax_dendro_row.spines['bottom'].set_visible(False)
            ax_dendro_row.spines['left'].set_visible(False)
        else:
            ax_dendro_row = None
        
        # 绘制热图
        sns.heatmap(
            clustered_data,
            ax=ax_heatmap,
            cmap=self.cmap,
            cbar_ax=cax,
            annot=False,  # 不显示单元格数值
            fmt='.2f',
            linewidths=0.5,
            square=True,
            xticklabels=True,
            yticklabels=True,
            vmin=self.vmin,
            vmax=self.vmax
        )
        
        # 调整颜色标尺的高度（如果指定了高度比例）
        if self.cbar_height_ratio is not None and cax is not None:
            # 获取当前颜色标尺的位置
            cbox = cax.get_position()
            
            # 计算新的位置
            new_bottom = cbox.y0 + (1 - self.cbar_height_ratio) / 2 * cbox.height
            new_height = cbox.height * self.cbar_height_ratio
            
            # 更新颜色标尺位置
            cax.set_position([cbox.x0, new_bottom, cbox.width, new_height])
        
        # 设置标签（使用自定义字体大小）
        ax_heatmap.set_title(self.title, fontsize=self.title_size, pad=20)
        ax_heatmap.set_xlabel(self.xlabel, fontsize=self.label_size)
        ax_heatmap.set_ylabel(self.ylabel, fontsize=self.label_size)
        
        # 调整刻度标签
        ax_heatmap.tick_params(axis='x', rotation=45, labelsize=self.tick_size)
        ax_heatmap.tick_params(axis='y', rotation=0, labelsize=self.tick_size)
        
        # 自动调整布局
        if not (self.cluster_rows and self.cluster_cols):
            plt.tight_layout()
        
        # 保存图片
        if save_path:
            self.save_heatmap(save_path)
        
        return self.fig, ax_heatmap
    
    def save_heatmap(self, save_path):
        """
        保存热图为PNG格式
        
        参数:
            save_path: 保存路径
        """
        if self.fig is None:
            raise ValueError("请先绘制热图")
        
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        
        if save_path.suffix.lower() != '.png':
            save_path = save_path.with_suffix('.png')
        
        self.fig.savefig(
            save_path,
            dpi=self.dpi,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none'
        )
        print(f"热图已保存至: {save_path}")
    
    def show(self):
        """显示热图"""
        if self.fig is not None:
            plt.show()
        else:
            print("请先绘制热图")


def main():
    """主函数示例"""
    
    print("=" * 60)
    print("专业热图数据可视化工具")
    print("=" * 60)
    
    # 示例1: 使用示例数据生成热图
    print("\n示例1: 生成示例数据并绘制热图 (默认设置)")
    print("-" * 60)
    
    visualizer = HeatmapVisualizer(
        title="基因表达热图分析 (viridis 默认)",
        xlabel="基因",
        ylabel="样本",
        cmap="viridis",
        dpi=300,
        cluster_rows=True,
        cluster_cols=True
    )
    
    # 生成示例数据
    visualizer.generate_sample_data(n_rows=20, n_cols=15)
    
    # 绘制并保存热图
    fig, ax = visualizer.plot_heatmap(save_path="heatmap_example_viridis.png")
    
    # 示例2: 使用不同的配色方案 + 大字体 + 宽颜色标尺
    print("\n示例2: 使用 coolwarm 配色 + 大字体 + 宽颜色标尺")
    print("-" * 60)
    
    visualizer2 = HeatmapVisualizer(
        title="基因表达热图分析 (coolwarm 大字体 宽标尺)",
        xlabel="基因",
        ylabel="样本",
        cmap="coolwarm",
        vmin=0,
        vmax=10,
        title_size=20,        # 大标题
        label_size=14,        # 大轴标签
        tick_size=10,         # 大刻度
        cbar_width_ratio=0.25,# 宽颜色标尺 (25%)
        cbar_height_ratio=0.9,# 高颜色标尺 (90%)
        dpi=300,
        cluster_rows=True,
        cluster_cols=True
    )
    
    visualizer2.generate_sample_data(n_rows=20, n_cols=15)
    visualizer2.plot_heatmap(save_path="heatmap_example_coolwarm.png")
    
    # 示例3: 使用 plasma 配色 + 小字体 + 窄颜色标尺
    print("\n示例3: 使用 plasma 配色 + 小字体 + 窄颜色标尺")
    print("-" * 60)
    
    visualizer3 = HeatmapVisualizer(
        title="基因表达热图分析 (plasma 小字体 窄标尺)",
        xlabel="基因",
        ylabel="样本",
        cmap="plasma",
        title_size=12,        # 小标题
        label_size=10,        # 小轴标签
        tick_size=6,          # 小刻度
        cbar_width_ratio=0.08,# 窄颜色标尺 (8%)
        cbar_height_ratio=0.6,# 短颜色标尺 (60%)
        dpi=300,
        cluster_rows=True,
        cluster_cols=True
    )
    
    visualizer3.generate_sample_data(n_rows=20, n_cols=15)
    visualizer3.plot_heatmap(save_path="heatmap_example_plasma.png")
    
    # 示例2: 使用自己的数据文件
    print("\n示例2: 使用CSV/Excel数据文件")
    print("-" * 60)
    print("如需使用自己的数据，请使用以下代码:")
    print("""
# 方法1: 使用viridis配色
visualizer = HeatmapVisualizer(
    data_file="your_data.csv",  # 或 "your_data.xlsx"
    title="自定义标题",
    xlabel="列标签",
    ylabel="行标签",
    cmap="viridis",
    cluster_rows=True,
    cluster_cols=True
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="custom_heatmap.png")
visualizer.show()

# 方法2: 使用蓝-白-红配色
visualizer = HeatmapVisualizer(
    data_file="your_data.csv",
    title="自定义标题",
    xlabel="列标签",
    ylabel="行标签",
    cmap="coolwarm",  # 蓝-白-红配色
    cluster_rows=True,
    cluster_cols=True
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="custom_heatmap.png")
visualizer.show()

# 方法3: 自定义颜色范围、字体大小和颜色标尺
visualizer = HeatmapVisualizer(
    data_file="your_data.csv",
    title="自定义标题",
    xlabel="列标签",
    ylabel="行标签",
    cmap="coolwarm",
    vmin=0,               # 颜色标尺最小值
    vmax=10,              # 颜色标尺最大值
    title_size=20,        # 标题字体大小
    label_size=14,        # 轴标签字体大小
    tick_size=10,         # 刻度标签字体大小
    cbar_width_ratio=0.2, # 颜色标尺宽度占主图比例 (0.05-0.3)
    cbar_height_ratio=0.8,# 颜色标尺高度占主图比例 (0.3-1.0, None表示自适应)
    cluster_rows=True,
    cluster_cols=True
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="custom_heatmap.png")
visualizer.show()

# 方法4: 不进行聚类
visualizer = HeatmapVisualizer(
    data_file="your_data.csv",
    title="自定义标题",
    xlabel="列标签",
    ylabel="行标签",
    cmap="viridis",
    cluster_rows=False,
    cluster_cols=False
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="custom_heatmap.png")
visualizer.show()

# 常用配色方案:
# 'viridis'   - 紫-蓝-绿-黄 (默认)
# 'plasma'    - 蓝-红-黄
# 'inferno'   - 黑-红-黄
# 'coolwarm'  - 蓝-白-红 (适合相关性)
# 'RdBu'      - 红-白-蓝
# 'RdYlBu'    - 红-黄-蓝
# 'PiYG'      - 紫-白-绿
# 'PRGn'      - 紫-白-绿
# 'Greens'    - 绿色渐变
# 'Blues'     - 蓝色渐变
# 'Reds'      - 红色渐变

# 字体大小参考:
# title_size  - 标题大小，建议 14-24
# label_size  - 轴标签大小，建议 10-16
# tick_size   - 刻度标签大小，建议 6-12

# 颜色标尺调整:
# cbar_width_ratio   - 颜色标尺宽度占主图比例 (建议 0.05-0.3)
# cbar_height_ratio  - 颜色标尺高度占主图比例 (建议 0.3-1.0, None=自适应)
    """)
    
    print("\n" + "=" * 60)
    print("程序执行完成！")
    print("生成的热图已保存为: heatmap_example.png")
    print("=" * 60)


if __name__ == "__main__":
    main()
