# GFF3 to BED 转换器

一个具有图形化界面的GFF3转BED格式转换工具。

## 功能特点

- 图形化用户界面，操作简单
- 支持选择多种特征类型进行转换(gene, mRNA, exon, CDS)
- 实时显示转换日志
- 自动生成BED文件

## 文件说明

- `gff3_to_bed_gui.py` - 主程序文件
- `start_converter.bat` - Windows启动脚本
- `requirements.txt` - 依赖包列表(主要使用Python标准库)

## 使用方法

### 方法一：使用启动脚本(推荐)

直接双击 `start_converter.bat` 文件即可启动程序。

### 方法二：使用Python命令

打开命令行,进入程序目录,运行:

```bash
python gff3_to_bed_gui.py
```

## 操作步骤

1. 点击"浏览"按钮选择GFF3输入文件
2. 程序会自动解析文件并显示检测到的基因信息
3. 程序会自动设置BED输出文件路径(也可手动修改)
4. 勾选需要转换的特征类型
5. (可选)点击"预览输出"查看转换结果
6. 点击"开始转换"或从预览窗口保存BED文件

## 主要功能

- **基因名称一致性保证**: 确保同一基因的所有特征使用相同的基因名称
- **智能基因识别**: 自动从GFF3的attributes字段中提取基因ID
- **预览功能**: 转换前可预览输出内容,避免保存不符合预期的结果
- **实时日志**: 显示转换进度和统计信息

## 格式说明

### GFF3格式

GFF3文件是以制表符分隔的文本文件,包含9列:

1. seqid - 序列标识符
2. source - 来源
3. type - 特征类型(gene, mRNA, exon, CDS等)
4. start - 起始位置(1-based)
5. end - 结束位置
6. score - 分数
7. strand - 链方向(+/-/.)
8. phase - 相位
9. attributes - 属性信息

### BED格式

BED文件是以制表符分隔的文本文件,本程序生成的BED格式包含4-5列:

1. geneID/transcriptID - 基因ID或转录本ID(格式: geneID 或 geneID/transcriptID)
2. start - 起始位置(1-based,保持GFF3原始坐标)
3. end - 结束位置
4. featureType - 特征类型(gene, mRNA, exon, CDS等)
5. phase - 相位(可选列,仅CDS特征有值)

## 注意事项

- BED输出使用1-based坐标系统,保持GFF3原始坐标
- 确保GFF3文件格式正确
- BED输出文件会被覆盖,请谨慎选择输出路径
- gene类型特征只显示geneID
- 其他特征类型(mRNA, exon, CDS等)显示 geneID/transcriptID 格式
- phase列仅对CDS特征有实际值,其他特征为空

## 示例

在程序目录中提供了示例文件:
- `基因成员.gff3` - GFF3格式示例输入
- `BED示例.txt` - BED格式示例输出

## 系统要求

- Python 3.6 或更高版本
- Windows 操作系统
- 无需安装额外依赖包

## 许可证

本工具为免费软件,可自由使用和修改。
