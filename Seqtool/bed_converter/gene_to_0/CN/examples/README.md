# 示例文件说明

本文件夹包含用于测试 BED/GFF3 坐标转换工具的示例文件。

## 文件列表

### 1. example.bed
- BED 格式示例文件
- 包含 3 个基因（AoNCED1, AoNCED2, AoNCED3）
- 包含 gene, mRNA, exon, CDS 等特征
- 包含相位信息
- 带有中文注释

### 2. example.gff3
- 标准 GFF3 格式示例文件
- 包含 3 个基因（AoNCED1, AoNCED2, AoNCED3）
- 位于不同的染色体（chr1, chr2）
- 包含完整的 GFF3 头信息
- 包含 gene, mRNA, exon, CDS 特征
- 包含正负链基因
- 属性字段包含 gene_id

## 使用示例

### 测试 BED 格式转换

1. 启动 `bed_converter.py`
2. 选择 "BED格式" 单选按钮
3. 浏览并选择 `examples/example.bed`
4. 点击 "加载文件" 查看原始内容
5. 点击 "转换" 进行坐标转换
6. 查看预览区域的结果
7. 检查统计信息（基因数和特征数）
8. 点击 "保存文件" 保存结果

### 测试 GFF3 格式转换

1. 启动 `bed_converter.py`
2. 选择 "GFF3格式" 单选按钮
3. 浏览并选择 `examples/example.gff3`
4. 点击 "加载文件" 查看原始内容
5. 点击 "转换" 进行坐标转换
6. 查看预览区域的结果
7. 检查统计信息（基因数和特征数）
8. 注意 GFF3 格式保持 1-based 坐标系统
9. 点击 "保存文件" 保存结果

## 预期结果

### BED 格式转换后

每个基因的第一个特征应该从坐标 0 开始。例如 AoNCED1:

```bed
AoNCED1	0	2958	gene
AoNCED1	2444	2958	mRNA
AoNCED1	2444	2958	exon
AoNCED1	992	995	exon
```

### GFF3 格式转换后

每个基因的第一个特征应该从坐标 1 开始（GFF3 标准）。例如 AoNCED1:

```gff3
chr1	Ensembl	gene	1	2958	.	+	.	ID=gene:AO001;gene_id=AoNCED1;Name=NCED1
chr1	Ensembl	mRNA	1	2958	.	+	.	ID=transcript:AO001.1;Parent=gene:AO001;gene_id=AoNCED1
chr1	Ensembl	exon	2445	2958	.	+	.	ID=exon:AO001.1.1;Parent=transcript:AO001.1;gene_id=AoNCED1
chr1	Ensembl	exon	993	996	.	+	.	ID=exon:AO001.1.2;Parent=transcript:AO001.1;gene_id=AoNCED1
```

## 验证要点

1. **坐标范围**: 检查每个基因内的坐标是否从 0（BED）或 1（GFF3）开始
2. **相对距离**: 验证特征之间的相对距离是否保持不变
3. **特征数量**: 确认转换后的特征数量与原始文件一致
4. **属性保留**: GFF3 文件的所有属性应完整保留
5. **注释行**: 注释行应完整保留
6. **统计信息**: 基因数和特征数应正确显示

## 格式验证

### BED 格式验证
- 至少包含 4 列
- 坐标为非负整数
- 基因 ID 相同的特征归为一组

### GFF3 格式验证
- 包含 9 列
- 起始位置 < 结束位置
- attributes 字段包含 gene_id
- 坐标为正整数（1-based）

## 故障排除

如果转换失败或结果不正确：

1. 检查文件格式是否符合要求
2. 确认选择了正确的文件类型（BED/GFF3）
3. 查看 BED 示例文件的格式对比
4. 检查 GFF3 文件中是否包含 gene_id 属性
5. 参考 FAQ 文档中的常见问题

## 自定义示例

您可以基于这些示例文件创建自己的测试文件：

- 修改基因 ID
- 调整坐标范围
- 添加更多特征类型
- 更新属性字段

确保文件格式与示例文件保持一致。
