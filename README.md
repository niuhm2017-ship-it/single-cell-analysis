# 单细胞分析管理项目模板

这是一个用于管理单细胞 RNA-seq 分析的标准项目模板，兼容 Seurat 和 Scanpy 两条常见分析路线。建议把原始数据、脚本、配置和结果分开管理，保证每次分析都能追踪来源、参数和输出。

## 项目结构

```text
single-cell-analysis-template/
├── configs/                 # 分析参数、样本信息、marker 基因表
├── data/                    # 数据目录，建议不提交大型原始数据
│   ├── raw/                 # 10x、FASTQ、表达矩阵等原始输入
│   ├── processed/           # 过滤、归一化、整合后的中间数据
│   └── metadata/            # 样本分组、批次、临床信息等
├── notebooks/               # 探索性分析和报告型 notebook
├── results/                 # 图表、表格、最终对象
│   ├── figures/
│   ├── tables/
│   └── objects/
├── scripts/                 # 可重复运行的分析脚本
└── TASKS.md                 # 项目任务清单
```

## 快速开始

1. 把原始数据放入 `data/raw/`，例如 10x Genomics 输出的 `filtered_feature_bc_matrix/`。
2. 在 `configs/samples.csv` 中登记样本名、分组、批次和数据路径。
3. 根据你使用的生态选择流程：
   - R/Seurat：运行 `scripts/seurat_pipeline.R`
   - Python/Scanpy：运行 `scripts/scanpy_pipeline.py`
4. 结果会输出到 `results/figures/`、`results/tables/` 和 `results/objects/`。

## 推荐分析流程

### 1. 数据导入

- 读取 10x 矩阵、h5ad、loom 或自定义表达矩阵。
- 合并样本前保留样本 ID、分组、批次等元数据。
- 记录原始细胞数和基因数，便于后续 QC 对比。

### 2. 质量控制

常用 QC 指标：

- 每个细胞检测到的基因数
- 每个细胞 UMI/count 数
- 线粒体基因比例
- 红细胞或核糖体基因比例
- 双细胞风险评分

建议把阈值写入 `configs/analysis_params.yaml`，不要散落在脚本里。

### 3. 归一化与高变基因

Seurat 常用：

- `NormalizeData`
- `FindVariableFeatures`
- `ScaleData`

Scanpy 常用：

- `sc.pp.normalize_total`
- `sc.pp.log1p`
- `sc.pp.highly_variable_genes`

### 4. 降维、邻接图和聚类

典型步骤：

- PCA
- 构建邻接图
- UMAP/t-SNE
- Leiden/Louvain 聚类
- 按样本、分组、批次检查聚类是否合理

### 5. 批次校正与整合

可选方法：

- Seurat integration
- Harmony
- Scanorama
- BBKNN
- scVI

批次校正前后都建议保留 UMAP 图和关键 marker 表达图。

### 6. 细胞类型注释

推荐组合：

- 经典 marker 基因
- 自动注释工具，如 SingleR、Azimuth、CellTypist
- 文献或公共数据库交叉验证

Marker 表建议维护在 `configs/marker_genes.csv`。

### 7. 差异表达与通路分析

常见比较：

- 不同细胞群之间
- 同一细胞类型内不同处理组之间
- 疾病组 vs 对照组

输出建议：

- 差异基因表
- 火山图
- DotPlot/Heatmap
- GO/KEGG/GSEA 富集结果

## 数据管理建议

- 不要把大型原始数据直接提交到 GitHub。
- 可以提交小型示例数据和配置文件。
- 大型数据建议使用 Git LFS、DVC、对象存储或机构服务器。
- 每次分析都记录软件版本、参数、输入路径和输出文件名。

## 环境建议

R/Seurat:

```r
install.packages("Seurat")
install.packages("tidyverse")
install.packages("patchwork")
```

Python/Scanpy:

```bash
conda create -n sc-analysis python=3.11
conda activate sc-analysis
pip install scanpy anndata pandas numpy matplotlib seaborn pyyaml
```

## 下一步

查看 `TASKS.md`，按项目阶段推进分析和管理工作。
