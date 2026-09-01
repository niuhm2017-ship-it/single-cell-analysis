# scRNA_PY1.ipynb 分析结果解读

来源 notebook：`E:\nhm\scRNA\scRNA_PY1.ipynb`

## Notebook 概况

- Kernel：Python，notebook 元数据为 Python 3.12.8；历史输出显示实际运行环境曾为 Python 3.8 的 `scRNA` conda 环境。
- 单元数：61
- 代码单元：47
- Markdown 单元：14
- 输出单元：47
- 图像输出：18 个

## 分析对象

该 notebook 使用 Scanpy/scverse 分析胎盘单细胞 RNA-seq 数据。代码从 `/mnt/e/nhm/scRNA/h5` 读取多个 10x Genomics `filtered_feature_bc_matrix.h5` 文件，并合并为一个 AnnData 对象。

样本映射中可见 14 个样本：

- `Pla_HDBR10917733`：`age_5week_Hrv100_1`
- `Pla_HDBR10701666`：`age_6week_Hrv43_1`
- `Pla_HDBR10142767`：`age_6week_Hrv43_2`
- `Pla_HDBR10142768`：`age_6week_Hrv43_3`
- `Pla_HDBR10701667`：`age_6week_Hrv43_4`
- `Pla_HDBR10917731`：`age_8week_Hrv99_1`
- `Pla_HDBR10142769`：`age_9week_Hrv46_1`
- `Pla_HDBR10142770`：`age_9week_Hrv46_2`
- `Pla_HDBR10701668`：`age_9week_Hrv46_3`
- `Pla_HDBR10917730`：`age_9week_Hrv98_1`
- `Pla_HDBR8624430`：`age_12week_H2_1`
- `Pla_HDBR8624431`：`age_12week_H2_2`
- `Pla_HDBR8715512`：`age_12_13week_H7a9_1`
- `Pla_HDBR8715514`：`age_12_13week_H7a9_2`

## 已完成的分析步骤

1. 导入核心库：`scanpy`、`anndata`、`numpy`、`pandas`、`scipy`、`scrublet` 等。
2. 批量读取 10x h5 文件并合并 AnnData。
3. 修复重复基因名：`adata.var_names_make_unique()`。
4. 添加样本注释：`sample` 和更易读的 `name`。
5. 检查关键 marker：`GAPDH`、`ACTB`、`CD4`、`CD8A`、`MKI67`、`PRDM1` 均存在。
6. 计算 QC 指标：线粒体基因、核糖体基因、血红蛋白基因比例。
7. 绘制 QC 小提琴图和散点图。
8. 过滤细胞和基因：`min_genes=100`、`min_cells=3`。
9. 使用 Scrublet 预测双细胞，并写入 `doublet_score`、`predicted_doublet`。
10. 保存原始 count 到 `adata.layers["counts"]`。
11. 归一化并 log 转换。
12. 选择 2000 个高变基因。
13. 运行 PCA、邻接图、UMAP。
14. 运行 Leiden 聚类，包括默认聚类和 `0.02`、`0.50`、`2.00` 多分辨率聚类。
15. 使用 marker gene dotplot 辅助细胞类型注释。
16. 进行 cluster 差异表达分析。
17. 保存整合后的 h5ad：`/mnt/storage/yxh/nhm/scRNA/placenta/combined_processed.h5ad`。

## 主要结果

### 样本整合

notebook 将 14 个不同孕周或样本来源的胎盘样本整合到同一个 AnnData 对象中，并通过 UMAP 按 `sample` 和 `name` 可视化。已有输出显示样本 ID 已成功映射为更易读的孕周/个体名称。

### 细胞类型初步注释

低分辨率聚类 `leiden_res_0.02` 被映射为三个一级细胞类型：

- cluster `0`：`B Cells`
- cluster `1`：`Monocytes`
- cluster `2`：`Erythroid`

但最后的计数结果显示映射后的标签为：

| 细胞类型 | 细胞数 |
| --- | ---: |
| Lymphocytes | 108,024 |
| Monocytes | 37,834 |
| Erythroid | 11,952 |

总计约 157,810 个细胞。说明 notebook 后续可能把 `B Cells` 合并或重命名为 `Lymphocytes`，最终结果以输出中的 `cell_type_lvl1` 为准。

### cluster 7 的差异表达结果

针对 `leiden_res_0.50` 的 cluster 7，前 5 个差异 marker 为：

| gene | score | logfoldchange | adjusted p |
| --- | ---: | ---: | ---: |
| CD247 | 111.58 | 4.60 | 0.0 |
| PTPRC | 107.71 | 3.15 | 0.0 |
| SRGN | 95.69 | 2.55 | 0.0 |
| CORO1A | 94.70 | 3.80 | 0.0 |
| NKG7 | 93.68 | 4.14 | 0.0 |

这些 marker 中 `CD247`、`PTPRC`、`NKG7` 指向免疫细胞，尤其是 T/NK 或广义淋巴细胞特征。cluster 7 很可能属于 lymphoid lineage，需要结合 dotplot 和 UMAP feature plot 进一步确认是 T cell、NK cell，还是混合状态。

## 需要注意的问题

### 1. 高变基因选择时机有问题

notebook 在归一化和 log 转换之后运行 `flavor="seurat_v3"` 的高变基因选择。输出中出现警告：`flavor='seurat_v3' expects raw count data, but non-integers were found.`

这说明 `seurat_v3` 方法期望输入原始 counts，但当前 `adata.X` 已经是归一化/log 后的数据。建议在 log 转换前基于原始 counts 做 HVG，或使用适合 log 数据的 HVG flavor。

### 2. 双细胞检测需要明确记录

Scrublet 已经生成 `doublet_score` 和 `predicted_doublet`，提取代码中可见后续执行了双细胞过滤。建议在报告中补充过滤前后细胞数和双细胞比例。

### 3. 细胞类型命名不一致

代码里曾将 `leiden_res_0.02` 映射为 `B Cells`、`Monocytes`、`Erythroid`，但最终统计显示 `Lymphocytes`、`Monocytes`、`Erythroid`。建议统一命名并补充注释依据。

### 4. 路径依赖较强

notebook 使用多个绝对路径，例如 `/mnt/e/nhm/scRNA/h5` 和 `/mnt/storage/yxh/nhm/scRNA/placenta/combined_processed.h5ad`。建议改为项目相对路径或配置文件路径，方便复现。

### 5. 大文件不适合直接放入 GitHub

本地目录存在 43GB BAM、3.64GB h5ad、3.04GB zip、多个数百 MB 到 1GB 的 h5 文件。建议放在外部存储、Git LFS、DVC 或对象存储中。

## 建议下一步

1. 把 `scRNA_PY1.ipynb` 拆成可复现脚本：加载合并、QC/Scrublet、归一化聚类、注释差异分析。
2. 把样本映射表写入 `configs/samples.csv`。
3. 把 QC 阈值、HVG 参数、Leiden resolution 写入 `configs/analysis_params.yaml`。
4. 修正 HVG 选择流程，避免在 log 数据上使用 `seurat_v3`。
5. 保存关键图表和差异基因表到 `results/`，并在 README 中记录主要结论。
