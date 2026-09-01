# 本地 scRNA 数据清单

本地数据目录：`E:\nhm\scRNA`

## 顶层内容

- `h5/`：多个样本的 10x Cell Ranger 输出，包含 `filtered_feature_bc_matrix.h5`、`raw_feature_bc_matrix.h5`、`molecule_info.h5`
- `placenta/`：胎盘分析目录，包含 Cell Ranger 输出、处理后的 h5ad 对象和结果
- `Human-Maternal-Fetal-Interface_MFI-main/`：参考或复用的母胎界面分析项目代码与 notebook
- `scRNA_PY1.ipynb`：当前重点分析 notebook，约 12.9MB
- `scRNA_py.ipynb`、`jihetuxing.ipynb`：其他分析 notebook
- `main.zip`：约 3.04GB，不建议放入普通 GitHub 仓库

## 文件类型概览

- 无扩展名 Cell Ranger 工作文件：约 2040 个
- `.csv`：约 437 个
- `.pdf`：约 331 个
- `.ipynb`：约 214 个
- `.h5`：约 47 个
- `.txt`：约 39 个
- `.pkl`：约 24 个
- `.gz`：约 15 个
- `.h5ad`：2 个
- `.bam`：2 个
- `.zip`：1 个

## 大文件

这些文件不适合直接提交到普通 GitHub 仓库：

- `placenta/data/cellranger/Pla_HDBR10142767/outs/possorted_genome_bam.bam`：约 43GB
- `placenta/processed_adata.h5ad`：约 3.64GB
- `main.zip`：约 3.04GB
- `placenta/combined_processed.h5ad`：约 1.12GB
- 多个 `molecule_info.h5`：约 54MB 到 1.01GB
- 多个 `raw_feature_bc_matrix.h5`：约 88MB 到 146MB
- 多个大型 notebook：约 50MB 到 103MB

## GitHub 管理建议

建议 GitHub 仓库保存：分析脚本、notebook 或提取脚本、配置文件、README、任务清单、小型结果表格或代表性图片。

建议不要保存：`.bam`、`.h5ad`、`molecule_info.h5`、大型 `.h5`、Cell Ranger 临时工作目录和大型压缩包。

大型数据可用 Git LFS、DVC、机构服务器、NAS、OneDrive/Dropbox 或对象存储管理。
