# MED12–SUPT6H–H2A.Z 分析代码

用于研究 MED12–SUPT6H 轴、H2A.Z 相关染色质组织、转录暂停释放及基因表达的分析脚本。

[English](README.md)

## 项目概述

本项目收录 APEX 质谱、ChIP-seq、PRO-seq、TT-seq、MNase-seq 和小鼠肝脏 RNA-seq 分析代码，涵盖差异分析、基因组信号定量、通路富集以及染色质与转录数据整合。

当前版本为已有脚本的整理集合。部分脚本需要配置本地路径并补充外部输入文件，尚未完成端到端流程验证。

## 目录结构

| 路径 | 内容 |
| --- | --- |
| `archive/MS/` | APEX 质谱分析 |
| `archive/CHIP-SEQ/` | ChIP-seq 注释及区域信号分析 |
| `archive/pro/` | PRO-seq 暂停指数分析 |
| `archive/TT/` | TT-seq 分析及 metagene 输入准备 |
| `archive/MNase/` | 核小体位置分析 |
| `archive/mouse/` | 小鼠 RNA-seq 差异分析 |
| `archive/mRNAseq/` | H2A.Z 与 RNA-seq 整合 |
| `archive/med12_disease/` | 外部数据集比较 |
| `archive/ratio/` | Metagene 可视化工具 |
| `notebooks/` | 已清除保存输出的 notebook |
| `scripts/` | Notebook 的纯代码导出 |
| `docs/` | 输入要求、已知问题及验证记录 |
| `SOURCE_MANIFEST.tsv` | 原始来源及 SHA-256 校验值 |

`archive/` 保留原始脚本，具体来源见校验清单。

## 环境要求

- Python 及 `requirements-analysis.txt` 中列出的依赖。
- 对应模块需要的 R 包，如 DESeq2、tidyverse、clusterProfiler 和 ChIPseeker。
- Shell 流程所需的 Linux 或 WSL 环境及测序分析工具。
- 与分析匹配的参考基因组注释、输入表格及信号文件。
- 部分绘图脚本要求安装 Arial 字体。

依赖列表尚未作为完整环境进行验证。详细要求见[输入与运行说明](docs/INPUTS_AND_EXECUTION.md)。

## 快速开始

```bash
git clone https://github.com/Berialim/MED12-SUPT6H-H2AZ-analysis.git
cd MED12-SUPT6H-H2AZ-analysis
python -m venv .venv
# 按当前操作系统和 shell 的方式激活 .venv。
python -m pip install -r requirements-analysis.txt
```

选择所需模块，准备输入文件并检查脚本中的路径与参数。R 包及外部命令行工具需要按模块另行安装。

## 使用示例

```bash
python archive/TT/plot_h2az_boxplot.py /path/to/results.txt \
  --density /path/to/density.bed -o /path/to/tt_h2az.pdf

python archive/pro/plot_pauseindex_h2az_boxplot.py /path/to/pauseindex.tsv \
  --density /path/to/density.bed \
  --ctrl-pause CTRL_pause --kd-pause KD_pause \
  -o /path/to/pro_h2az.pdf
```

上述命令使用已准备好的结果进行可视化，不包含原始测序数据处理。密度文件默认第 4 列为基因标识、第 7 列为 H2A.Z 信号；PRO-seq 输入表需包含 `gene` 列和指定的暂停指数列。

其他归档脚本可能包含固定路径或工作目录假设，运行前需逐项检查。数据处理请使用输入文件的工作副本。

## 数据获取

相关稿件列出的 GEO 编号为 **GSE325384、GSE325382、GSE325378 和 GSE325610**。本项目尚未独立核实这些数据的公开状态及其与测序类型的对应关系。

仓库不包含原始测序数据、参考基因组、质谱输入表或基因集文件。

## 文档与验证

- [输入与运行说明](docs/INPUTS_AND_EXECUTION.md)
- [已知问题与限制](docs/KNOWN_ISSUES.md)
- [验证记录](docs/VALIDATION.md)
- [来源校验清单](SOURCE_MANIFEST.tsv)

已检查来源文件完整性、Python 语法及 notebook JSON 格式，尚未重跑科学分析。归档暂停指数脚本的链方向处理和 MNase notebook 的分组导出存在已记录问题，解释结果前应先阅读相关说明。

## 引用

相关稿件：*MED12 recruits SUPT6H to H2A.Z nucleosomes for pause release*。

正式发表后补充完整引用信息及 DOI。

## 许可证

当前尚未指定许可证，代码使用许可待作者确定。
