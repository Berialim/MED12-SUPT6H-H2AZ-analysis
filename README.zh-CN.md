# MED12 SUPT6H H2A.Z 分析代码

根据 2026 年 10 月 17 日稿件整理，收录本机找到的部分对应脚本。英文 README 可作为 GitHub 首页。

代码按质谱、ChIP-seq、PRO-seq、TT-seq、MNase、小鼠 RNA-seq 和外部数据比较分类。原始文件位于 `archive/`，保留原样；`SOURCE_MANIFEST.tsv` 记录相对于本机 `D:/work/pby` 的来源及 SHA-256。

图号对应是根据代码内容推断的候选关系，并未证明这些文件就是最终论文图的生成版本。原始数据和输入表未打包，也未重新运行分析。请先阅读 `docs/KNOWN_ISSUES.md`。

特别需要处理：PRO-seq 负链坐标、MNase 四分位与二分组导出不一致、GSEA 排名与方法描述差异、DESeq2 阈值及 apeglm 收缩版本、鼠源 mm39 注释及对应输入。不要将此仓库描述为论文完整可复现代码。

整理时检索了 `D:/work/pby` 以及相关的 `D:/work/lailab2025`、`D:/work/med12_rpb1_apex2_20230923`；未遍历电脑全部目录。没有修改原始分析文件。GitHub 仓库：https://github.com/Berialim/MED12-SUPT6H-H2AZ-analysis 。仓库可见性沿用 GitHub 上的现有设置。

