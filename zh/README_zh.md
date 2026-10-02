# 室内传染病流行的随机游走 SIR 模型（中文版正文）

本目录是论文《室内传染病流行的随机游走 SIR 模型：再生数、离散个体，以及以暴发记录、接触记录和示踪测量为对照的
预先设定检验》（周易潇潇，Xiaoxiao Zhouyi，英国布里斯托大学工程数学与技术学院，2026）**正文**的简体中文版：
[`paper_zh.pdf`](paper_zh.pdf)。代码、数据、英文正文与补充材料见仓库根目录的 [`README.md`](../README.md)。

- 中文版逐句翻译英文正文（[`../paper.pdf`](../paper.pdf)，源文件 `../main.tex` 与 `../sections/m1_intro.tex` …
  `m8_discussion.tex`），保留每一个数值、公式、限定条件和否定性结果；两者如有出入，以英文正文为准。
- 源文件中每一条 `% src:` 注释与英文正文相同（保持英文），其中的路径相对于仓库根目录。
- 补充材料不翻译：正文以“补充材料”加带前缀 S 的编号引用英文补充材料
  （[`../supplement.pdf`](../supplement.pdf)），编号由 `xr-hyper` 从 `../supplement.aux` 读取。

## 目录内容

| 路径 | 内容 |
|---|---|
| `paper_zh.pdf` | 中文版正文 |
| `main.tex` | 前置部分（标题、作者、摘要、关键词）、数据与代码可得性、伦理与数据、致谢、参考文献 |
| `sections/m1_intro.tex` … `m8_discussion.tex` | 第 1–8 节，文件名与英文正文一一对应 |
| `macros.tex` | 符号宏（与 `../macros.tex` 逐一相同）；单位词与定理环境名为中文 |
| `refs.bib`、`plainnat-zh.bst` | 参考文献库（与 `../refs.bib` 相同）；中文引用格式（“和”“等”） |
| `data/v2_tab_holdout2.tex` | 第 7 节一张表的表体，只译出文字，数值与英文版逐字相同 |
| `figures/*_zh.pdf` | 正文所用 5 幅图的中文标注版本 |
| `code/zh_figure_labels.py` | 以正文用语重绘一幅图的标注（不重新计算任何数据） |

图的来源：`fig3_scene_maps_zh.pdf`、`fig_r0_mobility_zh.pdf`、`fig5_1_finite_size_zh.pdf` 与
`fig5_5_R_vs_mobility_zh.pdf` 与 `../figures/` 中的同名文件相同，由绘制英文图的同一脚本根据同一组存储结果写出；
`fig_hotspot_maps_zh.pdf` 由 `code/zh_figure_labels.py` 调用 `../code/plan_fig_theory.py` 的绘图函数重绘，只把
(b) 栏标注改为“传染效率 $q$”。

## 如何编译

需要 TeX Live 2025（XeLaTeX、BibTeX；Fandol 字体随 TeX Live 提供）。在本目录中运行

```sh
latexmk -xelatex main.tex
```

输出为 `main.pdf`。英文补充材料的标签文件 `../supplement.aux` 随仓库提供，因此无需先编译英文文档；若重新编译了
英文补充材料，中文版读取新的 `../supplement.aux`。重绘那幅图：在仓库根目录运行
`python zh/code/zh_figure_labels.py`（只写入本目录的 `figures/`）。
