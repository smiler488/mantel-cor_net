# mantel-cor_net

[English](README.md) | 中文

一个 R 脚本，用来画土壤和作物论文里常见的"相关性热图 + Mantel 连线"图。上三角是各解释变量之间的 Pearson 相关；每个响应变量（产量、磷肥利用率等）都和每个解释变量连一条线，表示两者的 Mantel 检验结果：颜色表示 p 值，线宽表示 |r|，虚线表示负相关。

这是我给一个磷肥试验做的图，想把不同土层的速效磷、根系和植株的磷性状、产量和 PUE 放在同一张图里。脚本只用 ggplot2 写成，布局的每个元素都可以自己改。

## 运行环境

R 4.0 或更高版本，需要 `ggplot2`、`vegan`、`ggnewscale`、`svglite`。第一次运行时，缺少的包会自动从 CRAN 安装。

## 使用

```bash
Rscript cor_mantel_overlay.R                    # 使用 demo.csv
Rscript cor_mantel_overlay.R mydata.csv fig3    # 输出 fig3.pdf、fig3.svg、fig3_mantel.csv
```

在脚本开头设置变量列表：

```r
vars    <- c("Year", "cultivar", "Treatment", "0-10cm AP", ...)   # 热图变量，也是连线的终点
targets <- c("Yield", "PUE", "SB")                                 # 响应变量
```

`demo.csv` 是 288 个小区 × 22 列的示例数据，格式就是脚本要求的格式，数值是为演示编造的。

## 计算内容

- **热图**：每对数值变量的 Pearson r 和 `cor.test` 的 p 值，`*` 表示 p < 0.05，`**` 表示 p < 0.01，`***` 表示 p < 0.001。涉及分类变量（品种、处理）的格子留空，因为对因子编码算 r 没有意义。
- **Mantel 连线**：对每个响应变量和每个解释变量，用 `vegan::mantel` 检验两个距离矩阵（Pearson，999 次置换，随机种子 123）。数值变量的距离取 |a − b|；分类变量水平相同记 0，不同记 1。
- `<前缀>_mantel.csv` 列出全部 Mantel r 和 p。

Mantel 检验的 p 值没有做多重比较校正。3 个响应变量 × 19 个解释变量，一共 57 次检验，解读单条连线之前，建议先对 CSV 里的 p 值做校正，比如 `p.adjust(..., "BH")`。

## 本次修订

- 分类变量原来被转成了任意的数字编码，导致热图里出现没有意义的 Pearson r，Mantel 检验用的也是编码之差。现在分类变量不进入热图，Mantel 检验改用"水平相同为 0、不同为 1"的距离。
- 线宽原来按 r 分档，而不是按 |r|，所以很强的负相关会被画成细线。
- 颜色图例原来写的是"0.001 < p ≤ 0.01"，但这一档其实也包含了 p ≤ 0.001 的情况。
- 响应变量节点的位置改为按个数自动排布，不再写死只适用于 Yield、PUE、SB 三个变量。
- 输入文件和输出前缀可以从命令行指定，缺少列时会给出明确的报错。
- Mantel 检验结果会额外导出为 CSV。
- 原 README 提到一个仓库里并不存在的辅助脚本，已删去。

## 引用

```bibtex
@software{deng_mantel_cor_net,
  author = {Deng, Liangchao},
  title  = {mantel-cor\_net: correlation heatmap with Mantel-test links in R},
  year   = {2025},
  url    = {https://github.com/smiler488/mantel-cor_net}
}
```

## 许可证

MIT，见 [LICENSE](LICENSE)。邓良超（Liangchao Deng）。
