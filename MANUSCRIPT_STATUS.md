# 书稿状态

状态标记：`planned` / `drafted` / `reviewed` / `verified`。

| 章 | 标题 | 草稿 | 审核 | 数值/推导验证 |
|---|---|---|---|---|
| 01 | 从量子态到可计算动力学 | drafted | reviewed | verified |
| 02 | Weyl 变换与 Wigner 函数 | drafted | reviewed | verified |
| 03 | 星积、Moyal 括号与算符对应 | drafted | reviewed | verified |
| 04 | 常见量子态、非经典性与观测量 | drafted | reviewed | verified |
| 05 | 二次哈密顿量的精确相空间流 | drafted | reviewed | verified |
| 06 | Kerr 振子：非线性与量子回复 | drafted | reviewed | verified |
| 07 | 截断维格纳近似的推导 | drafted | reviewed | verified |
| 08 | Bose--Hubbard 二聚体 | drafted | reviewed | verified |
| 09 | TWA 的有效范围与失效模式 | drafted | reviewed | verified |
| 10 | 从玻色模式到离散量子场 | drafted | reviewed | verified |
| 11 | 初态的 Wigner 抽样 | drafted | reviewed | verified |
| 12 | 经典场传播与数值可信度 | drafted | reviewed | verified |
| 13 | 观测量、关联函数与误差条 | drafted | reviewed | verified |
| 14 | 双分量凝聚体的模型与对称性 | drafted | reviewed | verified |
| 15 | Bogoliubov 模式与淬火协方差 | drafted | reviewed | verified |
| 16 | 相位扩散、模式耦合与预热化 | drafted | reviewed | verified |
| 17 | 条纹序、有限尺寸与超固态判据 | drafted | reviewed | verified |
| 18 | 开放系统、替代方法与研究工作流 | drafted | reviewed | verified |

## 全书构建

- 正文与附录：18 章、6 个附录，均为非占位初稿。
- 参考文献：41/41 条 DOI 与书目信息通过 Crossref 核验。
- 数值证据：10 个程序由 `code/run_all.py` 自动运行，10/10 通过，共 72 项递归机器检查。
- 静态审计：243 个标签无重复，10 幅证据图均有编号、图注和正文交叉引用。
- 完整 PDF：155 页；XeLaTeX/BibTeX 日志无未定义引用、重复标签、溢出框或警告。
- 视觉审核：封面、目录、图目录、核心推导、端到端案例、附录 A--F 与参考文献均已渲染抽查。

## 当前结论与后续边界

严格审阅中的 R01--R16 已在正文、程序或审核证据中关闭。本地修订稿已经达到
可发布的研究生教材级初稿标准。正式出版前仍需作者确定署名、机构和版权页，
邀请独立领域专家复核，并在获得真实参数后扩展 Raman SOC、开放非线性场和
MPS/实验交叉验证。当前版本已发布为 GitHub Release `v0.2.0-2026-07-14`。
