# 参考文献核验记录

核验日期：2026-07-14。

## 方法

- 按 `nature-academic-search` 的 citation verification 工作流，从 `bibliography/references.bib` 提取 DOI、题名、期刊与年份。
- 首选 Crossref T1 结构化元数据；本会话未提供论文检索 MCP 工具，因此改用只读 Crossref REST API。
- DOI 统一转为小写比较，题名忽略大小写与排版标记，期刊历史更名保留原发表时常用名称。
- 分类：`verified` / `mismatch` / `not_found` / `suspicious` / `manual_needed`。

## 汇总

- 总数：17。
- verified：17。
- mismatch：0。
- not_found：0。
- suspicious：0。
- manual_needed：0。

## 明细

| 键 | DOI | 状态 | 说明 |
|---|---|---|---|
| Wigner1932 | 10.1103/PhysRev.40.749 | verified | 题名、Physical Review、1932、40、749–759 匹配。 |
| Moyal1949 | 10.1017/S0305004100000487 | verified | 题名、1949、45、99–124 匹配；Crossref 显示期刊现名 Mathematical Proceedings，书稿保留历史常用刊名。 |
| Hillery1984 | 10.1016/0370-1573(84)90160-1 | verified | Physics Reports 106, 121–167 匹配。 |
| Steel1998 | 10.1103/PhysRevA.58.4824 | verified | Physical Review A 58, 4824–4835 匹配。 |
| Sinatra2002 | 10.1088/0953-4075/35/17/301 | verified | Journal of Physics B 35, 3599–3631 匹配。 |
| Polkovnikov2010 | 10.1016/j.aop.2010.02.006 | verified | Annals of Physics 325, 1790–1852 匹配。 |
| Blakie2008 | 10.1080/00018730802564254 | verified | 按 Crossref 将题名标准化为 “using c-field techniques”。 |
| Lin2011 | 10.1038/nature09887 | verified | Nature 471, 83–86 匹配。 |
| Li2012 | 10.1103/PhysRevLett.108.225301 | verified | Physical Review Letters 108, 225301 匹配。 |
| Zhai2015 | 10.1088/0034-4885/78/2/026001 | verified | Reports on Progress in Physics 78, 026001 匹配。 |
| Li2017 | 10.1038/nature21431 | verified | Nature 543, 91–94 匹配。 |
| Boninsegni2012 | 10.1103/RevModPhys.84.759 | verified | Reviews of Modern Physics 84, 759–776 匹配。 |
| Hohenberg1967 | 10.1103/PhysRev.158.383 | verified | Physical Review 158, 383–386 匹配。 |
| Gring2012 | 10.1126/science.1224953 | verified | Science 337, 1318–1322 匹配。 |
| Lindblad1976 | 10.1007/BF01608499 | verified | Communications in Mathematical Physics 48, 119–130 匹配。 |
| Drummond1980 | 10.1088/0305-4470/13/7/018 | verified | Journal of Physics A 13, 2353–2368 匹配。 |
| Vidal2003 | 10.1103/PhysRevLett.91.147902 | verified | Physical Review Letters 91, 147902 匹配。 |

## 引用边界

- 第14章引用 Raman SOC 的实验、相图与综述，不把均匀 `J` 模型冒充 SOC。
- 第16章用 Gring 等人的实验支持“预热化可在孤立一维量子气体中观测”，不据此证明本书模型必然热化。
- 第17章用超固态综述、SOC 条纹实验与低维长程序定理支持判据边界，不把玩具标度程序当作具体实验复现。
- 第18章用 Lindblad、positive-P 和 MPS 原始来源支持方法介绍。
