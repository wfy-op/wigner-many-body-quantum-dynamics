# 维格纳相空间方法与玻色多体量子动力学

本仓库包含中文研究生教材的 LaTeX 初稿、逐章审核记录、数值程序和可复核证据。
当前严格修订稿包括 18 章、6 个附录、10 幅程序生成的证据图和 41 条已核验参考文献；
最新本地构建为 155 页。

## PDF 与版本说明

本地构建会同时生成：

```text
output/pdf/wigner_manybody_dynamics_draft.pdf
output/pdf/维格纳相空间方法与玻色多体量子动力学_2026-07-14.pdf
```

[GitHub Release v0.1.0-2026-07-14](https://github.com/wfy-op/wigner-many-body-quantum-dynamics/releases/tag/v0.1.0-2026-07-14)
是严格审阅前的 126 页历史快照，不包含本轮对 BdG 约定、端到端空间场案例、
证据门禁、附录 E/F 和图目录的系统修订。本轮未重新推送或发布 GitHub Release；
因此当前修订稿以本仓库的本地 PDF 和证据文件为准。

## 构建

在项目根目录运行：

```powershell
pwsh -NoProfile -File .\build.ps1
```

构建脚本依次运行 XeLaTeX、BibTeX 和最终 PDF 生成步骤，并复制上述两个稳定文件名。

## 验证

运行全部实际数值基准：

```powershell
C:\ProgramData\anaconda3\envs\AI_group\python.exe .\code\run_all.py
```

运行书稿完整性审计与参考文献核验：

```powershell
C:\ProgramData\anaconda3\envs\AI_group\python.exe .\code\audit_manuscript.py
C:\ProgramData\anaconda3\envs\AI_group\python.exe .\code\verify_references.py
```

最近一次完整科学计算由 10 个程序给出 72 项递归机器检查，10/10 个程序均通过。
每项程序都把预先固定的机器可读阈值、逐项判据和
`validation.all_passed` 写入 JSON。这里的“通过”仅表示各自声明范围内的结构化
基准通过：均匀 Rabi 双分量场案例不是 Raman SOC 模拟，第17章的条纹标度程序
仍是统计玩具模型。

总入口不会把旧 JSON 冒充本轮结果。每个子程序启动前都会记录墙钟起点和旧
metrics 的存在性、`mtime_ns` 与 SHA-256，并向子程序传入唯一 `run_id`；子程序
必须重新写入统一 schema、该 `run_id`、脚本 SHA-256 与运行时间。返回后，入口
要求 metrics 的修改时间不早于本轮起点；若旧文件存在，还要求修改时间推进且
SHA-256 改变。入口随后递归读取简单布尔检查和富结构中的 `passed`，重算总判据
并与 `validation.all_passed` 严格比对。缺文件、陈旧文件、schema/run 元数据不符
或判据不一致都会使总入口以失败退出，即使子程序返回码为零。

汇总文件 [`data/test_summary.json`](data/test_summary.json) 保存每个脚本与 metrics
的 SHA-256、metrics schema、新旧文件状态和递归判据，同时在任何基准改写产物前
记录 Python/关键库版本、操作系统、Git 提交与工作区起始状态。Kerr 与二聚体曲线同时保存 Monte Carlo
标准误和独立子随机流稳定性结果；图中的阴影来自实际轨迹计算，而不是视觉示意。

## 目录约定

- `chapters/`：分章正文。
- `appendices/`：数学、抽样、数值、习题、文献导读和案例地图。
- `reviews/`：逐章审核、全书严格审阅与闭环记录。
- `code/`：产生书中数值证据的程序及统一门禁。
- `tests/`：证据入口的回归测试。
- `data/`：程序产生的轻量数据和审计账本。
- `figures/`：程序生成或原创的书稿图。
- `bibliography/`：可核验参考文献。
- `output/pdf/`：最终交付 PDF。

原始需求文件 `prompt_wigner.txt` 保持不变。
