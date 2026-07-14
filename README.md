# 维格纳相空间方法与玻色多体量子动力学

本目录包含中文研究生教材的 LaTeX 初稿、逐章审核记录、数值程序和构建产物。

## 下载 PDF

编译好的 2026-07-14 初稿可从 [GitHub Release](https://github.com/wfy-op/wigner-many-body-quantum-dynamics/releases/tag/v0.1.0-2026-07-14) 下载。附件名为 `维格纳相空间方法与玻色多体量子动力学_2026-07-14.pdf`。

## 构建

在项目根目录运行：

```powershell
pwsh -NoProfile -File .\build.ps1
```

成功后，稳定文件位于：

```text
output/pdf/wigner_manybody_dynamics_draft.pdf
```

当前完整初稿为 18 章、4 个附录、126 页。

## 验证

运行全部实际数值基准：

```powershell
C:\ProgramData\anaconda3\envs\AI_group\python.exe .\code\run_all.py
```

运行书稿完整性审计：

```powershell
C:\ProgramData\anaconda3\envs\AI_group\python.exe .\code\audit_manuscript.py
```

最新结果：9/9 数值基准通过，18/18 章有审核记录，17/17 条参考文献已核验。

## 目录约定

- `chapters/`：分章正文。
- `appendices/`：数学、抽样、数值和习题提示。
- `reviews/`：逐章审核记录。
- `code/`：产生书中数值证据的程序。
- `data/`：程序产生的轻量数据。
- `figures/`：程序生成或原创的书稿图。
- `bibliography/`：可核验参考文献。
- `output/pdf/`：最终交付 PDF。

原始需求文件 `prompt_wigner.txt` 保持不变。
