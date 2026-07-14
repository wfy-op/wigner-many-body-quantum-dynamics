# 严格修订稿构建报告

构建日期：2026-07-14。

## 交付物

- 主文件：`main.tex`
- 稳定 PDF：`output/pdf/wigner_manybody_dynamics_draft.pdf`
- 日期命名 PDF：`output/pdf/维格纳相空间方法与玻色多体量子动力学_2026-07-14.pdf`
- GitHub 附件 PDF：`output/pdf/Wigner_Phase-Space_Methods_and_Bosonic_Many-Body_Quantum_Dynamics_2026-07-14.pdf`
- 章节：18
- 附录：6
- PDF 页数：155
- PDF 字节数：1,204,542
- PDF SHA-256：`2EBC06710ECAA8833EF04F36075FDD43182339F9F9DE64CCFCAB8CFA06B9576E`
- 参考文献：41
- 程序生成的证据图：10

## 验证结果

- 干净来源提交：`393bf7bbfdc7ad87392c9eeaaaa40f6f3e0f57d2`。
- `code/run_all.py`：10/10 个程序通过，共72项递归机器检查；全部metrics为本轮新产物。
  `data/test_summary.json` 的 `environment_at_start.git` 记录上述提交、`dirty=false` 和空状态列表。
- `tests/test_run_all_freshness.py`：5/5 通过，包括“脚本退出 0 但不产生新 metrics 必须失败”。
- `code/audit_manuscript.py`：通过；18 章、6 附录、243 标签、41 引文键和 10 幅图闭合。
- `code/verify_references.py`：41/41 条 DOI 与书目信息通过 Crossref 核验。
- `pwsh -NoProfile -File .\build.ps1`：通过，生成三份同内容稳定 PDF。
- `main.log` 与 `main.blg`：无未定义引用、重复标签、溢出框或 LaTeX/BibTeX 警告。
- PDF 跨部分视觉抽查：通过，详见 `reviews/full_pdf_visual_review.md`。

## 可复现性说明

本报告记录“干净来源提交 → 全量数值复跑 → 静态审计 → 文献核验 →
完整构建 → PDF哈希”的闭环。通用名、中文日期名和 GitHub ASCII 日期名三份
PDF 的字节数与 SHA-256 完全一致。当前修订稿已发布为
[`v0.2.0-2026-07-14`](https://github.com/wfy-op/wigner-many-body-quantum-dynamics/releases/tag/v0.2.0-2026-07-14)；
`v0.1.0-2026-07-14` 保留为严格审阅前的126页历史快照。

MiKTeX 尝试写用户日志目录时会输出 `log4cxx ... 拒绝访问` 环境提示，但所有
构建命令退出码为 0，XeLaTeX、BibTeX、xdvipdfmx 和最终 PDF 均正常。该提示不属于书稿警告。
