# 严格修订稿构建报告

构建日期：2026-07-14。

## 交付物

- 主文件：`main.tex`
- 稳定 PDF：`output/pdf/wigner_manybody_dynamics_draft.pdf`
- 日期命名 PDF：`output/pdf/维格纳相空间方法与玻色多体量子动力学_2026-07-14.pdf`
- 章节：18
- 附录：6
- PDF 页数：155
- PDF SHA-256：将在最终干净来源复跑后锁定
- 参考文献：41
- 程序生成的证据图：10

## 验证结果

- `code/run_all.py`：10/10 个程序通过，共 72 项递归机器检查；全部 metrics 为本轮新产物。
- `tests/test_run_all_freshness.py`：5/5 通过，包括“脚本退出 0 但不产生新 metrics 必须失败”。
- `code/audit_manuscript.py`：通过；18 章、6 附录、243 标签、41 引文键和 10 幅图闭合。
- `code/verify_references.py`：41/41 条 DOI 与书目信息通过 Crossref 核验。
- `pwsh -NoProfile -File .\build.ps1`：通过，生成两份同内容稳定 PDF。
- `main.log` 与 `main.blg`：无未定义引用、重复标签、溢出框或 LaTeX/BibTeX 警告。
- PDF 跨部分视觉抽查：通过，详见 `reviews/full_pdf_visual_review.md`。

## 可复现性说明

最终证据提交将记录“干净来源提交 → 全量数值复跑 → 静态审计 → 文献核验 →
完整构建 → PDF 哈希”的闭环。旧的 GitHub Release v0.1.0-2026-07-14 是严格
审阅前的 126 页历史快照；本轮未重新推送或发布外部 Release。

MiKTeX 尝试写用户日志目录时会输出 `log4cxx ... 拒绝访问` 环境提示，但所有
构建命令退出码为 0，XeLaTeX、BibTeX、xdvipdfmx 和最终 PDF 均正常。该提示不属于书稿警告。
