# 初稿构建报告

构建日期：2026-07-14。

## 交付物

- 主文件：`main.tex`
- 稳定 PDF：`output/pdf/wigner_manybody_dynamics_draft.pdf`
- 章节：18
- 附录：4
- PDF 页数：126
- PDF SHA-256：`B9CD70C6DE7D6EED7CABB28ABFB9E8F7439CA1AB66ADA3C183BB4C738DE4E7C6`
- 参考文献：17
- 实际数值图：9

## 验证结果

- `C:\ProgramData\anaconda3\envs\AI_group\python.exe code\run_all.py`：9/9 通过。
- `C:\ProgramData\anaconda3\envs\AI_group\python.exe code\audit_manuscript.py`：通过。
- `pwsh -NoProfile -File .\build.ps1`：通过。
- `main.log` 与 `main.blg`：无未定义引用、重复标签、溢出框或 LaTeX/BibTeX 警告。
- PDF 跨部分视觉抽查：通过；扩写后的卷首使用指南另经逐页渲染检查，详见 `reviews/full_pdf_visual_review.md` 与 `reviews/how_to_read_review.md`。
- DOI 核验：17/17 通过，详见 `reviews/reference_verification.md`。

## 环境说明

MiKTeX 尝试写用户日志目录时会输出 `log4cxx ... 拒绝访问` 环境提示，但所有构建命令退出码为 0，XeLaTeX、BibTeX、xdvipdfmx 和最终 PDF 均正常。该提示不属于书稿警告。
