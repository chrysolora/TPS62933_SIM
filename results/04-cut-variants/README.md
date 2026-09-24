# 04 — 电感下铺铜三版涡流场对比（FEA）

| 文件 | 说明 |
|---|---|
| `figures/compare_top.png` | 三版顶视 \|J\| 并排（原 ParaView 图） |
| `figures/compare_top_light.png` | **可读重绘**（浅色标 log 1e3–5e7；上=顶层，下=底层） |
| `figures/compare_radial.png` | 三版径向 \|J\| 剖面（顶+底合并，起点=L2 中心） |
| `figures/fig_field_top_*.png` | 各版单独顶视图 |
| `reports/REPORT.md` | 子代理报告（方法/数值/局限） |
| `reports/metrics_*.json` | 各版实测指标（\|J\|max、分环带均值、峰值半径、面积占比） |
| `scripts/` | `gen_geo4.py`（建模）`make_sif4.py`（求解器）`extract4.py`（提铜层）`radial4.py`（剖面）`compare4.py`（对比）`replot4.py`（可读重绘） |

**核心数值**：顶层 \|J\|max FullCopper 4.89e7 / TopCutout 2.00e6 / DualCutout 2.64e6 A/m²；
电感下(r<3.25mm)均 \|J\| 5.3e5 → 0 → 0；合并均 \|J\| 2.58e5 → 1.75e4(−93%) → ~5.8e3(−98%)。

**局限**：板厚 1.0mm（估值）、µr=40（fea3 估值）、REGION 外接矩形近似、网格较 fea3 略粗。
**口径提醒**：本目录 FullCopper 为「顶+底两层」重建，与 `results/03`（仅顶层，\|J\|max 2.9e7）不同口径。
