> ⚠️ **已被取代（SUPERSEDED，2026-09-24）**：本目录结果基于**旧前提**（无参考地平面、无真实电缆/CM 通路、无上游 LM50-20B24 参数、源仅本板开关）。
> **更新前提后的三版辐射对比请看 [`results/09-radiated-b2/`](../09-radiated-b2/)**（含参考地平面 + 10cm 电缆 + CM 通路 + 上游 65kHz/150mV 源）。
> 本目录仅作历史记录，**结论不要引用**。

# 05 — 板级辐射发射仿真（openEMS FDTD）

| 文件 | 说明 |
|---|---|
| `figures/fig_rad_compare_freq.png` | **三版远场 vs 频率**（叠 CISPR 32 Class B @3m 限值线） |
| `figures/fig_rad_compare_hotspot.png` | 三版辐射**方向图热点并排**（50.7 MHz） |
| `figures/fig_farfield_vs_cispr.png` | FullCopper 单版远场 vs CISPR |
| `figures/fig_farfield_3d.png` | FullCopper 3D 远场方向图 |
| `reports/REPORT_rad.md` | 报告（方法/归一化校验/数值/局限） |
| `scripts/` | `board_b1.py`（B1 板级建模）`build_rad.py`（三版建模）`post_b1.py`/`post_compare.py`（出图） |

**核心数值**：E@3m（±2 谐波包络）FullCopper 24.4/22.1/12.0 dBµV/m @30/100/300MHz；
**TopCutout +2.2~2.5 dB，DualCutout +3.0~3.1 dB** → 电感下挖铜**恶化**辐射。

**局限**：绝对幅值不可信（无共模路径/真实电缆/器件寄生）；**FDTD 未收敛到 −40dB**（三版设定一致，相对结论成立）；
GND pour 形状近似；激励热回路面积为估值；SW dv/dt 源本版未加。
