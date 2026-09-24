# 共模治理方案研究：Y 电容到 PE vs 共模扼流圈

**前提**：`docs/07`（强制）。⚠️ 本研究的 Y 电容方案**新增了板侧 PE 通路**，**偏离 docs/07 的"板侧无 PE"实测前提**，属**假设性改进方案**，非当前实测配置。

## 结论
- **1 nF/2kV 到 PE：基本无效**。@0.805 MHz 仅 **+0.34 dB**；**@29.8 MHz 反而恶化 7.45 dB**（与电缆电感反谐振）。
- **要达标需 ~245 nF** → **50Hz 漏电 ≈17.7 mA**，**远超 Class I 的 3.5 mA 限值**，**安规不可行**。
- **共模扼流圈 1 mH**：0.81 MHz **−0.28 dB（无效）**、30 MHz **−51 dB（极有效）** —— 与 Y 电容趋势相反，**两者都救不了 0.81 MHz**。
- **🚩 结构性提示**：本模型把 CM 源当**理想电流源**（源阻抗 = C_p ≈19.6 kΩ@0.81MHz）→ 回路上加任何元件都改不了电流 → 结论"只能降 C_p"。**真实共模滤波有效** ⇒ **C_p=10pF 假设很可能过大**。结论受 C_p 支配（见 `results/07-conducted-cm/README.md` 的 C_p 敏感性表）。

## 产物
- `figures/fig_emi_cm_ycap_compare.png` — 基线 vs 1nF vs 最优 C_y（三工况 + 限值线）
- `figures/fig_emi_cm_ycap_sweep.png` — 最差裕量 vs C_y（三工况）
- `reports/REPORT_emi_cm_ycap.md`、`reports/emi_numbers_ycap.txt`
