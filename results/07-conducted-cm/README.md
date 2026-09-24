# 传导共模 (CM) EMI — 150 kHz–30 MHz

**前提**：`docs/07-input-conditions-and-assumptions.md`（强制）。板侧无 PE、仅 +24V/GND 两根线、无 Y 电容、塑料外壳无屏蔽。

## 结论（TL;DR）
- **CM 才是真正的瓶颈**：三工况**全部超出 CISPR 32 Class B**，满载最差 **−28.0 dB @ 0.81 MHz**。
- **CM 比 DM 高约 71 dB**（满载 DM 仅 13.2 dBµV / +39.9 dB 裕量）→ **π 滤波只治 DM，对 CM 无效**。
- CM 由 **SW 节点 dv/dt 经寄生电容 C_p 注入**，**基本不随负载电流变化**（空载 −25.2 / 半载 −29.3 / 满载 −28.0 dB）→ 印证 CM 由**电压跳变**主导，不是电流。
- C_p 敏感性：2 pF→−14.0 dB；10 pF→−28.0；50 pF→−42.0 dB ⇒ **即使最乐观 2 pF 仍超标 14 dB**。

## 图
- `figures/fig_emi_cm_cispr.png` — CM 电压 vs CISPR 32 B 限值（三工况）
- `figures/fig_emi_cm_dm_total.png` — CM+DM 合成 vs 限值
- `figures/fig_emi_cm_sensitivity.png` — C_p / tr / C_y 敏感性

## 需处理（工程建议）
共模扼流圈 / 降 C_p（缩短 SW 回路、加屏蔽/贴地）/ 必要时加 Y 电容（但本板无 PE，Y 电容通路受限）→ 专门治 CM。

## 局限
C_p 为估算（唯一闭合 CM 回路的元件）；SW 谱含仿真振铃伪像（已与理想 5 ns 梯形波并报交叉核对）；仅建模 SW 单一源（未含上游 65 kHz 侧与 12 V 输出侧 CM）；直接比 QP 限值（偏保守）。详见 `reports/REPORT_emi_cm.md`。
