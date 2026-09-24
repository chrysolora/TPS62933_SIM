# 传导共模 (CM) EMI — 150 kHz–30 MHz

**前提**：`docs/07-input-conditions-and-assumptions.md`（强制）。板侧无 PE、仅 +24V/GND 两根线、无 Y 电容、塑料外壳无屏蔽。

## 结论（TL;DR）
- **CM 才是真正的瓶颈**：三工况**全部超出 CISPR 32 Class B**，满载最差 **−28.0 dB @ 0.81 MHz**。
- **CM 比 DM 高约 71 dB**（满载 DM 仅 13.2 dBµV / +39.9 dB 裕量）→ **π 滤波只治 DM，对 CM 无效**。
- CM 由 **SW 节点 dv/dt 经寄生电容 C_p 注入**，**基本不随负载电流变化**（空载 −25.2 / 半载 −29.3 / 满载 −28.0 dB）→ 印证 CM 由**电压跳变**主导，不是电流。
- C_p 敏感性：2 pF→−14.0 dB；10 pF→−28.0；50 pF→−42.0 dB ⇒ **即使最乐观 2 pF 仍超标 14 dB**。

## 图（按要求只保留 1 张）
- `figures/fig_emi_cm_cispr.png` — **CM 电压 vs CISPR 32 B 限值（三工况）**（唯一保留图）
- （CM+DM 合成图、C_p/tr 敏感性图已移除；数值见 `reports/emi_numbers_cm.txt` 与报告 §4/§5）

## 🔁 更新（2026-09-24）
本目录的**黑箱 C_p=10 pF** 已被 **[`results/10-cm-pcb-geometry/`](../10-cm-pcb-geometry/) 取代**（用真实铜皮几何 + Elmer 静电场算 C_sw ≈ 0.2–0.4 pF）。
→ 本目录的"超限 25–29 dB"是**悲观上界**；按几何值约 **+28…+34 dB 改善**，FullCopper 满载**临界（−0.3dB）**、挖铜后转正。**结论引用请以 10 为准。**

## ⚠️ 重要局限：本传导 CM 模型**未结合 PCB 铜皮布局**
- 本模型用**集总寄生电容 C_p（10 pF，估）**代表"板↔大地"耦合，**是黑箱耦合，不是真实铜皮几何**。
- **三种铺铜方案（FullCopper/TopCutout/DualCutout）在本模型中无区别** —— 传导 CM 结果**不随铜皮变化**。
- 真正"结合 PCB 铜皮"的只有**辐射**部分（openEMS，见 `results/05-radiated/`、`results/06-rad-3d/`）。
- 若要让传导 CM 反映布局（SW 铜面积、回路面积、C_p 随布局变化）→ 需专门的**场提取/3D 电磁**建模，属后续工作。

## 需处理（工程建议）
共模扼流圈 / 降 C_p（缩短 SW 回路、加屏蔽/贴地）/ 必要时加 Y 电容（但本板无 PE，Y 电容通路受限）→ 专门治 CM。

## 局限
C_p 为估算（唯一闭合 CM 回路的元件）；SW 谱含仿真振铃伪像（已与理想 5 ns 梯形波并报交叉核对）；仅建模 SW 单一源（未含上游 65 kHz 侧与 12 V 输出侧 CM）；直接比 QP 限值（偏保守）。详见 `reports/REPORT_emi_cm.md`。
