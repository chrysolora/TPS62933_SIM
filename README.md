# TPS62933 Buck 仿真结果集

TI **TPS62933**（24V → 12V，最大 3A）降压电路的仿真成品结果，含**仿真方法**、**仿真过程**与**结果**。

> 说明：本仓库只收"**成品/可用**"结果，每个结果标注**方法、可信度与局限**。
> 尚未定稿的迭代（含寄生 DM-EMI v2、真实铺铜几何场仿真 fea3）另行标注为 **pending**。

---

## 目标电路（板上 U21）

| 项 | 值 | 来源 |
|---|---|---|
| 芯片 | TPS62933**F**DRLR（F 版，带 SS 脚） | 原理图 |
| 输入 / 输出 | 24V → 12V，Iout ≤ 3A | 需求 |
| 开关频率 fsw | 805 kHz（实测网表）/ 784 kHz（RT=26.7kΩ 公式） | 仿真 / 数据手册 |
| 输入滤波 | `24VIN—F1—[3×10µF]—L1(1µH)—[2×10µF‖100nF]—U21.VIN` + RC 阻尼 | 原理图（已核对） |

---

## 交付索引

| # | 结果 | 目录 | 状态 | 可信度 |
|---|---|---|---|---|
| 01 | 时域瞬态（纹波/开关/效率） | `results/01-stage1-transient/` | ✅ 交付 | ⚠️ **工程近似** |
| 02 | 频域传导 EMI（DM, CISPR 32 近似）+ 软启动 | `results/02-conducted-emi/` | ✅ **v2 定稿** | ⚠️ 含寄生后**贴限值**（见下） |
| 02b | DM 敏感性 + 无 π 滤波对比 | `results/02-conducted-emi/` | ✅ v2 | ⚠️ 趋势可信 |
| 03 | 电感下铺铜涡流场（FEA） | `results/03-field-eddy-copper/` | ✅ **v3 真几何** | ⚠️ 参数为估值 |
| 04 | **三版铺铜涡流对比**（FullCopper/TopCutout/DualCutout） | `results/04-cut-variants/` | ✅ | ⚠️ 相对可信 |
| 05 | **板级辐射发射**（openEMS FDTD）+ 三版对比 | `results/05-radiated/` | ✅ | ⚠️ 仅相对/热点 |
| 06 | **三版总汇总（涡流 vs 辐射）** | `docs/06-multi-variant-summary.md` | ✅ | ⚠️ 务读 |
| — | 方法/过程/局限文档 | `docs/` | ✅ | — |

### 最新结论（2026-09-24）

- **DM EMI（含寄生）**：v1 的 +34.4dB 是**假象**（缺走线电感）；v2 **标称 +15.9dB / 悲观 +0.7dB**（几乎贴 CISPR 32 Class B @4.03MHz）。**主导寄生 = 走线/回路电感**。**π 滤波挡掉 ~58dB**。CM 尚未正式建模（很可能是更大风险）。
- **场（真几何）**：板框 26.50×50.50mm，铜=板框−铺铜（占 81.2%）；|J|max≈2.9e7 A/m² @r≈4.7mm，径向单调衰减；>1e7 仅 0.073% 节点（无大面积斑点）。
- **三版涡流（fea4）**：电感下(r<3.25mm)均|J| FullCopper 5.3e5 → TopCutout **0** → DualCutout **0**；合并均|J| **−93% / −98%**。
- **三版辐射（rad）**：E@3m — TopCutout **+2.5dB**、DualCutout **+3.1dB**（**挖铜恶化辐射**）。
- ➡️ **总汇总（两者相反）**：`docs/06-multi-variant-summary.md`　**挖铜减损耗、但增辐射**。

**方法**：`docs/01-simulation-method.md`　**过程/时间线**：`docs/02-simulation-process.md`　**局限**：`docs/04-known-issues-and-limits.md`

---

## 一句话结论（TL;DR）

- **时域**：可用 ngspice 收敛，但**免费未加密模型有物理边界**（VIN 脚电容 >~30µF 不启动）→ 只能用缩版 π 滤波、跳过 14ms 软启动 ⇒ **结果定性可用、定量不可信**（详见 01）。
- **传导 EMI（DM）**：含已确认的真实 π 滤波时，DM **远低于 CISPR 32 Class B**；去掉 π 滤波则明显抬高。**但真正的瓶颈很可能是 CM（共模）**，本模型尚未覆盖。
- **场**：已打通「PCB → Gmsh → Elmer 3D 涡流 → ParaView 出图」链路；几何/磁芯参数为简化+估值，**用于趋势判断，不用于定量**。
- **三版铺铜（涡流+辐射）**：**挖铜大幅降涡流损耗（−93%/−98%），但恶化辐射（+2.5/+3.1 dB）——两者相反**。辐射绝对值不可信，相对差可用。见 `docs/06-multi-variant-summary.md`。

---

## 复现环境

- 计算主机：内网 `node249`（Ubuntu 24.04，12 核 / 15G，数据盘 RAID10 `/mnt/raid10`）
- 工具链：**ngspice 42**、**Elmer FEM 9.0**、**Gmsh 4.12**、**ParaView 5.11**、Python(numpy/scipy/matplotlib)
- 无头渲染：`LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a -s '-screen 0 800x600x24' pvbatch ...`

---

_生成于 2026-09-24。_
