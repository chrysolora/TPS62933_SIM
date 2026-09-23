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
| 02 | 频域传导 EMI（DM, CISPR 32 近似）+ 软启动 | `results/02-conducted-emi/` | ✅ 交付 | ⚠️ **理想化、偏乐观** |
| 03 | 电感下铺铜涡流场（FEA） | `results/03-field-eddy-copper/` | ⚠️ 链路验证 | ⚠️ **局部简化几何** |
| 04 | 方法/过程/局限文档 | `docs/` | ✅ | — |

**方法**：`docs/01-simulation-method.md`　**过程/时间线**：`docs/02-simulation-process.md`　**局限**：`docs/04-known-issues-and-limits.md`

---

## 一句话结论（TL;DR）

- **时域**：可用 ngspice 收敛，但**免费未加密模型有物理边界**（VIN 脚电容 >~30µF 不启动）→ 只能用缩版 π 滤波、跳过 14ms 软启动 ⇒ **结果定性可用、定量不可信**（详见 01）。
- **传导 EMI（DM）**：含已确认的真实 π 滤波时，DM **远低于 CISPR 32 Class B**；去掉 π 滤波则明显抬高。**但真正的瓶颈很可能是 CM（共模）**，本模型尚未覆盖。
- **场**：已打通「PCB → Gmsh → Elmer 3D 涡流 → ParaView 出图」链路；几何/磁芯参数为简化+估值，**用于趋势判断，不用于定量**。

---

## 复现环境

- 计算主机：内网 `node249`（Ubuntu 24.04，12 核 / 15G，数据盘 RAID10 `/mnt/raid10`）
- 工具链：**ngspice 42**、**Elmer FEM 9.0**、**Gmsh 4.12**、**ParaView 5.11**、Python(numpy/scipy/matplotlib)
- 无头渲染：`LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a -s '-screen 0 800x600x24' pvbatch ...`

---

_生成于 2026-09-24。_
