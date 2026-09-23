# 01 · 仿真方法

本文写清**每个结果是怎么算出来的**：模型、拓扑、参数、求解器设置、假设与出处。
（时间线/迭代见 `02-simulation-process.md`；局限见 `04-known-issues-and-limits.md`。）

---

## 0. 三层仿真总览

| 层 | 目的 | 工具 | 产物 |
|---|---|---|---|
| A. 电路时域 | 纹波、开关波形、效率、启动 | ngspice（PSPICE 模型） | `results/01-stage1-transient/` |
| B. 频域传导 EMI | CISPR 32 DM 近似、滤波效果对比、软启动 | Python 线性节点分析 + FFT | `results/02-conducted-emi/` |
| C. 电磁场 FEA | 电感下铺铜的涡流/屏蔽/损耗 | Gmsh + Elmer 3D 谐波涡流 + ParaView | `results/03-field-eddy-copper/` |

> **关键前提**：SPICE 网表**不含 PCB 铜皮** ⇒ 三种铺铜方案（电感下铺铜 / 挖顶层 / 全挖）在 A、B 层结果**完全相同**，差异只能靠 C 层 FEA 体现。

---

## 1. 公共输入

### 1.1 拓扑（已从原理图逐脚核对）

```
24VIN ──F1──▶ PPHV ──[C17‖C62‖C63 = 3×10µF/50V]──┐
                        │D138(SMF24A TVS) │D1(SS36)
                        ▼
        PPHV ──L1(1µH, FXL0420-1R0-M)──▶ PPHV_OUT_FILTER
                                              │  [C73‖C74 = 2×10µF] + [C72=100nF]
                                              │  [R85=100mΩ 串 C16‖C64 = 20µF]  ← RC 阻尼
                                              ▼
                                          U21.VIN
        EN: R87(383k)/R88(35.7k) 分压 → U21.EN
        RT: R89 = 26.7kΩ → GND  (⇒ fsw ≈ 2.09e4/26.7 ≈ 783.6 kHz)
```

- 输入 π 滤波器 = **前级 3×10µF 银行 ‖ L1(1µH) ‖ 后级 2×10µF+100nF**。
- **已核实**：3×10µF **在 EMI 滤波电感 L1 之前**（见 `docs/03-input-filter-topology.md`）。
- 原理图详情：`docs/03-input-filter-topology.md`；勘误：`docs/04-correction-input-filter.md`。

### 1.2 器件型号（用于取寄生）

| 位号 | 型号 | 值 | 备注 |
|---|---|---|---|
| C17/C62/C63/C73/C74 | `CL31A106KBHNNNE` | 10µF | Samsung 1206 X7R 50V |
| C16/C64 | `JVJ50v10M5x5` | 10µF | 阻尼支路 |
| C72 | `CC0603KRX7R9BB104` | 100nF | 0603 X7R |
| C18 | `0603CG2R7C500NT` | 2.7pF | — |
| L1 | `FXL0420-1R0-M` | 1µH | EMI 滤波电感 |
| L2 | `ZEMS0650-150M` | 15µH | 功率电感（场仿真对象） |
| R85 | — | 100mΩ | 阻尼 |
| U21 | `TPS62933FDRLR` | — | F 版，带 SS 脚 |

### 1.3 模型文件

- 可用（未加密）：`TPS62933P_TRANS.lib`（`.SUBCKT TPS62933P_TRANS BST EN FB GND PG RT SW VIN PARAMS: STEADY_STATE=0`，**无 SS 脚**）
- 不可用：`TPS62933_TRANS.LIB`（`**$ENCRYPTED_LIB`）

> ⚠️ 板上是 **F 版（有 SS 脚）**，但能用的模型是 **P 版（无 SS 脚）** ⇒ 14ms 软启动**无法仿真**。

---

## 2. A 层：电路时域（ngspice）

**目标**：三负载工况（空载 ≈0A / 半载 1.5A / 满载 3.0A）的纹波、开关波形、效率。

**求解设置**（收敛调参结果）：
```
reltol=1e-3  abstol=1e-7  vntol=1e-4  gmin=1e-9  trtol=7  maxstep=2n
.tran 2n 80u uic      (STEADY_STATE=1)
```

**输入滤波的"缩版"**：免费模型在 **VIN 脚电容 >~30µF 时不启动**（内部 UVLO/使能不置位）⇒ 真实 π（>40µF）无法仿真，改用缩版 π `20µF–1µH–10µF`。

**结果**（详见 `results/01-stage1-transient/REPORT_v2.md`）：
| 工况 | Vout | Vpp | η |
|---|---|---|---|
| 空载 | 11.99 V | 8.4 mV | 49% |
| 半载 1.5A | 12.20 V | 347 mV | 72.6% |
| 满载 3.0A | 11.29 V | 424 mV | 71.8% |

（**定性可用、定量不可信**——未真正稳态、能量不守恒、缩版滤波。）

---

## 3. B 层：频域传导 EMI（CISPR 32 近似）

**思路**：开关电流的谐波 × 滤波网络传递阻抗 → LISN 口电压 → dBµV → 比 CISPR 限值。

```
V_LISN(k·fsw) = I_src(k·fsw) × |Zt(k·fsw)|,   dBµV = 20·log10(V/1µV)
```

- **Zt**：从"开关侧噪声源电流"到"LISN 口电压"的**传递阻抗**，用 2 节点线性节点分析求解：
  - 节点 A（LISN 口 / PPHV）：3×10µF 银行 + (1µF+50Ω) + 50µH
  - 节点 B（PPHV_OUT_FILTER = U21.VIN）：2×10µF + 100nF + (R85 + C16‖C64)
  - L1(1µH + DCR) 连接 A–B
- **源**：斩波输入电流的梯形波谐波（幅值 = 负载电流，D≈0.5，tr/tf≈15ns 估值），并用 `results_v2/E_*.txt` 数据的 FFT 作**交叉校验**。
- **限值**：CISPR 32 Class B / Class A QP（dBµV）。

**结果**（详见 `results/02-conducted-emi/REPORT_emi.md`）：
- DM 三工况**全过 Class B**（最差裕量 +34.4dB）——但**过于理想**（未含走线/回路电感、源为估值）。
- **CM（共模）粗估 ≈80dBµV，超 Class B 20–45dB** ⇒ **瓶颈很可能是 CM**。
- **软启动**：解析式 0→12V/14ms（非仿真），隐含 Iss = Css·Vref/Tss = 100nF·0.8V/14ms ≈ 5.71µA。

---

## 4. C 层：电磁场 FEA（电感下铺铜涡流）

**目标**：看电感 L2(15µH) 正下方顶层铺铜中的**感应涡流**，评估"铺铜 / 挖铜"对屏蔽与损耗的影响。

**链路**：PCB 工程（立创 epru）→ 提取几何 → Gmsh 网格 → Elmer 3D 谐波 A-V 涡流 → ParaView 出图。

- **求解器**：`CoilSolver`（闭合线圈电流源）+ `WhitneyAVHarmonicSolver`（复值，MUMPS）+ `MagnetoDynamicsCalcFields` + `ResultOutputSolver(VTU)`。
- **频率**：fsw = **784 kHz**（= 2πf 进 Angular Frequency）。
- **几何**（简化）**：空气域 40×40×16mm；铜皮 20×20mm（35µm）；磁芯 6.5×6.5×3mm；线圈 5.2×5.2mm 闭合方环。
- **关键坑**：**35µm 铜必须用 2D 壳单元**，不能切 3D 实体四面体（否则网格畸变、数值斑点）。
- **出图量**：`|J| = sqrt(Jre²+Jim²)`，铜平面切片，对数色标 1e2–1e6 A/m²。

**已知估值**：磁芯 **µr≈40**（无 datasheet，一体成型粉芯典型值）、线圈激励电流 6.0A（≈15匝×0.4A 纹波，估值）、铜厚 35µm（=1oz）。
→ **趋势可用，定量不保证**；真实铺铜几何 + 板框版本（fea3）见后续。

---

## 5. 复现

各结果目录内均含可重跑脚本；环境见根 `README.md`。命令示例：

```bash
# A 层
ngspice -b results/01-stage1-transient/netlists/E_full.cir
# B 层
cd results/02-conducted-emi && python3 emi_model.py
# C 层
gmsh -3 build.geo -o mesh.msh && ElmerGrid 14 2 mesh.msh -autoclean && ElmerSolver case.sif
LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a -s '-screen 0 800x600x24' pvbatch render_final.py
```
