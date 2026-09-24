# TPS62933 — 板级辐射发射仿真（openEMS）：FullCopper 探路 + 三版铺铜对比

日期：2026-09-24　目录：`/mnt/raid10/sim-work/tps62933/rad/`（镜像 `docs/projects/tps62933/rad/`）

> ⚠️ **绝对幅值不可信**，本报告只支持「**相对比较 / 热点定位 / 量级**」。详见 §5。
> 本报告由主会话在子代理交付后**核对数值并补写**（子代理因工具输出被系统截断，未落地 REPORT 文件；图表与 npz 已核验可复现）。

---

## 1. 目的
在已完成的**传导 EMI（DM/CM）**、**涡流场（fea2/fea3）**之外，补上**板级辐射发射**：
1. 探路跑通 FullCopper（B1）；
2. 对**三版铺铜**（FullCopper / TopCutout / DualCutout）做**同源同网格**的辐射对比，回答「电感 L2 下方挖铜对辐射是改善还是恶化、差多少」。

## 2. 工具链与模型
- **openEMS**（源码编译，CSXCAD v0.7.0 + openEMS v0.37.0 + Python 绑定），已用官方 patch-antenna 例校验（Dmax 4.79 dBi、效率 95%）。
- **FDTD**，30 MHz–1 GHz 一次宽带。
- 板：FR4 **1.51 mm εr=4.5**；铜全部用**零厚度薄层 PEC**（不做实心体网格）；**底层整板 GND**、**顶层 12 个信号 pour + 顶层 GND**（形状按导出的 pour 多边形/矩形重建）。
- **10 cm 竖直输入电缆**（直线导体代表）。

### 激励（关键）
- 建**显式输入热回路**：顶层 PEC 条 5×0.3 mm + 两个 GND 过孔，回路面积 ≈ 5×1.51 mm²。
- 回路顶层串**小集中源**（0.25 mm 缝隙、R=50 Ω、Gauss 脉冲），**避免源本身当天线**（此前 5 mm 宽 port 会变成偶极子，已修）。
- **归一化**：传递函数 `H = E_raw / |i_f_tot|`（dt 因先相消）→ `E_3m = H · |I_loop| / 3`；`I_loop` 取真实梯形脉冲串（**Ipk=4.3 A、805 kHz、D=0.5、tr=tf=5 ns**）的精确谐波幅值。
- **自洽校验**：smoke_test 里 `E_norm.max = 1.264e-11` vs 由 Prad/Dmax 推出的 `1.002e-11`，比值 **1.26（≈1 dB）** → 归一化方法可信。
- **SW 节点 dv/dt 源：本版未加。**

## 3. 三版几何（唯一差异）
`rad/compare/build_rad.py --ver=full|top|dual`：**同一模型/源/网格**，只改铜皮——
在 L2（电感）正下方 6.6×7.0 mm 设**禁布区**：
- `full` = FullCopper（不挖）
- `top`  = TopCutout（**顶层**挖）
- `dual` = DualCutout（**双面**挖）

日志确认：`top`→top_gnd pieces 43（比 full 的 38 多 5 块 = 挖出的窗口）；`dual`→bot pieces 4（底层也开窗）。网格三者相当（~273–283k cells），`if_tot@100 MHz` 几乎不变（1.788/1.783/1.784e-17）。

## 4. 结果（E@3 m，各角度最大值；×±2 谐波避开 D=0.5 偶次零点）

| 版本 | 30 MHz | 100 MHz | 300 MHz |
|---|---|---|---|
| FullCopper | 24.4 | 22.1 | 12.0 |
| TopCutout  | 27.0 | 24.6 | 14.2 |
| DualCutout | 27.6 | 25.2 | 15.1 |

**Δ vs FullCopper（dB）**

| 版本 | 30 MHz | 100 MHz | 300 MHz |
|---|---|---|---|
| TopCutout  | **+2.53** | **+2.50** | **+2.19** |
| DualCutout | **+3.13** | **+3.13** | **+3.04** |

**定量结论**：**电感 L2 下方挖铜【恶化】辐射**——顶层挖空约 **+2.2~2.5 dB**，双面挖空约 **+3.0~3.1 dB**（全频段趋势一致）。
物理含义：**L2 下方的铜面在起屏蔽 / 回流作用**，挖空 ≈ 开了一个辐射口径。
小环（磁偶极子）解析式 `E = ηk²AI/(4π)` 交叉校核**同量级**。

## 5. 假设与局限（务必阅读）
1. **绝对幅值不可信**：未建真实走线 / 器件寄生 / **共模路径** / 机壳地；电缆为保守直线，**预期偏低**。30 MHz 处 FullCopper 最坏 ~24–26 dBµV/m，低于 CISPR 32 B@3 m（~40 dBµV/m）约 14 dB，但**此判读不可当合格结论**。
2. **FDTD 未收敛到 -40 dB 判据**：三版均在 `max timesteps` 到顶时结束（能量仍在 -0dB~-5dB 波动）。因三版设置**完全一致**，相对比较仍成立；但绝对/高频尾段可能未收敛。
3. **GND pour 形状为近似**：板框 − 信号 pour bbox + 外扩 clearance；精确 GND 多边形未导出。信号 pour 只用矩形/多边形近似。
4. **热回路面积/位置为估值**；**D=0.5 精确占空比**导致偶次谐波深零点（已用 ±2 谐波包络规避）。
5. 方向图近各向同性（Dmax≈0 dBi，小电尺寸 + 无真实地参考）。
6. 三版对比用**一致的 fast 网格（1.5 mm）**；0.8 mm 的 full 分辨率 FullCopper 单跑被中止（避免超时）。
7. **SW dv/dt 源、共模激励、真实电缆/参考地**均未加 → 下一步（B2）。

## 6. 产物（绝对路径，node249）
- `rad/b1/board_b1.py`、`post_b1.py`、`run` 日志 `fast_run.log`、`farfield_b1.npz`
- `rad/b1/fig_farfield_3d.png`、`rad/b1/fig_farfield_vs_cispr.png`
- `rad/compare/build_rad.py`、`post_compare.py`、`compare.npz`
- `rad/compare/ver_{full,top,dual}/{port.npz,nf2ff.npz,sim/}`
- `rad/compare/fig_rad_compare_freq.png`、`fig_rad_compare_hotspot.png`
- `rad/compare/run_{full,top,dual}.log`
- 本报告 `rad/REPORT_rad.md`

## 7. 下一步（B2）
精确导出三版 pour 多边形；加 SW dv/dt + 共模激励与真实电缆/参考地；器件寄生、真实过孔走线；网格收敛性验证与实测校准。
