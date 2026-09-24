# TPS62933 辐射三维可视化 (rad3d) — 方法 / 产物 / 局限

> ⚠️ **注意（2026-09-24）**：本目录可视化**衍生自旧前提**的 `rad/` 远场（无参考地平面/无电缆/无 CM 通路/无上游参数）。
> 因此**数值不可参考**，仅作**三维可视化演示**（板模型 + 方向图 + 传播动画）。当前有效辐射结果见 `results/09-radiated-b2/`。

本目录是 `rad/` 的**可视化衍生**（只读引用 `rad/`，未改动任何 `rad/` 文件）。
目标：把已有的板级 openEMS 辐射结果做成**立体直观图** + 一个**场传播动画**，并结合一块
**从板几何重建的 STEP 3D 板模型**。

---

## 0. 产物清单（绝对路径，均在 node249）

| 文件 | 内容 |
|---|---|
| `/mnt/raid10/sim-work/tps62933/rad3d/board.step` | **真实 STEP 三维板模型**（AP214，由 gmsh/OCC 从 `rad/geom/` 板几何生成）+ 器件方块 |
| `/mnt/raid10/sim-work/tps62933/rad3d/fig_pattern_3d.png` | 30 / 100 / 300 MHz **三维远场方向图**（三联，r ∝ 归一化场，板模型置原点） |
| `/mnt/raid10/sim-work/tps62933/rad3d/propagation.gif` (+ `.mp4`) | **近场 \|E\| 时间快照传播动画**（31 帧，0–1.88 ns，侧视 + 俯视双切面） |
| `/mnt/raid10/sim-work/tps62933/rad3d/fig_board3d.png` | **STEP 板模型 + 近场 \|E\| 切面**叠加三维视图（ParaView 渲染） |
| `.../field_slice.vtk` | 动画末帧 \|E\| 切面（VTK 结构化网格，供 ParaView） |
| `.../board.geo`, `board.stl` | gmsh 几何源与 STL（ParaView 无 STEP reader，用 STL 渲染） |
| `.../prop/sim/*.h5` | 原始 openEMS 时域 E 场 dump（`E_sliceV.h5`, `E_sliceH.h5`） |
| 脚本 | `gen_board_geo.py` `make_pattern3d.py` `build_prop.py` `make_gif.py` `export_field_vtk.py` `pv_board3d.py` |

---

## 1. STEP 板模型 `board.step`

- 工具：**gmsh 4.12 OpenCASCADE**（`SetFactory("OpenCASCADE")`），脚本 `gen_board_geo.py`
  → `board.geo` → `gmsh board.geo -0 -o board.step`（AP214，480 KB，多实体 Compound）。
- 几何来源：`rad/geom/board_outline.csv` / `components.csv`。
  - **板框**：X ∈ [−6.223, 20.277] mm，Y ∈ [−32.974, 17.526] mm（= 26.50 × 50.50 mm，**精确**），
    沿 Z 拉伸 FR4 厚度 **1.51 mm**（**精确**），板中心 z ∈ [−0.755, +0.755] mm。
  - **器件**：`components.csv` 每个器件按中心 (x,y) 与 `angle_deg` 旋转生成方块，置于顶层。
- **L2（ZEMS0650-150M）**：footprint 取 **6.6 × 7.0 mm** —— 来源是 `rad/b1/board_b1.py`
  中为 L2 定义的 **keep-out 区域**（精确，非猜测），高度 3.0 mm 为**估值**。
- ⚠️ **器件尺寸为估值**：`components.csv` 只给中心坐标+封装型号，无 bbox。各器件方块按**封装代号**
  推出的行业标准尺寸（EIA 0603/0805/1206/1210、SOT-563、SMA、2410 fuse…），**高度多为估值**。
  完整对照表见 `gen_board_geo.py` 的 `COMP_SIZE`。**这些尺寸非厂家原始 STEP 数据。**
- ⚠️ **本 STEP 是从导出的板几何重建的，不是厂家原始机械 STEP。**

---

## 2. 三维远场方向图 `fig_pattern_3d.png`

- 输入：`rad/b1/nf2ff_b1.npz`（E_norm[36,37,401], freq, theta, phi, Dmax）、`rad/b1/port_b1.npz`（if_tot）。
- **归一化口径沿用 `rad/b1/post_b1.py`（未改）**：
  - `H(f,θ,φ) = E_norm / |i_f_tot(f)|`　[Ohm/m]（dt 因子相消）
  - `E_3m = H · |I_loop(f)| / 3`　[V/m]
  - `I_loop` = 输入热环电流实波形的精确谐波谱：**Ipk=4.3 A, fsw=805 kHz, D=0.5, tr=tf=5 ns**。
- 显示：球面半径 `r = E/E_max`（**归一化线性场**，只反映**形状**），面颜色 = 相对最大值的 dB。
  板模型以示意矩形置于原点（**非等比例**，仅示意朝向）。
- 本次取最近谐波结果：

| 目标 | 实际谐波 | E_max @3m |
|---|---|---|
| 30 MHz | 37×805k = 29.8 MHz | **29.1 dBµV/m** |
| 100 MHz | 124×805k = 99.8 MHz | **26.5 dBµV/m** |
| 300 MHz | 373×805k = 300.3 MHz | **−31.2 dBµV/m** |

  （300 MHz 偏低是因为 tr=5 ns 的梯形谱在 1/(π·tr)≈64 MHz 后按 −20 dB/dec 滚降；与 `rad/REPORT_rad.md`
  中 E 随频率下降的趋势一致。三版/限值线口径见 `rad/compare/post_compare.py`。）
- **方向图近各向同性（Dmax≈0 dBi）** —— 与 `rad/REPORT_rad.md` 结论一致（小电尺寸 + 无真实地参考）。

---

## 3. 近场传播动画 `propagation.gif`

做法（**重新跑了一次短时 FDTD**，脚本 `build_prop.py`）：

- 模型 = `rad/b1/board_b1.py` 的 **b1「fast」粗网格版**（薄层铜用 **PEC 薄层**，不切体网格），
  频率设置同 `rad/b1`（`SetGaussExcite(f0=515 MHz, fc=485 MHz)`，即 30 MHz–1 GHz）。
- 网格 **75×89×43 = 273 504 单元**；**dt ≈ 1.93×10⁻¹⁴ s (19.3 fs)**（与 b1 fast 一致）。
- `NrTS = 100 000` → 模拟 **≈1.88 ns**；`OverSampling = 8`
  → dump 间隔 = Nyquist/8 = **62.5 ps** → **31 帧**。
- dump：**时域 E 场**（`dump_type=0`, `file_type=1`, `dump_mode=2`）在两张切面上：
  - `E_sliceV`：x–z 平面（y = −6.3 mm，穿过源）→ **侧视**
  - `E_sliceH`：x–y 平面（z = +3 mm，板面上方）→ **俯视**
- 渲染：matplotlib 出 PNG 序列 → ffmpeg 合成 `propagation.gif`（7 fps）与 `.mp4`。
  配色为 **turbo 亮基色带 + dB 标度（相对全局最大，量程 −85…0 dB）**，
  **未使用 log 0.1–3e7 那种量程**（会压成全黑）。
- 运行开销：219 s @ 131 MCells/s（12 线程）。

---

## 4. 板模型 + 场叠加 `fig_board3d.png`

- ParaView 5.11 无头渲染（`LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a -s '-screen 0 1200x900x24' pvbatch pv_board3d.py`）。
- ParaView **无 STEP reader** → 用 `board.stl`（gmsh 从同一 OCC 几何导出）显示板模型；
  叠加动画末帧（t=1.875 ns）的近场 `field_slice.vtk`（`E_dB` 标量，Turbo 色带）。

---

## 5. 局限（**必读**）

1. **绝对幅值不可信。** 沿用 b1 归一化口径（`H=E/|i_f_tot|`, `E_3m=H·|I_loop|/3`）只能给**相对/形状**
   信息；绝对电平依赖 `I_loop` 与模型理想化程度，**不可用于合规判定**。方向图 r 已进一步归一化。
2. **仅板级 + 一根 10 cm 简化电缆**（b1 模型），**没有真实共模路径 / 线束 / 机壳 / 地平面参考**。
   现实中辐射常由电缆共模主导，本可视化**不含**这些路径。
3. **传播动画用粗网格 + 短时**：只跑到 ≈1.88 ns；因 openEMS 该 Gaussian 激励本身长达数 ns
   （远场运行里 openEMS 请求 5.9 ns 脉冲），模拟结束时**源仍在驱动**，动画显示的是
   **近场随源建立并向外扩展**的过程，而非一个已分离的辐射波包。
4. **该结构是弱辐射体（Dmax≈0 dBi）**：近场能量集中在开关环附近，`|E|` 在离源几 mm 外即比峰值低数十 dB。
   故动画中场的可见范围很小（−85 dB 量程下 ~10–30 mm）。这是**物理事实**，不是渲染问题。
5. **所有器件尺寸（尤其高度）为估值**（封装代号推得）；L2 footprint 取自 b1 keep-out（精确），高度为估值。
6. **STEP 是重建模型，非厂家原始 STEP。**
7. 频率/谐波/CISPR 32 Class B@3m 限值线口径与 `rad/compare/post_compare.py` 一致，本目录**未做合规对比**，
   只做可视化。

---

*生成于 node249；未修改 `rad/`、`fea*/`、`emi*/`、`results*/epro/`；未 push、未外发。*
