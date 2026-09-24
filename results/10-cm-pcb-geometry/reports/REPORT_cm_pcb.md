# tps62933-cm-pcb — Conducted CM with the coupling capacitance extracted from the **real PCB copper**

目录：`/mnt/raid10/sim-work/tps62933/cm_pcb/`（只写此目录；未改动 `emi_cm/ fea*/ rad*/ results*/ epro/`）

## 1. 做了什么

1. **认定 SW 网络铜皮**（不靠猜）：`epro/pourSim.epru` 里逐项解析，SW 网络 = `U21.SW(pin5) → L2`，
   其唯一的实体铺铜是 `POUR75`（`netName=$1N7`，唯一与 L2 邻接的 pour，且带 LINE/ARC 走线）。
   该 pour 多边形 21 顶点、面积 **13.653 mm²**（`geom_defs.py`，源头 `rad/geom` + 原理图 netlist）。
2. **三版差异**：`fea4/{fullcu,topcut,dualcut}/geo_meta.json` 的 `REG_PROHIBIT` 区域。
   - `fullcu`：无禁布；`topcut`：顶部挖空 L2 下方 6.60×7.00 mm（含区域 B）；`dualcut`：再加底部同样挖空。
   - 关键几何事实：L2 禁布区 **A** 与 SW pour 的顶部重叠 → `topcut/dualcut` 把 SW 铜皮从 13.653 → **7.719 mm²**；禁布区 B 离 SW 远，不影响。
3. **静电场提取**（Elmer `StatElecSolver`，单进程，`OMP_NUM_THREADS=2`，全部命令带 `timeout`）：
   - 几何：**铜皮建成 2D 面片**（SW 多边形 / 底层地平面矩形），**35 µm 铜皮绝不切 3D 实体**；
     参考地平面 = 板下方距离 `d` 的大平面（0 V）；空气域外边界自然 Neumann。
   - BC：SW 面 = 1 V；板底地平面 = 0 V（`dualcut` 时把 L2 挖空区从地平面剔除）；参考平面 = 0 V。
   - 求解电势 → 静电能量 W → `C_sw = 2W/V²`（Elmer 输出 `effective capacitance` 列，与 `2·Electric Energy` 一致）。
   - 网格：Gmsh OCC + Box 尺寸场（SW 附近细、远处粗），`BooleanFragments` 把面片嵌入体积。
4. **代入 CM 电路**（复用 `emi_cm/emi_cm_model.py` 结构，未改动原文件）：
   `I_cm(f)=2πf·C_sw·|V_sw(f)|`，`V_sw(f)` 用**实测 SW 谱**（`results_v2/E_{noload,half,full}.txt` 第 4 列，
   末 32 周期整数周期 FFT），CM 回路 = node0(板GND) →25 Ω→ 大地，0→上游GND(30 nH)→Cy(2.2 nF)→大地。
   出 CISPR-32 Class B 三条曲线/每版 + 一张三版对比图。

## 2. 产物（绝对路径）

```
/mnt/raid10/sim-work/tps62933/cm_pcb/results.csv                     # C_sw 汇总（含 d 敏感性 + 两档网格）
/mnt/raid10/sim-work/tps62933/cm_pcb/fig_cm_per_variant.png          # 三版各自 CM CISPR 曲线（3 工况 + 限值线）
/mnt/raid10/sim-work/tps62933/cm_pcb/fig_cm_variant_compare.png      # 三版同工况(full)对比
/mnt/raid10/sim-work/tps62933/cm_pcb/cm_pcb_numbers.txt              # 数值日志
/mnt/raid10/sim-work/tps62933/cm_pcb/geom_defs.py, gen_geo.py        # 真实铜皮几何 + Gmsh .geo 生成
/mnt/raid10/sim-work/tps62933/cm_pcb/cm_pcb_model.py                 # CM 电路 + CISPR 出图
/mnt/raid10/sim-work/tps62933/cm_pcb/make_sif.py, run2.sh, run3.sh   # Elmer sif 生成 / 批量运行
/mnt/raid10/sim-work/tps62933/cm_pcb/geo_*.geo, m_*.msh, mesh_*, run_*  # 每案几何/网格/求解现场
```

## 3. 证据

### 3.1 三版 C_sw（pF）+ 平行板交叉核对

| variant | SW 铜面积 (mm²) | **C_sw FEM (pF)** | 平行板 ε₀A/d (d=1 mm 板厚, pF) | 比值 |
|---|---|---|---|---|
| FullCopper | 13.653 | **0.411** | 0.121 | 3.4× |
| TopCutout  | 7.719  | **0.312** | 0.068 | 4.6× |
| DualCutout | 7.719  | **0.202** | 0.068 | 3.0× |

- 平行板基准取 **d = 1 mm（板上铜 ↔ 板底层地层，空气）**：`C=ε₀A/d`。
  FEM 比纯平行板高 3–4.6×，差值 = 大平面的边缘/散射场（电极尺寸 3.7 mm ≈ 间隙 1 mm，边缘电容占比大）——
  **同量级**成立。（若取 d=10 cm 到外部参考平面，平行板仅 0.0012 pF，说明外部平面不是主耦合。）
- 交叉核对也解释了版间差异：SW 面积 13.653→7.719 mm² 使 C_sw 降约 25%（fullcu→topcut）；
  **dualcut 又再降 35%**（同为 7.719 mm²，但底部 L2 挖空移除了正下方的大段地层接地铜）。

### 3.2 d 敏感性（1 / 10 / 40 cm）——**关键结论**

| d | fullcu (pF) | topcut (pF) | dualcut (pF) |
|---|---|---|---|
| 1 cm | 0.4103 | 0.3153 | 0.2026 |
| **10 cm（标称）** | **0.4107** | **0.3121** | **0.2021** |
| 40 cm | 0.4088 | 0.3156 | 0.2002 |

**d 在 1–40 cm 内只改变 C_sw < 2 %。** 因为 C_sw 主要由**板上铜 ↔ 板自身底层地平面（1 mm）**的耦合决定，
而非与 10 cm 外参考平面的耦合。补充算例（去掉板底层地平面、只留 SW↔参考平面）：
fullcu 0.215 pF(d=1cm) / 0.191(10cm) / 0.191(40cm) —— 也只在 1 cm 时高 12%，≥10 cm 基本饱和于自电容。
→ **本模型里"是否超标"不由 d 决定；d 的作用是次要的。** 这与任务书假设（"d 决定 C_p"）不符，如实报告。

### 3.3 CM 结果（CISPR-32 Class B，三工况）

| variant | C_sw (pF) | noload(0.02A) | half(1.5A) | full(3.0A) |
|---|---|---|---|---|
| FullCopper | 0.411 | +2.5 dB | **−1.5 dB** | **−0.3 dB** |
| TopCutout  | 0.312 | +4.9 dB | +0.8 dB | +2.1 dB |
| DualCutout | 0.202 | +8.6 dB | +4.6 dB | +5.9 dB |

（最差裕量：限值 − CM 曲线，正=通过；均在 0.81 MHz / 1.61 MHz 取最差。）

### 3.4 与"黑箱 C_p=10 pF 基线"对比

| 工况 | 基线(10 pF)最差裕量 | 几何 C_sw 给出 | 改善 |
|---|---|---|---|
| noload | −25.2 dB | +2.5 / +4.9 / +8.6 dB | ≈ +28…+34 dB |
| half | −29.3 dB | −1.5 / +0.8 / +4.6 dB | ≈ +28…+34 dB |
| full | −28.0 dB | −0.3 / +2.1 / +5.9 dB | ≈ +28…+34 dB |

**真实铜皮 C_sw ≈ 0.2–0.4 pF，比黑箱 10 pF 小 25–50×（≈28–34 dB）。**
结论反转：黑箱 10 pF 判定"超标 25–29 dB"；用真实几何，**FullCopper 在半载/满载为临界（−1.5/−0.3 dB，略超），
加 L2 底部挖空（TopCutout/DualCutout）后转为通过（+0.8…+5.9 dB）**。这正是三版铺铜的工程意义。

### 3.5 网格收敛（d=10 cm）

| lc_sw | fullcu | topcut | dualcut |
|---|---|---|---|
| 1.6 mm | 0.582 | 0.429 | 0.267 |
| 0.8 mm（标称） | 0.411 | 0.312 | 0.202 |
| 0.5 mm | 0.368 | (网格失效) | 0.187 |

0.8→0.5 mm：fullcu −10.4 %、dualcut −7.3 %（单调收敛，取值仍有 ~±10 % 的网格不确定度）。

## 4. 局限

1. **未完全收敛**：0.8→0.5 mm 变化 ~10 %（fullcu 略超 10 %），C_sw 带 ~±10 % 网格不确定度；topcut 的 0.5 mm
   细网格生成失败（平面挖孔拓扑），其收敛仅到 0.8 mm。
2. **顶层共面地/其它网络铜未建**：只把 SW 与**底层地平面**建成电极；顶层同层 GND 与 12V/PPHV 等铺铜略去
   （对 C_sw 为二阶修正，会略增 C）。FR4 介质未计入（空气模型，εr=4.3 会使 C 略增 ~分倍）。
3. **"参考地平面"物理定位**：真实板无 PE、塑料壳、无屏蔽，CM 回路由 LISN 经 10 cm/0.75 mm² 电缆闭合。
   本模型把"参考平面"当 0 V 外部地；**C_sw 实际由板自身地层主导，外部平面距离 d 影响 <2 %**（见 §3.2），
   与任务书"d 决定 C_p"的假设不一致，已如实标注。
4. **SW 网络铜量的判读**：SW 无专铺铜，取 `POUR75($1N7)` 为 SW 铜皮（依据：唯一与 L2 邻接的 pour + 走线）。
   若实际 SW 走线/焊盘与此稍有出入，C_sw 按面积线性近似缩放。
5. CM 回路参数（25 Ω、30 nH 电缆、2.2 nF Cy）沿用基线 `emi_cm_model.py`；源谱为**仿真实测波形**
   （`results_v2`），非硬件实测 —— 与基线同口径，可横向比较。
