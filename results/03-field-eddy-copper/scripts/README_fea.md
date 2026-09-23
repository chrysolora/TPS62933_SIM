# TPS62933 L2 电感下铺铜 — 3D 涡流场仿真（链路验证）

> 阶段2 · 单例 + 粗网格链路验证。**不做**三方案完整对比 / 参数扫描 / 高精度。
> 目标：验证「PCB 几何 → Gmsh 网格 → Elmer 3D 涡流 → ParaView 出图」链路可用。

## 一、结论

**链路跑通 ✅**。用 Elmer 9.0 的 `WhitneyAVHarmonicSolver`（3D 谐波 A-V 涡流）+ `CoilSolver`（闭合线圈电流源）
成功求解了「屏蔽一体成型电感 L2 在其正下方顶层铜皮中感应的 497 kHz 涡流场」，
导出铜皮上的 `current density re/im`、`joule heating` 等场量，并用
`xvfb-run + pvbatch` 无头渲染出 PNG（`fea_field_top.png` / `fea_field_3d.png`）。

关键中间结果（铜皮单元场量范围，粗网格）：
- `current density re` ≈ ±9.9e6 A/m²，`current density im` ≈ ±1.3e7 A/m² → 存在明显涡流
- `joule heating` ≈ 0 … 2.6e9 W/m³（谐波实部，含符号）
- 铜皮 497 kHz 趋肤深度 δ = √(2/(ωμσ)) ≈ **94 µm**

## 二、从 PCB 提取的几何与叠层

原始工程：`/mnt/raid10/sim-work/tps62933/epro/pourSim.epru`（立创EDA Pro epru，单文件含全部文档）。
解析脚本：`parse_epru.py` / `extract_l2.py` / `geom_pour.py`（见本目录）。

- **叠层（LAYER_PHYS，3 套 PCB 完全一致）**：**2 层板**（顶层 TOP + 底层 BOTTOM），
  - 铜层：TOP(zIndex 1) / BOTTOM(zIndex 2)，`thickness=0`（未显式给出 → **铜厚按 1 oz = 35 µm 估值**）
  - 介质：FR4 core `59.449 mil ≈ 1.510 mm`，εr = **4.5**
  - 阻焊：top `0.394 mil`，εr = 3.3
- **L2 器件**：`ZEMS0650-150M`，DeviceName 确认，封装 6.5×6.5 mm。
  - 坐标：`x=245, y=-660`（mil 坐标，TOP 层 layerId=1），旋转 90°
  - 主电感 **L1 = FXL0420**（另一颗，位于 645,-1250），本任务只关注 L2
- **L2 正下方铺铜**：顶层 L2 所在区域被 3 个铜皮轮廓覆盖（`POUR` net）：
  - `POUR24 / $1N17`、`POUR2 / $1N22`、`POUR3 / PPHV`
- **三套 PCB 的差异**：`POUR` 轮廓（网名/层）三例一致；差异在 `POURED` 填充结果多边形
  （填充区域数 79 vs 80，以及局部少量多边形不同）——即「电感下是否保留铜皮填充」。
  POURED 填充数据为局部坐标、且未带 net/层标签，未逐点多边形对齐（属本次不做完整对比的边界内）。
  → 本次链路验证按 **「电感下铺铜」** 这一例建模：L2 正下方顶层铜皮**保留**为整片平面。

## 三、磁芯参数（µr）

- **未获取到 ZEMS0650-150M 官方 datasheet**（无外网检索接口，未查到该型号磁芯 µr / tanδ）。
- **估值（明确标注为估值）**：
  - 一体成型屏蔽电感用金属/合金粉芯（distributed gap），有效相对磁导率 **µr ≈ 40**
    （此类 6.5 mm 尺寸、~15 µH 产品典型 µr 区间 20–60，取中值 40）。
  - 本次只关心**铜皮涡流损耗**，磁芯损耗暂不建模（磁芯 `Electric Conductivity = 0`）。
  - **依据**：一体成型粉芯为降低饱和/损耗采用分布式气隙，有效 µr 远低于铁氧体（数百~数千），
    典型落在数十量级。**请以实测/datasheet 为准。**

## 四、仿真模型（简化 + 粗网格）

Gmsh 4.12.1 建模（`build.geo`，单位 m），OpenCASCADE，`BooleanFragments` 保证体共形：

| 体 | 几何 | 材料 |
|---|---|---|
| Air | 20×20×8 mm 域 | µr=1, σ=0 |
| Copper | 16×16 mm × 35 µm，z∈[-35,0]µm（L2 正下方顶层铺铜） | µr=1, **σ=5.8e7** |
| Core | 6.5×6.5×3.0 mm 屏蔽磁芯块 | **µr=40(估值)**, σ=0 |
| Coil | 5.2×5.2 方环、截面 0.6×0.5 mm 的**闭合回路**（置于磁芯内，代表绕组通流） | µr=1, σ=1(仅用于 CoilSolver) |

- 网格：Gmsh 粗网格 → 1701 节点 / 7507 四面体（`ElmerGrid 14 2 mesh.msh -autoclean`）
- 激励：`CoilSolver` 闭合线圈，`Desired Coil Current = 6.0 A`（= 估算安匝 N·I，N≈15 匝 × 纹波幅值≈0.4 A，**估值**）
- 频率：`Angular Frequency = 3.1227e6`（= 2π × 497 kHz）
- 边界：外边界 `AV = 0`（远场）
- 求解器：`CoilSolver` → `WhitneyAVHarmonicSolver`（复值，MUMPS 直接法）→ `MagnetoDynamicsCalcFields` → `ResultOutputSolver(VTU, ASCII)`

## 五、复现步骤

```bash
cd /mnt/raid10/sim-work/tps62933/fea/case1
gmsh -3 build.geo -o mesh.msh                 # 生成网格
ElmerGrid 14 2 mesh.msh -autoclean           # 转 Elmer 网格
ElmerSolver case.sif                         # 求解（输出 mesh/case_t0001.vtu）
LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a -s '-screen 0 800x600x24' pvbatch render.py   # 出图
```

## 六、关键坑位（供后续复用）

1. **闭合线圈必须在 Elmer `Component` 段定义**（`Master Bodies` / `Coil Closed` / `Coil Normal` / `Desired Coil Current` / `Coil Type = String "Stranded"`），
   只在 Solver 段写 `Coil Closed` 会被忽略（CoilSolver 会把整个求解域当成线圈）。
2. `CoilSolver` 需 `Calculate Elemental Fields = Logical True`，否则 Whitney 报
   `Elemental current requested but not found: CoilCurrent e`。
3. `WhitneyAVHarmonicSolver` 需 `Use Elemental CoilCurrent = Logical True` 才能用 CoilSolver 的预计算电流密度。
4. 谐波复值系统用 **`Linear System Solver = Direct` + `Direct Method = MUMPS`** 才收敛
   （BiCGStab+ILU1 会 "Too many iterations"）。
5. Elmer 默认 VTU 为 `encoding="raw"` 二进制，ParaView 读失败 → 加 **`Ascii Output = Logical True`**。
6. 线圈材料给 σ=1（非 0），否则 CoilSolver 矩阵奇异出 NaN。
7. 无头渲染必须 `xvfb-run + LIBGL_ALWAYS_SOFTWARE=1 pvbatch`。

## 七、局限（务必知悉）

- **单例 + 粗网格**：仅验证链路，**未做**「铺铜 / 挖顶层 / 全挖」三方案对比，**未做**网格无关性/收敛性验证。
- **几何简化**：L2 用「磁芯块 + 单匝方环线圈」等效，未还原真实多层绕组与焊盘；
  铺铜按整片平面近似（未逐点多边形切割 POURED 填充）。
- **磁芯参数为估值**（µr=40、无磁芯损耗），铜厚按 1 oz 估值；线圈安匝为估值。
- 涡流损耗绝对值**不可用于定量结论**，仅证明场分布与出图链路成立。
- 三套 PCB 的「铜皮挖除」差异未在几何上落实（POURED 为局部坐标未对齐）。

## 八、文件清单

```
fea/
├── README.md               本文件
├── fea_field_top.png       ★ 交付图：铜皮平面俯视切片，按 joule heating 着色
├── fea_field_3d.png        3D 视角切片
├── case.sif                Elmer 求解输入
├── build.geo               Gmsh 几何
├── render.py               ParaView/pvbatch 出图脚本
├── case1/                  完整运行目录（mesh/、run.log、*.vtu、*.png）
├── parse_epru.py / extract_l2.py / geom_pour.py / ...  PCB 解析脚本
└── pourSim.epru.bak        原始工程备份
```
