# 板级电磁**辐射**仿真 — 准备阶段报告 (`rad/REPORT_rad_prep.md`)

主机：`node249`（Ubuntu 24.04, 12 core / 15 GB, `/mnt/raid10/sim-work/tps62933/`）
时间盒：1.5 h（实际约 2 h，见"局限")
范围：**只做安装 + 验证 + 几何导出**，未做正式板级远场仿真。

---

## 1. 装了什么 / 怎么装

**结论：Ubuntu 24.04 官方源没有 openEMS；`thliebig/openEMS-ppa` 在本机被屏蔽（HTTP 403）。
最终采用【源码构建】方式，从 github.com/thliebig 编译 CSXCAD + openEMS（含 python3 绑定）。**

| 组件 | 版本 / 来源 | 安装位置 |
|---|---|---|
| CSXCAD (C++ lib) | git `bd2c133` → v0.7.0（源码构建） | `/usr/local/lib/libCSXCAD.so`, `/usr/local/include/CSXCAD` |
| openEMS 求解器 | git `65f8771` → v0.37.0（源码构建） | `/usr/local/bin/openEMS` |
| CSXCAD python | pip wheel `CSXCAD-0.7.0...` | `~/.local/lib/python3.12/site-packages/CSXCAD` |
| openEMS python | pip wheel `openems-0.37.0...`（含 `_nf2ff`） | `~/.local/lib/python3.12/site-packages/openEMS` |
| h5py | pip 3.16.0 | 同上 |

构建期使用的库（运行时 `openEMS` 自报）：`HDF5 1.10.10`、`tinyxml 2.6.2`、`fparser`、`boost 1.83`、`VTK 9.1.0`、`CGAL`。

### 安装步骤（可复现）
```bash
# 1) 基础依赖
sudo apt-get install -y cmake g++ gfortran git swig python3-dev python3-pip python3-venv \
  libhdf5-dev libtinyxml-dev libfparser-dev libboost-all-dev libcgal-dev \
  libjsoncpp-dev libpng-dev libdouble-conversion-dev liblz4-dev liblzma-dev \
  libtiff-dev libtbb-dev libutfcpp-dev libexpat1-dev libeigen3-dev libglew-dev libgl2ps-dev
# 2) 源码
git clone --depth 1 https://github.com/thliebig/CSXCAD.git
git clone --depth 1 https://github.com/thliebig/openEMS.git
# 3) 编译安装 CSXCAD -> /usr/local
cd CSXCAD && mkdir cmb && cd cmb && cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release
make -j12 && sudo make install && sudo ldconfig
# 4) 编译安装 openEMS -> /usr/local
cd ../../openEMS && mkdir cmb && cd cmb && cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release
make -j12 && sudo make install && sudo ldconfig
# 5) python 绑定
cd ../../CSXCAD/python && CSXCAD_INSTALL_PATH=/usr/local python3 -m pip install --break-system-packages .
cd ../../openEMS/python && CSXCAD_INSTALL_PATH=/usr/local OPENEMS_INSTALL_PATH=/usr/local \
   python3 -m pip install --break-system-packages .
```

### 两个必须记录的坑（均已绕过）
1. **VTK 依赖冲突**：`libvtk9-dev` 与已安装的 `python3-paraview` 声明冲突（会连带卸载 ParaView）。
   由于 `CSXCAD`/`openEMS` 的 CMake 里 `find_package(VTK REQUIRED)`，无 VTK 无法配置。
   **绕过**：用 `apt-get download libvtk9-dev` 取 .deb，**`dpkg -x` 就地解包到 `/`**（不改 apt 数据库、可删可回滚），
   并把 `/usr/lib/x86_64-linux-gnu/cmake/vtk-9.1/vtk-config.cmake` 里的 `WRAP_PYTHON` 由 ON 改 OFF（原文件已备份 `.bak`）。
2. **缺 VTK 第三方模块**：`find_package(VTK COMPONENTS IOGeometry IOPLY)` 需要 jsoncpp/png/double-conversion/lz4/lzma/tiff/tbb/utfcpp 等 dev 包（上面第 1 步已一并装）。

> 注：`AppCSXCAD` GUI 未安装成功（需完整 Qt+VTK GUI 栈），但**命令行求解器 + NF2FF 已完全可用**，对无头仿真无影响。

---

## 2. 官方小例子跑通（证据在 `rad/smoke_test/`）

### 2.1 官方教程 `Simple_Patch_Antenna.py`（推荐基线）
- 运行：`MPLBACKEND=Agg python3 Simple_Patch_Antenna.py`
- **求解器跑通**：`FDTD 49x47x45 = 103635 cells`，`Speed 204.5 MCells/s`，`-52.89 dB` 收敛，仿真约 8 s。
- **远场出结果**（NF2FF）：
  - `f_res = 2.430 GHz`
  - `Dmax = 4.792 dBi`（贴片天线量级正确 ✓）
  - `Prad/Pin = 95.3%`（辐射效率物理合理 ✓）
- 出图：`rad/smoke_test/pattern_farfield.png`（远场方向图，θ 扫描，已完成渲染）。

### 2.2 自定义偶极子 `dipole_nf2ff.py`（2.4 GHz 半波偶极 + NF2FF）
- `68921→103243 cells`，收敛，`Prad/Pin` 一致，`Dmax ≈ 1.66 dBi`。
- 出图：`rad/smoke_test/dipole_sim/pattern_2p4GHz.png`。

### ⚠️ 关于 Prad 绝对值的单位约定（重要，正式仿真前须处理）
官方例子里 `nf2ff_res.Prad[0] = 5.56e-25 W`、`port.uf_inc ≈ 7.64e-12 V`（≈ 时间步长 dt）。
openEMS 的端口/远场频域谱**按 DFT 约定未归一化**；`Pin = 0.5·|uf_inc|²/R = 5.83e-25 W` 与之自洽，
故 **`Prad/Pin=95%` 是正确的**，但绝对值需按激励谱 **除以 `|uf_inc(f)|²`** 才是真实瓦特
（`5.56e-25 / (7.64e-12)² ≈ 9.5e-3 W`）。正式板级仿真必须做这一步归一化，否则报告 W 会错 20 个数量级。

---

## 3. 起始几何导出（`rad/geom/`）

复用 `fea3/parse_lib.py`（`load_pcbs` / `parse_path`），脚本：`rad/geom/export_geom.py`。
坐标系：**米，原点 = L2 中心 (245, -660) mil**，与 fea3 一致。**未修改 epro/ 任何文件（只读）。**

| 产物 | 内容 | 关键数值 |
|---|---|---|
| `board_outline.csv` | 板框多边形 (mil + m) | **板框 26.50 × 50.50 mm ✓（与背景一致）** |
| `pours_top.csv` | 顶层 13 个铜 Pour（name/net/bbox） | 含 GND? 顶层无 GND 覆盖全区；`POUR4 net=GND x[-5,1220] y[-2010,70] mil` 实为底层 GND |
| `pours_bottom.csv` | 底层 1 个 Pour = **GND 整板地平面** | 覆盖整板 → 天然参考地 |
| `poured_polygons.json` | 实际铜填充多边形（含挖空），117 个 | 各层 `pourFill` 路径点 |
| `poured_summary.csv` | 每个填充多边形 bbox + 面积(mm²) | 供 openEMS `AddPolygon` |
| `components.csv` | 29 个器件放置（x,y,angle,layer,器件型号,Channel ID） | 供后续定位 L2/U21 |
| `board_openems_template.py` | openEMS 板级几何/网格/NF2FF **模板骨架** | 板框+FR4+底层 GND+30MHz–1GHz 频段已搭好，源待接 |

板框（米，原点=L2）：`x[-6.22, 20.28] mm, y[-32.97, 17.53] mm`；FR4 1.51 mm；铜 35 µm。

> **几何注意事项（沿用 fea3 教训）**：铜 35 µm **不要切成 3D 实体体网格**；openEMS 里用
> **薄层导体（PEC 或高 kappa 薄片）** 表示顶层/底层铜即可，避免网格畸变。
> 器件设计号（L2/U21）在 epru 里通过 `Channel ID`/ATTR 关联，本次先导出坐标+器件型号，设计号映射列为下一步。

---

## 4. 正式板级远场仿真的可行性评估（RAM / 网格 / 耗时）

**结论：可行，但需先解决"源模型"和"绝对值归一化"两件事。**

- **频段 30 MHz–1 GHz**：FDTD 单次宽带高斯脉冲即可覆盖（`SetGaussExcite` 到 1 GHz）。30 MHz 端波长 10 m，结构电小，无额外代价；1 GHz 端 λ=300 mm。
- **网格量级（估算）**：
  - 空气区按 λ_min/20（1 GHz → 15 mm）→ 外边界取板框 ±3λ_min 左右（~150×200 mm 量级），空气网格很粗。
  - board 区域细网格 ~0.2–0.5 mm 以解析铜边沿/走线 → 这是主要开销。
  - **混合网格估计 3M–20M cells**（细化铜边沿会显著抬升上界）。
- **RAM（估算）**：openEMS 压缩 SSE 约 **数十~~150 B/cell** → **3M cells ≈ 0.5–1.5 GB；20M cells ≈ 3–8 GB**，均在 15 GB 内，建议 20M cells 以内。
- **耗时（估算）**：本机实测 ~200 MCells/s（12 线程）。`T ≈ cells × N_TS / speed`。
  - 5M cells × ~2–3 万时间步 → **约 10–20 min/run**；
  - 20M cells → **约 40–80 min/run**。
- **风险点**：1 GHz 与 30 MHz 的动态范围（收敛判据 `EndCriteria` 要足够小）；低端场量小易被数值底噪淹没，建议分段/加大激励或延长时间步。

### 下一步（建议）
1. 把开关节点 dv/dt、输入/输出**电流环**作为真实源接入 `board_openems_template.py`（缺参数需从原理图/实测确认，**不猜测**）。
2. 顶层铜 Pour + 关键回路走线用薄层 PEC 建出；GND 底层已是整板参考。
3. 明确 Prad 绝对值归一化流程（除 `|uf_inc|²`）。
4. 先跑 **1 GHz 单点**最小验证，再扩到 30 MHz–1 GHz 宽带扫描。

---

## 5. 局限
- 时间 ~2 h（略超 1.5 h 盒），主要耗在 VTK 依赖冲突排查与源码编译。
- `AppCSXCAD` GUI 未装（无头环境无需求）。
- 几何仅导出坐标/多边形，**未映射器件设计号**、**未接入辐射源**、**未做网格收敛性验证**。
- Prad 绝对值归一化尚未在脚本中实现（已在 §2 说明依据）。

## 6. 产物清单（绝对路径）
```
/mnt/raid10/sim-work/tps62933/rad/REPORT_rad_prep.md
/mnt/raid10/sim-work/tps62933/rad/smoke_test/Simple_Patch_Antenna.py
/mnt/raid10/sim-work/tps62933/rad/smoke_test/patch_run.log
/mnt/raid10/sim-work/tps62933/rad/smoke_test/pattern_farfield.png   <- 官方远场证据图
/mnt/raid10/sim-work/tps62933/rad/smoke_test/patch_nf2ff.h5
/mnt/raid10/sim-work/tps62933/rad/smoke_test/dipole_nf2ff.py
/mnt/raid10/sim-work/tps62933/rad/smoke_test/dipole_sim/pattern_2p4GHz.png
/mnt/raid10/sim-work/tps62933/rad/geom/export_geom.py
/mnt/raid10/sim-work/tps62933/rad/geom/board_outline.csv
/mnt/raid10/sim-work/tps62933/rad/geom/pours_top.csv
/mnt/raid10/sim-work/tps62933/rad/geom/pours_bottom.csv
/mnt/raid10/sim-work/tps62933/rad/geom/poured_polygons.json
/mnt/raid10/sim-work/tps62933/rad/geom/poured_summary.csv
/mnt/raid10/sim-work/tps62933/rad/geom/components.csv
/mnt/raid10/sim-work/tps62933/rad/geom/board_openems_template.py
/mnt/raid10/sim-work/tps62933/rad/logs/   (apt/csxcad/openems 的 cmake/make/install 日志)
```
