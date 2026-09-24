# REPORT_cispr_split.md — TPS62933 传导 CISPR 32 Class B 独立单面板图

**任务**：把原来"一张图塞 DM+CM × 3 工况 × 3 变体"的大图，拆成**每张只有一个面板、元素最少**的图。
**只画图，未重做仿真**——所有数值来自已有结果/模型的只读复用。
**限值**：CISPR 32 / EN 55032 **Class B**，150 kHz–30 MHz，**QP + AV 双限值**（AV 限值 = QP − 10 dB）。

**产物目录**：`/mnt/raid10/sim-work/tps62933/cispr_split/`
- `cispr_dm_no_pi.png`（DM，无 π 滤波）
- `cispr_dm_with_pi.png`（DM，有 π 滤波）
- `cispr_cm.png`（CM）
- `plot_cispr_split.py`（生成脚本，只写本目录）
- `numbers_cispr_split.txt`

---

## 0. 通用口径与数据来源（精确路径）

| 项 | 值 | 来源（精确路径） |
|---|---|---|
| DM 模型函数 | — | `dm_recheck/dm_recheck.py` 前缀（`lines = []` 之前，只读 exec 复用）；该脚本为 `emi_v2/emi_model_v2.py` 的 bit-for-bit 复现 |
| DM 源谱 | 解析梯形波，tr=tf=10 ns，D=0.5，Ipk=3.0 A | 同上模型 |
| 真实铜皮 L_hot | **2.70 nH**（microstrip over GND plane，l=9.3 mm × w=3.0 mm × h=1.43 mm） | `cispr_final/numbers_cispr_final.txt`（L_field 用 Elmer 2D 交叉验证 = 2.824 nH，±15% 内 → 等级 A） |
| CM 模型函数 | `V_cm = 25 Ω · 2πf · C_p · |V_sw|` | `cm_redo/cm_redo_model.py`（逐字复用） |
| CM 源谱 v(sw) | FFT 末 32 周期，列 4（v(sw)） | `results_v2/E_full.txt`（本图用满载；另 `E_noload.txt`/`E_half.txt`） |
| C_p | **2.416 pF** | `cm_redo/numbers_cm_redo.txt` / `cm_audit`（真实铜皮 2 层，core 1.43 mm FR4，top coplanar GND） |
| π 段定义 | L1(1 µH FXL0420-1R0-M) + PPHV_OUT_FILTER 组 (C73‖C74=20 µF + C72=100 nF)；30 µF 银行在 L1 前，**不属 π 段** | `cispr_final_nopi/numbers_nopi.txt` |

**只读声明**：`emi_v2/`、`results_v2/`、`dm_recheck/`、`cm_redo/`、`cm_fix/`、`cispr_final*/` 一律**未改动**；本脚本仅 `read/exec` 其模型，输出只写 `cispr_split/`。

**显示工况**：三张图统一取 **满载 3.0 A**（最恶劣的真实工况）。原大图含 noload/half/full 三档，本拆分图各只放**该工况 1 条噪声曲线 + QP/AV 限值线**（图例 3 项），符合"元素最少"要求。

---

## 1. `cispr_dm_no_pi.png` — DM，无 π 滤波

- **曲线**：DM 噪声（整段 π 移除 = `nopi_b`，红实线+圆点）
- **限值**：CISPR 32 Class B QP（黑实线）+ AV（黑虚线）
- **最差裕量**：**QP −25.1 dB @ 4.03 MHz**；AV −35.1 dB @ 4.03 MHz（图中红圈标注）
- **数据来源**：`dm_recheck/dm_recheck.py` 前缀 + 本目录 `plot_cispr_split.py`（`Zt_mode('nopi_b')`）
- **可信度：B**（相对趋势可信；绝对电平依赖假定 tr=tf=10 ns 与源谱；L_hot=2.70 nH 为 A）

## 2. `cispr_dm_with_pi.png` — DM，有 π 滤波

- **曲线**：DM 噪声（完整 π，`pi`，蓝实线+圆点）
- **限值**：同上 QP + AV
- **最差裕量**：**QP +25.1 dB @ 4.03 MHz**（过限值）；AV +15.1 dB @ 4.03 MHz（图中红圈标注）
- **数据来源**：同上，`Zt_mode('pi')`；**与 `cispr_final/numbers_cispr_final.txt` 的 full 行 (+25.1 / +15.1) 逐位一致**
- **可信度：B**

### π 滤波器对 DM 的压制量（满载，4.03 MHz 最差点）
| 对比 | 发射差 | 裕量提升 |
|---|---|---|
| 有 π − 无 π-a（只短路 L1） | **−47.5 dB** | **+47.5 dB** |
| 有 π − 无 π-b（整段 π 移除） | **−50.3 dB** | **+50.3 dB** |

即：仅去掉 L1 就让 DM 从 +25.1 dB 掉到 −22.3 dB（差 47.5 dB）；整段 π 移除再掉 2.8 dB。

## 3. `cispr_cm.png` — CM（满载）

- **曲线**：CM 噪声（绿实线+圆点，`full`，C_p=2.416 pF）
- **限值**：CISPR 32 Class B QP（黑实线）+ AV（黑虚线）
- **最差裕量**：**QP −16.0 dB @ 0.81 MHz**（超限值）；AV −26.0 dB @ 0.81 MHz（图中红圈标注）
- **数据来源**：`cm_redo/cm_redo_model.py` 逐字复用；v(sw) 取自 `results_v2/E_full.txt` 列 4；C_p=2.416 pF
- **可信度：B**（集总 V_cm 近似）；其中 **C_p=2.416 pF 为 A、V_sw 源为 A**

### 为什么没有"CM 有 π / 无 π"两张图（该组合不适用，未硬造）
π 段元件（L1 + 电容）**都回流到板内 GND**；CM 环路是
`SW → C_p(板↔大地) → LISN(25 Ω) → 输入电缆 → 源`，**π 段没有任何元件落在该环路里**。
因此在当前口径下 **CM 曲线对 π / 无 π-a / 无 π-b 完全相同**，拆成"有/无 π"两张图是**无意义的重复**，故只出**一张 CM 图**。（此为 09/`cispr_final_nopi` 已确认的结论。）

---

## 4. 可信度汇总（沿用 09 分级）

| 结论 | 等级 | 依据 |
|---|---|---|
| DM 满载 **有 π 通过**（QP +25.1 / AV +15.1 dB @4.03 MHz） | **B** | 相对可信；绝对值依赖假定 tr=tf=10 ns + 源谱；**AV 检测未做 QP/AV 加权**，仅幅度对限值比较 |
| DM 满载 **无 π 严重超标**（QP −22.3…−25.1 dB） | **B** | 同上 |
| **π 段 DM 压制 ≈ 47.5–50.3 dB** | **B** | 同一模型 Zt 比值；随 L1/布局变化 |
| **CM 满载超标**（QP −16.0 @0.81 MHz，AV −26.0） | **B** | 集总模型；趋势可靠 |
| C_p = 2.416 pF、V_sw 源谱、L_hot = 2.70 nH | **A** | 真实铜皮/几何抽取，已交叉验证 |
| 半载 CM（未达稳态） | **C** | 上游标记：SPICE half 未达稳态，仅参考（本图**未**采用） |

**能信 / 不能信**：
- **能信**：DM 有 π 裕量充裕、无 π 严重超标；π 段压制约 47–50 dB；CM 三工况均超标、最差满载。
- **不能信（绝对）**：DM/CM 的**绝对 dBµV 数字**——受 tr/tf 与集总近似影响；且**AV 结论未做检测器加权**，真机 AV 可能比此处更接近限值。

---

## 5. 局限

1. **只画图，未重算**：DM 复用 `dm_recheck` 模型、CM 复用 `cm_redo` 模型，二者分级沿用 09，本次未新增物理验证。
2. **只显示满载**：为满足"图例 ≤3 / 元素最少"，单图只放满载曲线；noload/half 数字见 `numbers_cispr_split.txt` 与既有报告。
3. **CM 起始频率**：CM 谐波从 fsw=805 kHz 起（无 150–805 kHz 数据点），图中 QP 限值线在该段仅作参考。
4. **AV 检测器未加权**：模型用 CW 谐波幅度直接对 QP/AV 限值比较，未做 CISPR 检波器加权，AV 结论偏乐观/需实测确认。
5. 半载 CM 未采用（C 级，未达稳态）。
