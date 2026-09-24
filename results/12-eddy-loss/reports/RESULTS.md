# tps62933-loss-redo — P_eddy(W) 三版×三工况 + 温升粗估

## 0. 结论摘要（诚实版）
- 无法给出**可信的绝对** P_eddy：现有 fea4（3D, CoilSolver/AV 双势方法）**激励标定不可信**
  （colil 场自洽性不一致 + 三个版本 bottom 层结果出现物理上不可能的 160x 差异）。
- 可给出**相对趋势**（定性）：电感 L2 正下方（top 层）铜皮挖空后，该层涡流损耗下降 ~85–92%。
- 数量级：即便按 fea4 数值取，实际工况（电感纹波 0.19–0.30 A 幅值）下 P_eddy ≈ 0.1–0.6 µW 量级；
  即使按“实心短路盘”解析上界（最保守）也只 ~mW 量级 → 相对 DCR 铜损(~0.1–0.45 W)占比 <0.1% →
  **对温升无可测影响**。挖铜省损耗/降温的收益在此问题上可以忽略。

## 1. 激励电流来源（SPICE）
- 来源文件：results_v2/E_{noload,half,full}.txt（ngspice wrdata，6 个信号 × (time,val)）。
  列映射：v(vout)=c1, v(sw)=c3, i(vin)=c5, **i(l2)=c7**（0-based 列 6=time,7=值），v(fb)=c9, v(vin)=c11。
  已用 t=0 值核对：noload→0.0200A, half→1.4999A, full→2.9999A（=各工况负载电流）。
- 取末段(后25%)稳态做 FFT：
  | 工况 | Imean(A) | Iac_rms(A) | 纹波pkpk(A) | fsw基波 | 基波幅值(A) |
  |noload| ~0.04 | 0.065 | 0.22 | 796 kHz | 0.074 |
  |half  | ~2.24 | 0.52  | 2.06 | 54 kHz (!) | 0.57 (未settle，见局限) |
  |full  | ~4.0  | 0.15  | 0.60 | 779 kHz | 0.187 |
- 涡流损耗由**纹波(AC)**驱动；DC 不产生涡流。取 fsw 基波幅值作激励是合适的（三角纹波总损耗≈1.64×基波）。
- 注意：half/noload 的 Imean 与负载电流不符（half 2.24 vs 1.5A），E_half 存在 54kHz 慢振荡 → 该工况激励不可信。

## 2. fea4 数据再处理（唯一含完整几何的 3D 模型）
- fea4/{fullcu,topcut,dualcut}/mesh/case_t0001.vtu 中**已含 Elmer “joule heating” 场**（此前只出 |J| 图，从未积分）。
- 本任务用 tetra 体积加权在铜层 z 带内积分（top: z∈[-35µm,0]；bot: z∈[-1.035,-1.0]mm；空气 σ=0 不贡献）。
- 独立复核：用 ∫ ½ρ|J|² (ρ=1/σ=1.724e-8) 自行积分，与 Elmer joule heating 一致到 ~20–40%（验证积分/量纲正确；
  关键：正确公式是 ∫ρ|J|²，**不是** ∫σ|J|²——后者会差 σ²≈3.4e15 倍）。

### P_eddy @ fea4 激励 (Desired Coil Current = 6.0 A, f=784 kHz, µr(core)=40)
| 版本 | top层 (W) | bot层 (W) | vs FullCopper(top) |
|FullCopper| 3.94e-4 | 2.22e-4 | — |
|TopCutout | 3.25e-5 | 1.44e-6 | −91.8% |
|DualCutout| 5.51e-5 | 2.23e-5 | −86.0% |
(独立 ∫ρ|J|²/2: fullcu top 3.14e-4 / bot 1.50e-4; topcut 2.48e-5/1.07e-6; dualcut 4.44e-5/1.80e-5 —— 同趋势)

### 不可信证据（关键）
1) **bottom 层跨版本矛盾**：TopCutout 的 bottom 铜皮完整（体积 4.23e-8 同 fullcu），但 bottom损耗
   比 fullcu 低 154x；而 DualCutout（把 bottom 也挖了）反而比 TopCutout 高 15x。物理上不可能 → 结果不自洽。
2) **B 与 J 不自洽**：铜层内 |B|max≈4.6 mT，但 “current density re” max 仅 1.76e7 A/m²；
   由 J≈σωBr/2 与 B 对应的预期 J 大 ~1–2 个量级 → 符合已知的 CoilSolver/AV 双势激励标定缺陷
   （fea5 失败的同一根因，只是 fea4 表现较隐蔽）。
3) **网格收敛失败**：把 fullcu 网格加密(近铜 lc 0.4→0.2mm)后，BiCGStabl 迭代 1000 步残差仍 0.06（"Too many iterations"）
   → 加细即不收敛，无法给出 <10% 的收敛证据。

## 3. 自检：2D 解析基准（验证“损耗后处理+求解器”）
- 采用 2D planar (Az) 谐波涡流，回避 CoilSolver：薄铜板(6mm×35µm)置于均匀切向 AC 场 B0=1mT。
- Elmer 计算 P=0.01256 W/m；解析薄板解 P=σω²B0²wt³/24=0.01509 W/m；**比值 0.83**
  （t/δ=0.475，薄板极限偏高、屏蔽使之偏低，0.83 合理）→ **2D 求解器 + ∫½ρ|J|²(Az) 后处理链已验证**。
- 结论：2D 路线在方法上可行（若激励接对，可给出可信绝对值）。

## 4. 尝试的真实 2D 模型（未成功，如实记录）
- 建了全几何：线圈双腿(±2.0–2.6mm, z 0.9–1.5mm)、磁芯(±3.25mm,µr=40)、
  top/bot 铜皮（可挖 x<xc=3.3mm 中心区）、空气；外边界 Potential=0。
- **失败点**：以 BodyForce “Current Density” 施加已知线圈电流密度时，2D 谐波求解器报
  “Material 1 electric conductivity / relative permeability” 未使用 → 铜 σ 未被采纳 → 无涡流，
  解退化为病态静磁解（Az~1e9，P 无意义）。**激励未能接通**。
- 由于时间盒用尽，未继续修复（下一步：2D 单线圈 CoilSolver 或改用 “Current Density”/“jfix”正确关键字验证）。

## 5. 温升粗估
- 实际工况 P_eddy(纹波激励) ≈ (I_ripple/6A)² × P_eddy@6A × 1.64：
  full: top(fullcu)=3.94e-4×(0.187/6)²×1.64≈6.2e-7 W；noload≈9.7e-8 W；half 不可信。
- 对比 DCR 铜损 ~0.1–0.45 W → P_eddy 占比 <0.1% → 温升贡献 ≲ 0.01°C 量级，**可忽略**。
- 甚至用最保守解析上界(实心短路盘 ~0.5W@6A)外推，也只 ~1mW 量级，仍 <1% 铜损。

## 6. 局限 / 未决
- fea4 激励标定缺陷未修复 → 绝对 P_eddy 不可用；跨版本对比仅定性。
- 2D 真实模型激励未接通 → 无替代可信绝对值。
- E_half SPICE 未稳态；fsw 取 784–805kHz（fea4 用 784kHz）。
- 磁芯 µr=40 为估值（真实铁氧体更高，会增大损耗上界）；未做完整热仿真。
