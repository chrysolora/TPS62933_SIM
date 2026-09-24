# 25 — 传导 CISPR 拆分单图（DM 无π / DM 有π / CM）

来源：只读复用 `emi_v2/emi_model_v2.py`（经 `dm_recheck/` bit-for-bit 复现）与 `cm_redo/cm_redo_model.py`；**未重做仿真**。满载 3.0A，CISPR 32 Class B QP/AV。

- **DM 无 π**：QP **−25.1 dB @4.03 MHz**（AV −35.1）
- **DM 有 π**：QP **+25.1 dB @4.03 MHz**（AV +15.1）→ π 段相对"整段移除"压制 **+50.3 dB**
- **CM（满载）**：QP **−16.0 dB @0.81 MHz**（AV −26.0）；C_p=2.416 pF
- **"CM 有/无 π"为不适用组合**（π 元件回流板内 GND，不在 CM 环内）→ 未出图
- 可信度：DM/CM 结论 **B**（依赖 tr=tf=10ns、集总近似；AV 未做检波器加权）；L_hot/C_p/V_sw 源 **A**
