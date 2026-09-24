# 20 — 传导 EMI 正式版（DM + CM 统一口径，CISPR 32 Class B）

**结论：差模(DM) 稳过；共模(CM) 全超。**
- DM 满载：QP **+25.1 dB** / AV **+15.1 dB**（保守 L_hot=4.44nH 时仍 +21.3/+11.3）→ 过
- CM 满载：QP **−16.0 dB** / AV **−26.0 dB** → 超（C_p 2.42–2.57 pF 范围内均超）
- 口径：DM 用**本板真实铜皮几何** L_hot=2.70 nH（微带）｜CM 用真实几何 C_p=2.416 pF
- **L_hot 已由磁场场提取交叉验证**：Elmer MagnetoDynamics2D → 2.824 nH（与微带公式 +4.6%，互差<±15%）→ **L_hot 等级 B→A**
- 等级：L_hot **A**｜C_p **A**｜V_sw 源 **A**｜DM/CM 裕量 **B**（依赖假定 SW 边沿 tr=tf=10ns / 集总 V_cm 模型）
- 图：`figures/fig_conducted_cispr_final.png`（QP：DM 实线 + CM 虚线 + 限值；AV 副图）
- 局限：2D 截面×段长的 L（忽略弯折）；铺铜净空未抠出（已知缺陷）；半载非稳态仅参考。

脚本/报告：`scripts/`、`reports/`。真源：`docs/board-facts.md`。
