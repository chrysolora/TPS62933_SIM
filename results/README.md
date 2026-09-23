# 结果目录说明

| 目录 | 内容 | 状态 / 可信度 |
|---|---|---|
| `01-stage1-transient/` | ngspice 时域：纹波、开关波形、效率、网表 | ⚠️ 工程近似（见 `docs/04` A 节） |
| `02-conducted-emi/` | 频域 DM 传导 EMI（CISPR 32）+ 软启动 | ⚠️ 理想化偏乐观（见 `docs/04` B 节） |
| `03-field-eddy-copper/` | 电感下铺铜涡流场 FEA 图 + 脚本 | ⚠️ 链路验证级（见 `docs/04` C 节） |

每个目录内保留其**原始报告**（`REPORT_*.md` / `README.md`），本仓库根的 `docs/` 为**汇总方法/过程/局限**。
