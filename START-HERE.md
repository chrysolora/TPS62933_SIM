# START HERE — 这个仓库怎么看

> 一句话：**这是个 TPS62933（24V→12V/3A）降压电路的完整仿真结果集**。
> **先读 [`docs/board-facts.md`](docs/board-facts.md)（本板唯一真源）和 [`docs/09-credibility-rubric.md`](docs/09-credibility-rubric.md)（可信度分级 A/B/C），否则容易把"过程数据"当成"结论"。**

---

## 1. ✅ 核心结果（按重要性，看这 6 个）

| 看什么 | 目录 | 一句话结论 | 等级 |
|---|---|---|---|
| **传导 EMI（最终）** | [`results/20-conducted-cispr-final/`](results/20-conducted-cispr-final/) | **差模 DM 稳过（QP +25 dB）；共模 CM 全超（QP −16 dB）** | A/B |
| **辐射 EMI** | [`results/13-rad-converged/`](results/13-rad-converged/) | **三版铺铜在辐射上"没有可分辨差异"** | A−/B |
| **板级热（最终）** | [`results/24-thermal-final/`](results/24-thermal-final/) | **板均温 73℃；U21 结温 120~139℃，<150℃ 通过** | A/B |
| **真实铜皮重算 DM** | [`results/18-dm-real-copper/`](results/18-dm-real-copper/) | 按**本板真实铜皮**量出 L_hot≈2.7 nH → **DM 确实稳过** | A−/B+ |
| **铺铜涡流损耗** | [`results/12-eddy-loss/`](results/12-eddy-loss/) | **涡流损耗 ≪ 铜损（<0.1%）→ "挖铜省损耗"可忽略** | B |
| **CM 修正** | [`results/17-cm-fix/`](results/17-cm-fix/) | CM 超限方向确凿；幅度取决于建模（1–16 dB） | A−/B |

**📌 只想看两张图**：
- 传导 EMI → `results/20-conducted-cispr-final/figures/fig_conducted_cispr_final.png`
- 热 → `results/24-thermal-final/figures/fig_thermal_final_map.png`

**补充目录**：`results/25-cispr-split/`（传导拆成 3 张单图）｜`results/26-thermal-3x3/`（三铺铜×三负载温度场）｜`results/27-geom-verify/`（铺铜几何核验）｜`results/23-thermal-v2/`（热步骤1留档）

---

## 2. 📄 三个必读文档

1. [`docs/board-facts.md`](docs/board-facts.md) — **本板唯一真源**（拓扑、器件、层叠、几何、修正记录）。
2. [`docs/09-credibility-rubric.md`](docs/09-credibility-rubric.md) — **可信度分级 A/B/C**。
3. [`docs/07-input-conditions-and-assumptions.md`](docs/07-input-conditions-and-assumptions.md) + [`docs/08-upstream-supply-LM50-20B24.md`](docs/08-upstream-supply-LM50-20B24.md) — 前提与上游电源参数。

---

## 3. 结论速览（五个维度，三版铺铜对照）

| 维度 | 结论 | 来源 |
|---|---|---|
| **传导 DM** | **稳过** Class B（QP +25 / AV +15 dB） | 20 / 18 |
| **传导 CM** | **全超** ~16 dB（π 滤波无效，真瓶颈） | 20 / 17 |
| **辐射** | 三版**无可分辨差异**；绝对量偏悲观（无外壳/屏蔽） | 13 |
| **涡流屏蔽** | 全铺铜电感下 \|J\| 降 **93–98%**（相对量） | 12 |
| **涡流损耗** | **可忽略**（<0.1% 铜损） | 12 |
| **散热** | 挖铜 → 板均温几乎不变，热点略集中 | 24 / 26 |

> **一句话总账**：**全铺铜在"屏蔽 + 散热"上占优；挖铜只在"CM 寄生电容"上占优；其它维度基本无差。**

---

## 4. 尚未完成

- **辐射的绝对电平**（现在是"相对可信、绝对偏悲观"）
- **实测验证**：满载输入功率 / tr,tf / 实际传导发射
- **磁涡流三版重做**（进行中）

---

*更新于 2026-09-24。已被取代的旧版目录已移除，保留在 git 历史中。*
