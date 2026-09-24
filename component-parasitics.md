# component-parasitics.md — TPS62933F 输入板 关键器件寄生/损耗参数（带来源）

> 依 `skills/simulate-from-board-copper` + `board-facts.md` §7 GAPS 建立。
> 来源级别：**D**=器件数据手册/厂商官网 ｜ **E**=权威应用手册/厂商曲线典型值（绝对值不可引用，须标注）｜ **X**=理论坏布局上界（**禁止**当本板值）。
> 建立日期：2026-09-24（子代理联网抓取）。**未从权威来源确认者一律进 §GAPS，不猜测、不填典型值。**
> 抓取工具：web_fetch（DuckDuckGo HTML 检索）+ curl + pypdf（PDF 文本抽取）。web_search API 不可用。

---

## 参数表

| 器件位号 | 型号 | 参数 | 数值 + 单位 + 条件 | 级别 | 来源(URL/页) |
|---|---|---|---|---|---|
| C17,C62,C63,C73,C74 | CL31A106KBHNNNE | 容值/容差/耐压/介质/封装 | **10 µF · ±10% · 50 V · X5R · 1206**（3.20×1.60×1.60 mm，T=1.60mm） | **D** | Samsung SPEC (Reference Sheet) `mlccsupplier.com/wp-content/uploads/2025/07/SpecSheet_CL31A106KBHNNN.pdf`；交叉 `product.samsungsem.com/mlcc/CL31A106KBHNNNE.do` |
| ″ | ″ | 损耗角正切 DF (tanδ) | **≤ 0.125 max**（1 kHz；可靠性判据） | **D** | 同上（Samsung SPEC，C 节 Reliability test） |
| ″ | ″ | **ESR @100 kHz** | **未给出（手册无数值）** | **GAPS** | Sony/ Samsung 手册仅定性解释 ESR/ESL/SRF 概念（§2-8-1/2/3），**无曲线/数值** → 需 Samsung **SimSurfing / S2P** 模型 |
| ″ | ″ | **ESL** | **未给出** | **GAPS** | 同上 |
| ″ | ″ | **SRF** | **未给出** | **GAPS** | 同上 |
| C16,C64 | JVJ50v10M5x5 | 器件类型 | **SMD 铝电解电容（Aluminum Electrolytic，非聚合物）** | **D** | LCSC/JLCPCB 商品参数（`item.szlcsc.com/48680110.html`；`jlcpcb.com/partdetail/jieerrui-JVJ50v10M5x5/C46550416`） |
| ″ | ″ | 容值/容差/耐压/封装 | **10 µF · ±20% · 50 V · SMD D5×L5.4 mm**；−55…+105 °C；2000 h@105 °C | **D** | 同上 |
| ″ | ″ | 额定纹波电流 | **18 mA @120 Hz** | **D** | 同上 |
| ″ | ″ | **ESR** | **未给出**（铝电解手册未列 ESR 数值） | **GAPS** | 需厂家 tanδ 曲线或实测 |
| ″ | ″ | **ESL** | **未给出** | **GAPS** | — |
| C1 | HGC1210R5476M250NSVK | **容值/耐压/介质/封装（确认）** | **47 µF · ±20% · 25 V · X5R · 1210** (MLCC) | **D** | LCSC `lcsc.com/product-detail/C7432791.html`；GlobalSpec `datasheets.globalspec.com/.../hgc1210r5476m250nsvk`；JLCPCB `partdetail/Chinocera-HGC1210R5476M250NSVK/C7432791` |
| ″ | ″ | ESR / ESL | 未给出 | **GAPS** | Chinocera 手册未列 |
| C49,C51,C72 | CC0603KRX7R9BB104 | 容值/耐压/介质/封装 | **100 nF · 50 V · X7R · 0603**（Yageo CC 系列） | **E** | 依 Yageo 型号命名规则解码（未逐项核对 Yageo 手册 PDF；见 GAPS） |
| ″ | ″ | ESR / ESL | 未给出 | **GAPS** | 0603 X7R 手册通常不列 ESR/ESL |
| C18 | 0603CG2R7C500NT | 容值/介质/封装/耐压 | **2.7 pF · C0G(NP0) · 0603 · 50 V** | **E** | 依型号命名规则解码（未逐项核对手册；见 GAPS） |
| ″ | ″ | ESR / ESL | 未给出 | **GAPS** | — |
| L1 | FXL0420-1R0-M | 电感值/容差 | **1 µH · ±20%** | **D** | LCSC/JLCPCB 数据手册：`lcsc.com/datasheet/C167203.pdf`（元数据）/ `jlcpcb.com/partdetail/178586-FXL0420_1R0M/C167203` |
| ″ | ″ | **DCR** | **27 mΩ**（typical） | **D** | 同上（JLCPCB/LCSC 规格：`1uH 27mΩ 6A 7A Molded Inductor ±20% SMD,4.4x4.2mm`） |
| ″ | ″ | 额定电流 Irms | **6 A** | **D** | 同上 |
| ″ | ″ | 饱和电流 Isat | **7 A** | **D** | 同上 |
| ″ | ″ | **SRF** | **未公布** | **GAPS** | 数据手册未载 |
| L2 | ZEMS0650-150M | 电感值/容差 | **15 µH · ±20%**（一体成型/合金电感・7×6.6×4.8 mm） | **D** | LCSC `item.szlcsc.com/43235556.html`（C41415624）；规格书 PDF 经 atta.szlcsc.com |
| ″ | ″ | **DCR** | **85 mΩ** | **D** | 同上 |
| ″ | ″ | 额定电流 | **4 A** | **D** | 同上 |
| ″ | ″ | 饱和电流 Isat | **5 A** | **D** | 同上 |
| ″ | ″ | **磁芯材料 / µr** | **手册未给出** | **GAPS** | 仅知“一体成型/合金(molding)电感”；µr 无依据（board-facts 现用 40=**E**） |
| U21 | TPS62933FDRLR | RdsON 高侧 (RDSON_HS) | **76 mΩ** @ TJ=25°C, VBST−SW=5 V | **D** | TI **SLUSEA4D** 数据手册 §8.5 Electrical Characteristics（pdf，`ti.com/lit/ds/symlink/tps62933.pdf`） |
| ″ | ″ | RdsON 低侧 (RDSON_LS) | **32 mΩ** @ TJ=25°C | **D** | 同上 |
| ″ | ″ | 静态电流 Iq（**TPS62933F**） | **125 µA**（Non-switching，EN=5 V, VFB=1 V）。注：非 F 版（PFM）为 12 µA；本板为 **F 版=125 µA** | **D** | 同上 |
| ″ | ″ | 结温上限 | **TJ(max)=150 °C**（工作）；热关断 TSHDN=165 °C | **D** | 同上 §8.5 / §7.1 |
| ″ | ″ | 开关频率 fSW | **500 kHz（RT 悬空）· 310 kHz@RT=71.5 kΩ · 2100 kHz@RT=9.09 kΩ**；本板 RT=26.7 kΩ → 约 **784 kHz** | **D**(表)/**E**(26.7k 插值) | TI SLUSEA4D OSCILLATOR FREQUENCY 表 |
| ″ | ″ | **SW 上升/下降时间 tr/tf** | **手册未给** → 暂用 **5 ns**（TI 典型，**E**） | **E** | board-facts §4（须实测替换） |
| D1 | SS36 (LCSC C7420367) | 器件身份 | 肖特基整流 3 A / 60 V（SMC/DO-214AB），串联输入二极管 | **D**(身份)/**GAPS**(数值) | 网表/netlist（board-facts §3）；有效数据手册 PDF 未取到（见 GAPS） |
| ″ | ″ | **Vf @额定** | **未取到有效手册值** | **GAPS** | 见 GAPS |
| ″ | ″ | **结电容 Cj** | **未取到有效手册值** | **GAPS** | 见 GAPS |

---

## GAPS（抓不到 / 需人确认）

1. **Samsung CL31A106KBHNNNE 的 ESR@100kHz / ESL / SRF**：Samsung 规格书只有定性说明（§2-8-1..3），**无任何数值/曲线**。→ 需从 Samsung **SimSurfing** 取 S2P/阻抗-频率模型（本次未能联网导出其 S2P，`weblib/SimSurfing` 为 JS 交互，PDF 不含曲线数据）。**勿用“典型 MLCC”补脑。**
2. **JVJ50v10M5x5（SMD 铝电解）的 ESR / ESL**：铝电解手册未列 ESR 数值（仅给纹波 18mA@120Hz）。→ 需厂家 tanδ-频率曲线或**实测**；ESL 通常需实测。
3. **HGC1210R5476M250NSVK 的 ESR / ESL**：Chinocera 手册未列。→ 同 MLCC，需厂家阻抗模型/实测。
4. **CC0603KRX7R9BB104 / 0603CG2R7C500NT 的容值与封装身份**：本次未取到 Yageo 及其它原厂**手册 PDF 原文**，仅按型号命名规则解码为 **100nF/50V X7R 0603** 与 **2.7pF/50V C0G 0603**（标 **E**）。→ 建议**核对原厂手册**确认（尤其耐压码 9=50V、104=100nF、2R7=2.7pF）。二者的 ESR/ESL 亦未公布（GAPS）。
5. **FXL0420-1R0-M 的 SRF**：数据手册未公布（仅得 DCR 27mΩ / 6A / 7A）。→ 需原厂手册或实测。**另**：DCR/Irms/Isat 取自 JLCPCB/LCSC 规格条目（源自厂家手册元数据），级别记 **D**，但**未逐行核对原厂 PDF 表格**，建议复核。
6. **ZEMS0650-150M 的磁芯材料与 µr**：一体成型/合金电感，手册未给 µr。board-facts 现用 **µr=40（E，无依据）**。→ 需原厂手册或实测（影响涡流 |J| 线性比例，形态不变）。
7. **SS36 的 Vf 与结电容 Cj**：多次尝试（onsemi/Diodes/MCC/MDD/LCSC PDF）均被 Cloudflare/JS 拦截，未取到**任何有效 SS36 数据手册原文**（曾误取 SK12–SK16=1A 文档，非本器件，已弃用）。→ 需人工提供 SS36 原厂手册（onsemi SS36=M BRS360 系列 / MDD / 台半 TSC），或实测 Vf@额定与 Cj。
8. **TPS62933F 的 SW 上升/下降时间 tr/tf**：TI SLUSEA4D **未给**该项规格。→ 仍为 board-facts 的 **E=5 ns**，**须实测替换**（VHF 段包络关键量）。

---

## 附：本次抓取的方法/局限

- `web_search`（Brave）**无 API key，不可用**；改用 `web_fetch` 抓 `html.duckduckgo.com/html/?q=` 做检索，再抓目标页/PDF。
- 部分厂商站点（onsemi、Mouser、TME、LCSC 商品页）为 Cloudflare/JS 渲染，`curl`/`web_fetch` 均取不到正文 → 相关器件（SS36、FXL0420 原厂 PDF）落入 GAPS。
- 成功用 `curl + pypdf` 抽取 PDF 正文的：TI **SLUSEA4D**（TPS62933 系列，51 页）、Samsung **CL31A106KBHNNNE SPEC**（38 页镜像）。
- **未改动**任何 `emi*/fea*/results*/rad*/cm_*/loss_redo/` 目录；未 git push、未动 Gitea、未对外发消息。
