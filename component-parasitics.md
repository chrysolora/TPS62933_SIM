# component-parasitics.md — TPS62933F 输入板 关键器件寄生/损耗参数（带来源）

> 依 `skills/simulate-from-board-copper` + `board-facts.md` §7 GAPS 建立。
> 来源级别：**D**=器件数据手册/厂商官网 ｜ **E**=权威应用手册/厂商曲线典型值（绝对值不可引用，须标注）｜ **X**=理论坏布局上界（**禁止**当本板值）。
> 建立日期：2026-09-24（子代理联网抓取）。**未从权威来源确认者一律进 §GAPS，不猜测、不填典型值。**
> 抓取工具：web_fetch（DuckDuckGo HTML 检索）+ curl + pypdf（PDF 文本抽取）。web_search API 不可用。
> 二轮补充工具：Brave/DDG HTML 检索（含 captcha 回退）＋ **r.jina.ai 阅读器**（穿透 Cloudflare/JS 页，可读 LCSC PDF）＋ Samsung EC-Datasheet 内嵌 JSON 解析。

---

## 参数表

| 器件位号 | 型号 | 参数 | 数值 + 单位 + 条件 | 级别 | 来源(URL/页) |
|---|---|---|---|---|---|
| C17,C62,C63,C73,C74 | CL31A106KBHNNNE | 容值/容差/耐压/介质/封装 | **10 µF · ±10% · 50 V · X5R · 1206**（3.20×1.60×1.60 mm，T=1.60mm） | **D** | Samsung SPEC (Reference Sheet) `mlccsupplier.com/wp-content/uploads/2025/07/SpecSheet_CL31A106KBHNNN.pdf`；交叉 `product.samsungsem.com/mlcc/CL31A106KBHNNNE.do` |
| ″ | ″ | 损耗角正切 DF (tanδ) | **≤ 0.125 max**（1 kHz；可靠性判据） | **D** | 同上（Samsung SPEC，C 节 Reliability test） |
| ″ | ″ | **ESR @100 kHz** | **6.5 mΩ**（25 °C, 0 Vdc；三星官方 EC 模型，**详见 §ESR/ESL/SRF**） | **D(推导)** | Samsung EC 模型 `weblib.samsungsem.com/mlcc/mlcc-ec.do?partNumber=CL31A106KBHNNN` |
| ″ | ″ | **ESL** | **≈0.34 nH**（官方模型反推；详见 §ESR/ESL/SRF 注①） | **D(推导)** | 同上 |
| ″ | ″ | **SRF** | **≈1.74 MHz**（\|Z\| 最小点） | **D(推导)** | 同上 |
| C16,C64 | JVJ50v10M5x5 | 器件类型 | **SMD 铝电解电容（Aluminum Electrolytic，非聚合物）** | **D** | LCSC/JLCPCB 商品参数（`item.szlcsc.com/48680110.html`；`jlcpcb.com/partdetail/jieerrui-JVJ50v10M5x5/C46550416`） |
| ″ | ″ | 容值/容差/耐压/封装 | **10 µF · ±20% · 50 V · SMD D5×L5.4 mm**；−55…+105 °C；2000 h@105 °C | **D** | 同上 |
| ″ | ″ | 额定纹波电流 | **18 mA @120 Hz** | **D** | 同上 |
| ″ | ″ | **ESR @120 Hz** | **≈18.6 Ω**（=tanδ/(2πfC)，tanδ=0.14@120 Hz；**详见 §ESR/ESL/SRF**） | **D(推导)** | jieerrui JVJ series datasheet（`item.szlcsc.com/48680110.html`） |
| ″ | ″ | **ESL** | **未给出** | **GAPS** | — |
| C1 | HGC1210R5476M250NSVK | **容值/耐压/介质/封装（确认）** | **47 µF · ±20% · 25 V · X5R · 1210** (MLCC) | **D** | LCSC `lcsc.com/product-detail/C7432791.html`；GlobalSpec `datasheets.globalspec.com/.../hgc1210r5476m250nsvk`；JLCPCB `partdetail/Chinocera-HGC1210R5476M250NSVK/C7432791` |
| ″ | ″ | ESR / ESL | 未给出 | **GAPS** | Chinocera 手册未列 |
| C49,C51,C72 | CC0603KRX7R9BB104 | 容值/耐压/介质/封装 | **100 nF · 50 V · X7R · 0603**（Yageo GP X7R 系列） | **D** | Yageo MLCC 规格书 `lcsc.com/datasheet/C14663.pdf`（LCSC C14663，GP X7R 6.3–250 V） |
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
| D1 | SS36 (LCSC C7420367) | 器件身份 | 肖特基整流 3 A / 60 V，串联输入二极管；**封装=DO-214AC (SMA)**（⚠ 与 board-facts “SMC/DO-214AB” 不符） | **D** | 珠海宏嘉诚 HJC **SS32–SS320** 手册 `lcsc.com/product-detail/C7420367.html` |
| ″ | ″ | **Vf @3 A** | **0.70 V (max)** @ IF=3.0 A, Ta=25 °C | **D** | 同上 |
| ″ | ″ | **结电容 Cj** | **300 pF (typ)** @ VR=4.0 V, 1 MHz | **D** | 同上 |

---

## ESR / ESL / SRF

> 本轮（2026-09-24 二轮）专攻上一轮留下的 ESR/ESL/SRF/损耗缺口。
> **突破点**：Samsung MLCC 的 EC-Datasheet 是 JS 交互页，但**其频率模型数据（\|Z\|、ESR vs f）内嵌在页面 JSON 中**，可直接解析 → 拿到 CL31 的官方模型数值。
> 级别：**D**=原厂手册/厂商官网；**D(推导)**=用原厂模型/手册数据按明确公式换算（方法随值注明）；**E**=权威典型（绝对值不可引用）；**X**=禁止。
> 本地模型/档案：`node249:/mnt/raid10/sim-work/tps62933/cap_models/`。

| 器件位号 | 型号 | 参数 | 数值 + 条件 / 频率 | 级别 | 来源 URL |
|---|---|---|---|---|---|
| C17,C62,C63,C73,C74 | CL31A106KBHNNNE | **ESR @100 kHz** | **6.5 mΩ**（25 °C, 0 Vdc；对官方模型 ESR-f 曲线 log-log 插值） | **D(推导)** | Samsung EC Data Sheet 模型 `https://weblib.samsungsem.com/mlcc/mlcc-ec.do?partNumber=CL31A106KBHNNN`（存档 `cap_models/CL31A106KBHNNN_samsung_ECmodel.json`；页面 `.../CL31A106KBHNNN_samsung_ECdatasheet_view.html`） |
| ″ | ″ | ESR @1 kHz（参考） | 357 mΩ（同模型；与手册 DF≤0.125@1 kHz 量级一致） | D(推导) | 同上 |
| ″ | ″ | ESR @谐振点（参考） | 5.0 mΩ（=\|Z\|min，即模型在 SRF 处 ESR） | D(推导) | 同上 |
| ″ | ″ | **SRF** | **≈1.74 MHz**（\|Z\| 最小点；电抗 X 过零点 1.85 MHz，两法一致） | **D(推导)** | 同上 |
| ″ | ″ | **ESL** | **≈0.34 nH**（由 f≥300 MHz 段 \|Z\|≈2πf·ESL 反推，300 MHz–1 GHz 区间；见下方注 ①） | **D(推导)** | 同上 |
| C16,C64 | JVJ50v10M5x5 | **tanδ (max)** | **0.14** @ 120 Hz, 20 °C（JVJ 系列 50 V 档“≤0.14”；本型号 JVJ50V10M5X5 行） | **D** | jieerrui **JVJ series** datasheet（本件 LCSC 页 `https://item.szlcsc.com/48680110.html`；存档 `cap_models/jvj_series.txt`） |
| ″ | ″ | 额定纹波电流 / 尺寸 | 18 mA @120 Hz, 105 °C；D5×L5.4 mm | D | 同上 |
| ″ | ″ | **ESR @120 Hz** | **≈18.6 Ω**（**ESR = tanδ/(2πfC) = 0.14/(2π·120·10 µF)**；此为 120 Hz 低频 ESR 上限） | **D(推导)** | 同上（公式） + 同上（tanδ 来源） |
| ″ | ″ | ESR @100 kHz / ESL | **未给**（铝电解手册仅到 tanδ@120 Hz） | **GAPS** | — |
| C1 | HGC1210R5476M250NSVK | ESR / ESL / SRF | **未给** | **GAPS** | 遂宁宏明华瓷 **HC202002** 通用系列规格书（LCSC C7432791，`https://www.lcsc.com/datasheet/C7432791.pdf`，22 页）**全篇无 ESR/ESL/SRF** |
| C49,C51,C72 | CC0603KRX7R9BB104 | 身份（确认） | **100 nF · 50 V · X7R · 0603**（Yageo General Purpose X7R 规格书；型号规则 104=100nF/9=50V 吻合） | **D** | Yageo MLCC 规格书 `https://www.lcsc.com/datasheet/C14663.pdf`（LCSC C14663）→ 正文 “General Purpose & High Cap. X7R 6.3 V to 250 V, 100 pF to 47 µF” |
| ″ | ″ | ESR / ESL / SRF（数值） | **未给数值**（手册仅提供“Impedance / ESR vs frequency”**图** Fig.7–12，无数据表；另有 D.F. 规格但无 ESR 数字） | **GAPS** | 同上 |
| C18 | 0603CG2R7C500NT | ESR / ESL | **未取到原厂手册**（厂商与手册均未确认） | **GAPS** | — |
| D1 | SS36 (LCSC C7420367) | **Vf @3 A** | **0.70 V (max)** @ IF=3.0 A, Ta=25 °C | **D** | 珠海宏嘉诚(HJC) **SS32–SS320** 数据手册 `https://www.lcsc.com/product-detail/C7420367.html`（本板件厂商，存档 `cap_models/SS36_C7420367_HJC_spec.txt`） |
| ″ | ″ | **结电容 Cj** | **300 pF (typ)** @ VR=4.0 V, 1 MHz | **D** | 同上 |
| ″ | ″ | 封装（订正） | **DO-214AC (SMA)** ⚠ 与 board-facts 所记 “SMC/DO-214AB” 不符（见回报，须由主会话改 board-facts） | **D** | 同上 |
| ″ | ″ | 其它 | IR 0.5 mA @25 °C；IFSM 80 A；RθJ-A 80 °C/W | D | 同上 |

**注 ①（CL31 ESL 的用法提醒）**：0.34 nH 由官方模型高频段 \|Z\| 斜率反推，**低于常见 1206 表列典型（约 0.5–1 nH）**——系模型法结果，用于 VHF/EMI 段时应作**敏感性**，不可当唯一真值。若仿真敏感，建议以 S2P 文件或实测交叉验证。

**注 ②（JVJ ESR 的用法提醒）**：18.6 Ω 仅为 **120 Hz** 低频 ESR；本板 buck 工作在 ~784 kHz，铝电解高频 ESR 显著更低但**手册未给**，须实测或厂家高频曲线方可引用。

---

## GAPS（抓不到 / 需人确认）

1. **JVJ50v10M5x5（SMD 铝电解）的高频 ESR 与 ESL**：本轮已由 JVJ 系列手册 tanδ=0.14@120 Hz 换算出 **ESR@120 Hz≈18.6 Ω**（见 §ESR/ESL/SRF），但**高频（~784 kHz）ESR 与 ESL 手册未给**。→ 需厂家高频曲线或**实测**（ESL 通常只能实测）。
2. **HGC1210R5476M250NSVK 的 ESR / ESL / SRF**：宏明华瓷 HC202002 规格书（22 页）**未列**。→ 需厂家阻抗模型/实测。
3. **CC0603KRX7R9BB104 / 0603CG2R7C500NT 的 ESR / ESL（数值）**：身份已确认（CC0603KRX7R9BB104 = 100 nF/50 V/X7R/0603，Yageo 规格书，**D**）；但 Yageo 手册只给阻抗/ESR-频率**图**（Fig.7–12）无数据表，0603CG2R7C500NT（2.7 pF/50 V/C0G/0603）**厂商与手册均未确认**。→ 二者的**数值型 ESR/ESL** 均需原厂模型/实测。
4. **FXL0420-1R0-M 的 SRF**：数据手册未公布（仅得 DCR 27mΩ / 6A / 7A）。→ 需原厂手册或实测。**另**：DCR/Irms/Isat 取自 JLCPCB/LCSC 规格条目（源自厂家手册元数据），级别记 **D**，但**未逐行核对原厂 PDF 表格**，建议复核。
5. **ZEMS0650-150M 的磁芯材料与 µr**：一体成型/合金电感，手册未给 µr。board-facts 现用 **µr=40（E，无依据）**。→ 需原厂手册或实测（影响涡流 |J| 线性比例，形态不变）。
6. **TPS62933F 的 SW 上升/下降时间 tr/tf**：TI SLUSEA4D **未给**该项规格。→ 仍为 board-facts 的 **E=5 ns**，**须实测替换**（VHF 段包络关键量）。

> **本轮已解决并从 GAPS 移除**：① CL31A106KBHNNNE 的 ESR@100kHz / ESL / SRF（Samsung 官方 EC 模型，见 §ESR/ESL/SRF）；② JVJ50v10M5x5 的 ESR@120 Hz（tanδ 换算）；③ SS36 的 Vf@3A=0.70 V、Cj=300 pF（本板件厂商 珠海宏嘉诚 HJC 手册）。

---

## 附：本次抓取的方法/局限

- `web_search`（Brave）**无 API key，不可用**；改用 `web_fetch` 抓 `html.duckduckgo.com/html/?q=` 做检索，再抓目标页/PDF。
- 部分厂商站点（onsemi、Mouser、TME、LCSC 商品页）为 Cloudflare/JS 渲染，`curl`/`web_fetch` 均取不到正文 → 相关器件（SS36、FXL0420 原厂 PDF）落入 GAPS。
- 成功用 `curl + pypdf` 抽取 PDF 正文的：TI **SLUSEA4D**（TPS62933 系列，51 页）、Samsung **CL31A106KBHNNNE SPEC**（38 页镜像）。
- **未改动**任何 `emi*/fea*/results*/rad*/cm_*/loss_redo/` 目录；未 git push、未动 Gitea、未对外发消息。
