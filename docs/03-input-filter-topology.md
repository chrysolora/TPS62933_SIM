# TPS62933F 板真实网表（从 epru 原理图逐脚重建）

> 方法：解析 `pourSim.epru` 的 SCH_PAGE（COMPONENT/SYMBOL/PIN/LINE + 全局网络标签），union-find 求连通性。
> 脚本：`node249:/tmp/nets4.py`。**非猜测**。日期 2026-09-24。

## 拓扑
```
24VIN --F1(保险丝)--> n_f1 --D1(SS36, 串联肖特基, A=n_f1,K=PPHV)--> PPHV   [PPHV ~= 24-Vf]
PPHV:   C17||C62||C63 (3x10uF/50V, CL31A106KBHNNNE) --> GND
        L1 (1uH, FXL0420-1R0-M) --> PPHV_OUT_FILTER
PPHV_OUT_FILTER: C73||C74 (2x10uF) --> GND ; C72 (100nF) --> GND
                 R85 (100mOhm, 1206) --> n_rc ; n_rc: C16||C64 (2x10uF/50V, JVJ50v10M5x5) --> GND   [RC 阻尼]
                 R87 (383kOhm, 0603WAF3833T5E) --> EN ; EN --> R88 (35.7kOhm) --> GND   [EN/UVLO 分压]
                 U21.VIN (pin3)
D138 (SMF24A TVS): PPHV <-> GND
U21 = TPS62933FDRLR (SOT583, 8pin)
  pin1 RT  --> R89 (26.7kOhm, 0603WAF2672T5E) --> GND        [fsw 设定; datasheet 200kHz-2.2MHz]
  pin2 EN  <-- R87/R88
  pin3 VIN = PPHV_OUT_FILTER
  pin4 GND
  pin5 SW  --> L2 (15uH, ZEMS0650-150M) --> 12V
  pin6 BST --> C51 (100nF) --> SW
  pin7 SS  --> C49 (100nF) --> GND                            [注意: 板上 pin7=SS, 非 PG]
  pin8 FB  <-- R14 (140k, 从12V) / R86 (10k, 到GND) ; C18 (2.7pF, Cff, 12V->FB)
12V: C1 (唯一输出电容; LCSC HGC1210R5476M250NSVK, 疑似 47uF/25V 待核) --> GND ; R1(44.2k)+LED9 指示灯
```

## 关键结论
- **fsw 由 RT=26.7kΩ 决定 → ≈784 kHz**（旧仿真 RT 悬空 → 497kHz 错误）
- 输入滤波 = π(3x10uF / L1 1uH / 2x10uF+100nF) **+ RC 阻尼(R85 100mΩ + C16||C64 20uF)**
- **输出电容只有 C1**，不是 5x10uF+10uF=60uF（10uF 全在输入侧）
- 模型变体：仿真用 TPS62933P(PG/FSS)，板为 TPS62933F(SS) — 存在差异

## 旧版阶段1网表的错误
1. RT 未接电阻（悬空）→ fsw 错误
2. 输入滤波未建模（含 R85/C16/C64 阻尼）
3. Cout 假设 60uF（实为 C1 单颗）
4. 引脚/变体：PG vs SS
5. D1 串联压降未建模
6. 满载效率 103%（模型不守恒，无效）
