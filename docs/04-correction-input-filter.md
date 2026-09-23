# 更正：输入 EMI 滤波拓扑（2026-09-24）

> 用户指出报告 §3/§6 的输入滤波描述有误。经**从 epru 原理图逐脚连通性重建**（脚本 union-find：LINE 线段 + 元件引脚 + 全局网络标签），确认如下。

## 实际拓扑（确认）
```
24VIN --F1(保险丝)--> PPHV
   PPHV:  C17,C62,C63 = 10uF/50V (SAMSUNG CL31A106KBHNNNE, 1206) -->  GND    ← 输入电容组(电感前)
        L1 = 1uH   [PPHV <-> PPHV_OUT_FILTER]         ← π 滤波电感
   PPHV_OUT_FILTER: C73,C74 = 10uF --> GND ;  C72 = 100nF --> GND            ← 滤波后电容组
        --> U21(VIN) = TPS62933FDRLR
```

## 报告中的错误
- 报告称输入滤波 = `C1 47uF + L1 1uH + C51 100nF`。**错误**。
- 实测：**C1 = 12V 输出侧电容**（C1: 12V↔GND），非输入电容。
- 输入电容实为 **10uF ×3（C17/C62/C63）**，位于 **PPHV**（L1 之前）。
- L1(1uH) 是 **PPHV↔PPHV_OUT_FILTER** 之间的 π 滤波电感，两侧各有一组电容。

## 连通性要点（重建结果，节选）
| 器件 | 脚1 | 脚2 | 备注 |
|---|---|---|---|
| F1 | 24VIN | (PPHV) | 保险丝 |
| C17/C62/C63 | PPHV | GND | 3×10uF 输入电容 |
| L1 | PPHV_OUT_FILTER | PPHV | 1uH 滤波电感 |
| C73/C74 | PPHV_OUT_FILTER | GND | 2×10uF |
| C72 | PPHV_OUT_FILTER | GND | 100nF |
| U21 | VIN=PPHV_OUT_FILTER | GND | TPS62933F |
| C1 | 12V | GND | 输出侧电容 |
| L2 | 12V | (SW) | 15uH 主电感 |
| R14 | (FB) | 12V | 反馈上 |

（C16/C64 两颗 10uF 同区，一端接 GND，另一端未解析干净，可能并于 PPHV 一组。）

## 根因与修正
- 根因：初版网表**按元件值+TI参考设计臆断拓扑**，未做原理图级连通性验证。
- 修正：EMI 评估须用上述**真实 π 滤波**重算。CISPR 频域评估时，滤波网络 = [C=30uF] - [L=1uH] - [C=20uF+100nF] (+ LISN 50uH/50Ω)。
