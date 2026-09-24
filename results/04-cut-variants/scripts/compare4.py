import json, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Noto Sans CJK JP", "WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
import matplotlib.image as mpimg

ROOT = sys.argv[1]
VAR = [("fullcu", "FullCopper (铜皮完整)"),
       ("topcut", "TopCutout (顶层挖空)"),
       ("dualcut", "DualCutout (顶+底层挖空)")]

fig, axes = plt.subplots(1, 3, figsize=(12, 5.6))
for ax, (v, t) in zip(axes, VAR):
    p = os.path.join(ROOT, v, "fig_field_top_%s.png" % v)
    ax.imshow(mpimg.imread(p)); ax.set_title(t, fontsize=11); ax.axis("off")
fig.suptitle("L2 下方铺铜涡流 |J| 顶视对比（同色标 log 0.1–3e7 A/m^2）", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "compare_top.png"), dpi=130)
print("wrote compare_top.png")

fig2, ax2 = plt.subplots(figsize=(8, 6))
for v, t in VAR:
    m = json.load(open(os.path.join(ROOT, v, "metrics.json")))
    b = m["comb"]["bins"]
    r = [0.5 * (lo + hi) for lo, hi, n, mean, mx in b]
    y = [max(mean, 1e-1) for lo, hi, n, mean, mx in b]
    ax2.semilogy(r, y, "o-", label=t)
ax2.axhline(1e7, color="r", ls=":", lw=1, label="1e7")
ax2.set_xlabel("半径 r (mm, 起点=电感 L2 中心)")
ax2.set_ylabel("平均 |J| (A/m^2)")
ax2.set_title("三版径向 |J| 剖面（顶+底层铜皮合并）")
ax2.grid(True, which="both", alpha=0.3)
ax2.legend(fontsize=9)
fig2.tight_layout()
fig2.savefig(os.path.join(ROOT, "compare_radial.png"), dpi=130)
print("wrote compare_radial.png")
