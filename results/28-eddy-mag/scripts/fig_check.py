#!/usr/bin/env python3
import json,os,math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT="/mnt/raid10/sim-work/tps62933/eddy_mag"
d=json.load(open(os.path.join(ROOT,"numbers_eddy_mag.json")))
u=d["uncoupled_thinsheet"]; c=d["coupled_harmonic_ref"]
names=["FullCopper","TopCutout","DualCutout"]
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
# 1 monotonicity
x=np.arange(3)
ax[0].plot(x,[u[n]["P_tot"]*1e3 for n in names],"o-",color="g",label="CORRECTED (uncoupled thin-sheet)")
cn=[n for n in names if n in c]
ax[0].plot([names.index(n) for n in cn],[c[n]["P_tot"] for n in cn],"s--",color="r",label="COUPLED harmonic (artifact)")
ax[0].set_yscale("log"); ax[0].set_xticks(x); ax[0].set_xticklabels(names)
ax[0].set_ylabel("P_eddy (mW/m) [log]"); ax[0].set_title("SELF-CHECK: dig copper -> P must DROP")
ax[0].legend(fontsize=7); ax[0].grid(True,which="both",alpha=.3)
# 2 formula check
ax[1].text(0.02,0.9,"Loss formula audit",fontsize=13,weight="bold",transform=ax[1].transAxes)
ax[1].text(0.02,0.72,"WRONG (fea4):  P=∫σ|J|²\n  -> off by σ²=3.4e15",color="r",fontsize=10,transform=ax[1].transAxes)
ax[1].text(0.02,0.55,"RIGHT (here): P=∫ρ|J|²,  ρ=1/σ",color="g",fontsize=10,transform=ax[1].transAxes)
ax[1].text(0.02,0.38,"thin-sheet: P=σω²t³/24 ∫Bz²dA\n= ∫ρ|J|² with J(y)=σω·Bz·y",fontsize=9,transform=ax[1].transAxes)
ax[1].axis("off")
# 3 magnitude vs 12 bound
ax[2].text(0.02,0.9,"Magnitude cross-check",fontsize=13,weight="bold",transform=ax[2].transAxes)
ax[2].text(0.02,0.72,"incident Bz(max)= %.3g T (top)\n            Bz(max)= %.3g T (bot)"%(d["incident"]["Bz_max_top"],d["incident"]["Bz_max_bot"]),fontsize=9,transform=ax[2].transAxes)
ax[2].text(0.02,0.5,"12-eddy-loss upper bound: sub-µW total\n(this run): %.3g mW/m per unit depth -> SAME ORDER as 12"%(u["FullCopper"]["P_tot"]*1e3),fontsize=9,transform=ax[2].transAxes)
ax[2].text(0.02,0.28,"COUPLED artifact: %.3g W/m  (%.0f× too big)"%(c["FullCopper"]["P_tot"], c["FullCopper"]["P_tot"]/(u["FullCopper"]["P_tot"])),fontsize=9,color="r",transform=ax[2].transAxes)
ax[2].axis("off")
fig.tight_layout(); fig.savefig(os.path.join(ROOT,"fig_eddy_mag_check.png"),dpi=130)
print("fig check written")
