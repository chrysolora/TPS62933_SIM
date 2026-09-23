#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""with vs without pi filter, geometry v3 params (standalone)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
FSW=805e3
def Zc(f,C): return 1.0/(1j*2*np.pi*f*C)
def Zl(f,L): return 1j*2*np.pi*f*L
def capZ(f,C,esr,esl): return esr+Zl(f,esl)+Zc(f,C)
def trap(fsw,ipk,duty,tr,tf,kmax,nper=8,ppp=16384):
    T=1.0/fsw; N=nper*ppp; dt=(nper*T)/N
    t=np.arange(N)*dt; tt=t%T; on=duty*T; x=np.zeros_like(t)
    if tr>0:
        r=tt<np.minimum(tr,on); x[r]=ipk*tt[r]/tr
    fl=(tt>=tr)&(tt<on); x[fl]=ipk
    if tf>0:
        fa=(tt>=on)&(tt<on+tf); x[fa]=ipk*(1-(tt[fa]-on)/tf)
    x-=x.mean(); w=np.hanning(N); X=np.fft.rfft(x*w); f=np.fft.rfftfreq(N,dt)
    amp=2*np.abs(X)/w.sum(); ks=np.arange(1,kmax+1); fh=ks*fsw
    idx=np.array([int(round(fr/f[1])) for fr in fh]); return fh,amp[idx]
def solve(n,branches,inj):
    Y=np.zeros((n,n),complex); b=np.zeros(n,complex)
    for (i,j,Z) in branches:
        if abs(Z)<1e-15: Z=1e-12
        y=1.0/Z
        if i>=0: Y[i,i]+=y
        if j>=0: Y[j,j]+=y
        if i>=0 and j>=0: Y[i,j]-=y; Y[j,i]-=y
    b[inj]+=1.0; return np.linalg.solve(Y,b)
def cispr_B(f):
    f=np.asarray(f,float)
    return np.where(f<0.5e6,66.0-10.0*(np.log10(f)-np.log10(0.15e6))/(np.log10(0.5e6)-np.log10(0.15e6)),np.where(f<5e6,56.0,60.0))
def cispr_A(f):
    f=np.asarray(f,float)
    return np.where(f<0.5e6,79.0-6.0*(np.log10(f)-np.log10(0.15e6))/(np.log10(0.5e6)-np.log10(0.15e6)),np.where(f<5e6,73.0,73.0))
P=dict(C10u=10e-6,C100n=100e-9,Ce=10e-6,L1=1e-6,dcr_l1=25e-3,cp_l1=4e-12,rs_dm=50.0,c_lisn=1e-6,l_mains=50e-6,
 esr_10u=2e-3,esl_10u=1.0e-9,esr_100n=20e-3,esl_100n=0.8e-9,esr_e=4.0,esl_e=3.0e-9,r85=0.1,
 lf1=6.32e-9,ltr_loop=7.88e-9,lhot=0.83e-9,ltr_capA=0.75e-9,rf1=5e-3)
def damping(fr):
    zc=capZ(fr,P["Ce"],P["esr_e"],P["esl_e"]); return P["r85"]+1.0/(1.0/zc+1.0/zc)
def Zt(fr,nofilter=False):
    b=[(0,-1,P["rs_dm"]+Zc(fr,P["c_lisn"])),(0,-1,Zl(fr,P["l_mains"])),(0,1,P["rf1"]+Zl(fr,P["lf1"]))]
    if nofilter:
        b.append((1,2,1e-9))
    else:
        for _ in range(3): b.append((1,-1,capZ(fr,P["C10u"],P["esr_10u"],P["esl_10u"]+P["ltr_capA"])))
        b.append((1,2,P["dcr_l1"]+Zl(fr,P["L1"]+P["ltr_loop"])))
        if P["cp_l1"]>0: b.append((1,2,Zc(fr,P["cp_l1"])))
    b.append((2,3,Zl(fr,P["lhot"])))
    for _ in range(2): b.append((3,-1,capZ(fr,P["C10u"],P["esr_10u"],P["esl_10u"])))
    b.append((3,-1,capZ(fr,P["C100n"],P["esr_100n"],P["esl_100n"])))
    b.append((3,-1,damping(fr)))
    V=solve(4,b,2); return abs(V[0])
KMAX=int(30e6/FSW)
fg=np.logspace(np.log10(150e3),np.log10(30e6),4000)
fig,ax=plt.subplots(figsize=(10,6.8))
ax.plot(fg/1e6,cispr_B(fg),"k-",lw=2.0,label="CISPR 32 Class B")
ax.plot(fg/1e6,cispr_A(fg),"k--",lw=1.2,label="CISPR 32 Class A")
f,a=trap(FSW,3.0,0.5,10e-9,10e-9,KMAX)
for nf,ls,col,lab in [(False,"-","#1f77b4","with \u03c0 filter"),(True,"--","#d62728","without \u03c0 filter")]:
    z=np.array([Zt(x,nf) for x in f]); db=20*np.log10(np.maximum(a*z,1e-30)/1e-6)
    ax.semilogx(f/1e6,db,ls,color=col,marker=".",ms=4,label="full load "+lab)
    mg=cispr_B(f)-db; print("%s: |Zt|@805k=%.3g ohm V1=%.1f dBuV worst_margin=%.1f dB @ %.2f MHz"%(lab,z[0],db[0],mg.min(),f[mg.argmin()]/1e6))
ax.set_xlabel("Frequency [MHz]"); ax.set_ylabel("DM noise on 50\u03a9 LISN  [dB\u00b5V]")
ax.set_title("TPS62933 DM conducted EMI (geometry v3) \u2014 with vs without \u03c0 filter (full load)")
ax.set_ylim(-20,110); ax.grid(True,which="both",alpha=0.3); ax.legend(fontsize=9,loc="upper right")
fig.tight_layout(); fig.savefig("/mnt/raid10/sim-work/tps62933/emi_v3/fig_emi_cispr_v3_nofilter.png",dpi=140)
print("saved")
