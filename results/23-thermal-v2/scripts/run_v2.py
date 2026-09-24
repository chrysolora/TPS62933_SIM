#!/usr/bin/env python3
"""TPS62933 board thermal v2 driver.
- builds meshes (3 refinement levels)
- ElmerGrid convert + BC auto-detect
- h(DeltaT) natural-convection outer iteration (main case)
- runs convergence + sensitivity cases
- collects numbers/convergence json
Usage: run_v2.py <stage>   stage in {mesh, setup, hconv, run, all}
"""
import os, sys, json, glob, subprocess, shutil, re
import numpy as np

ROOT = "/mnt/raid10/sim-work/tps62933/thermal_v2"
os.chdir(ROOT)
MM = 1e-3
Lx, Ly, Lz = 26.50, 50.50, 1.49   # mm
Lc = (Lx*Ly*1e-6) / (2*(Lx+Ly)*1e-3)  # char length m = area/perimeter

MESHES = {"L1": (2.5, 0.60), "L2": (1.6, 0.35), "L3": (1.0, 0.18)}


# ---------- air properties / natural convection ----------
_T = np.array([250, 300, 350, 400, 450.0])
_NU = np.array([11.44e-6, 15.89e-6, 20.92e-6, 26.41e-6, 32.39e-6])
_K = np.array([0.0223, 0.0263, 0.0300, 0.0338, 0.0373])
_PR = np.array([0.720, 0.707, 0.700, 0.690, 0.686])
G = 9.81


def h_natural(dT, Tamb_C):
    """horizontal-plate natural convection correlation (Incropera), E-level.
    Nu=0.54 Ra^0.25 (laminar) / 0.15 Ra^(1/3) (turbulent); L=area/perimeter."""
    if dT <= 0.5:
        return 5.0
    Tf = Tamb_C + dT/2.0 + 273.15
    Tf = min(max(Tf, 250), 450)
    k = np.interp(Tf, _T, _K)
    nu = np.interp(Tf, _T, _NU)
    pr = np.interp(Tf, _T, _PR)
    beta = 1.0 / Tf
    Ra = G * beta * dT * Lc**3 / nu**2 * pr
    if Ra < 1e7:
        Nu = 0.54 * Ra**0.25
    elif Ra < 1e11:
        Nu = 0.15 * Ra**(1.0/3.0)
    else:
        Nu = 0.15 * Ra**(1.0/3.0)
    return float(Nu * k / Lc)


# ---------- mesh build ----------
def build_meshes():
    for name, (lmax, lmin) in MESHES.items():
        if os.path.exists(f"mesh_{name}.msh"):
            print("exists", name); continue
        print("building", name, flush=True)
        r = subprocess.run(["python3", "build_mesh_v2.py", str(lmax), str(lmin), f"mesh_{name}"],
                           capture_output=True, text=True, timeout=900)
        if r.returncode != 0:
            print(r.stdout[-2000:], r.stderr[-2000:]); raise SystemExit("mesh fail " + name)
        print(name, "done", flush=True)


# ---------- ElmerGrid + BC detect ----------
def setup_elmer(name):
    src = os.path.abspath(f"mesh_{name}.msh")
    out = f"elmer_{name}"
    if not os.path.isdir(os.path.join(out, "mesh")):
        os.makedirs(out, exist_ok=True)
        r = subprocess.run(["ElmerGrid", "14", "2", src, "-out", "mesh", "-autoclean"],
                           cwd=out, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            print(r.stdout[-1500:], r.stderr[-1500:]); raise SystemExit("elmgrid fail")
    return detect_bc(os.path.join(out, "mesh"))


def detect_bc(mdir):
    N = {}
    for ln in open(os.path.join(mdir, "mesh.nodes")):
        p = ln.split()
        if len(p) < 5:
            continue
        try:
            N[int(p[0])] = (float(p[2]), float(p[3]), float(p[4]))
        except ValueError:
            pass
    groups = {}
    for ln in open(os.path.join(mdir, "mesh.boundary")):
        p = ln.split()
        if len(p) < 6:
            continue
        try:
            bidx = int(p[1]); nodes = [int(t) for t in p[5:]]
        except ValueError:
            continue
        groups.setdefault(bidx, []).extend(nodes)
    # ---- body detection (ElmerGrid renumbers bodies to 1..N) ----
    comps_xyz = {"U21": (8.25e-3, 26.88e-3), "D1": (14.10e-3, 6.56e-3),
                 "L2": (6.22e-3, 33.74e-3), "L1": (16.38e-3, 18.75e-3)}
    bnodes = {}
    for ln in open(os.path.join(mdir, "mesh.elements")):
        p = ln.split()
        if len(p) < 5:
            continue
        try:
            b = int(p[1]); nds = [int(t) for t in p[3:]]
        except ValueError:
            continue
        bnodes.setdefault(b, []).extend(nds)
    bcen = {}
    for b, nds in bnodes.items():
        pts = np.array([N[n] for n in nds if n in N])
        bcen[b] = pts.mean(axis=0) if len(pts) else None
    board_b = min([b for b in bcen if bcen[b] is not None], key=lambda b: abs(bcen[b][2]))
    bmap = {}
    for nm, (tx, ty) in comps_xyz.items():
        best, bd = None, 1e9
        for b, c in bcen.items():
            if b == board_b or c is None:
                continue
            d = (c[0]-tx)**2 + (c[1]-ty)**2
            if d < bd:
                bd, best = d, b
        bmap[nm] = best
    info = {}
    for k, ns in groups.items():
        pts = np.array([N[n] for n in ns if n in N])
        info[k] = {"n": len(ns), "cen": pts.mean(axis=0) if len(pts) else None,
                   "zmin": pts[:, 2].min() if len(pts) else None,
                   "zmax": pts[:, 2].max() if len(pts) else None}
    zmax_all = max(v["zmax"] for v in info.values())
    zmin_all = min(v["zmin"] for v in info.values())
    outer = max(info, key=lambda k: info[k]["n"])
    # terminals: small groups centred near terminal xy, top face
    terms = {}
    for nm, (tx, ty) in [("TERM_IN", (4.95e-3, 46.69e-3)), ("TERM_OUT", (5.00e-3, 4.76e-3))]:
        best, bd = None, 1e9
        for k, v in info.items():
            if k == outer or v["cen"] is None:
                continue
            d = (v["cen"][0]-tx)**2 + (v["cen"][1]-ty)**2
            if d < bd:
                bd, best = d, k
        terms[nm] = best
    return {"outer": outer, "term_in": terms["TERM_IN"], "term_out": terms["TERM_OUT"],
            "board_body": int(board_b), "bodymap": {k: int(v) for k, v in bmap.items()},
            "n_outer": info[outer]["n"], "zrange": [float(zmin_all), float(zmax_all)]}


def link_mesh(name, casedir):
    """copy (hardlink) mesh into casedir so results are not shared."""
    os.makedirs(casedir, exist_ok=True)
    m = os.path.join(casedir, "mesh")
    if os.path.islink(m):
        os.remove(m)
    src = os.path.abspath(os.path.join(f"elmer_{name}", "mesh"))
    if os.path.isdir(m):
        shutil.rmtree(m)
    shutil.copytree(src, m, copy_function=os.link)
    for f in glob.glob(os.path.join(m, "*.vtu")):
        os.remove(f)


# ---------- solve ----------
def solve(casedir):
    r = subprocess.run(["ElmerSolver", "case.sif"], cwd=casedir, capture_output=True,
                       text=True, timeout=1800, env={**os.environ, "OMP_NUM_THREADS": "2"})
    return r.returncode, r.stdout[-3000:]


def read_vtu(casedir):
    import re as _re
    fs = glob.glob(os.path.join(casedir, "mesh", "*.vtu"))
    def num(f):
        m = _re.search(r"_t(\d+)\.vtu$", f); return int(m.group(1)) if m else -1
    f = sorted(fs, key=num, reverse=True)[0]
    raw = open(f, "rb").read(); txt = raw.decode("latin1")
    npts = int(_re.search(r'NumberOfPoints="(\d+)"', txt).group(1))
    i = txt.index("<AppendedData"); base = txt.index("_", i) + 1
    def arr(off, dt, n):
        p = base + off
        ln = int(np.frombuffer(raw[p:p+4], dtype="<u4")[0])
        return np.frombuffer(raw[p+4:p+4+ln], dtype=dt)[:n]
    toff = int(_re.search(r'Name="temperature"[^>]*offset="(\d+)"', txt).group(1))
    poff = int(_re.search(r'<Points>\s*<DataArray[^>]*offset="(\d+)"', txt, re.S).group(1))
    T = arr(toff, "<f8", npts)
    P = arr(poff, "<f8", npts*3).reshape(-1, 3)
    return P, T


def surf_mean_T(casedir, mdir, bcidx):
    """area-ish mean surface temp over a boundary index (node-average)."""
    N = {}
    for ln in open(os.path.join(mdir, "mesh.nodes")):
        p = ln.split()
        if len(p) < 5:
            continue
        try:
            N[int(p[0])] = (round(float(p[2]), 9), round(float(p[3]), 9), round(float(p[4]), 9))
        except ValueError:
            pass
    ns = []
    for ln in open(os.path.join(mdir, "mesh.boundary")):
        p = ln.split()
        if len(p) < 6:
            continue
        try:
            if int(p[1]) != bcidx:
                continue
            ns += [int(t) for t in p[5:]]
        except ValueError:
            continue
    P, T = read_vtu(casedir)
    cmap = {}
    for j in range(len(P)):
        key = (round(P[j, 0], 9), round(P[j, 1], 9), round(P[j, 2], 9))
        cmap[key] = j
    vals = [T[cmap[N[n]]] for n in ns if n in N and (N[n] in cmap)]
    return float(np.mean(vals)) if vals else None


# ---------- orchestration ----------
def build_case_to(casedir, mesh, Iout, h, Tamb, fin, rad, thscale=1.0):
    bcm = json.load(open(f"bc_{mesh}.json"))
    env = {**os.environ, "THETA_SCALE": str(thscale),
           "BODY_BOARD": str(bcm["board_body"])}
    for k, v in bcm["bodymap"].items():
        env["BODY_" + k] = str(v)
    args = ["python3", "make_case_v2.py", str(Iout), str(h), str(Tamb), casedir,
            str(bcm["outer"]), str(bcm["term_in"]), str(bcm["term_out"]),
            "fin_on" if fin else "fin_off", ("nord" if not rad else "0.9")]
    r = subprocess.run(args, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(r.stdout, r.stderr); raise SystemExit("make_case fail " + casedir)
    link_mesh(mesh, casedir)
    return json.load(open(os.path.join(casedir, "heat_sources.json")))


def h_iterate(casedir, mesh, Iout, Tamb, fin=True, rad=True, h0=10.0, thscale=1.0):
    bcm = json.load(open(f"bc_{mesh}.json"))
    mdir = os.path.join(f"elmer_{mesh}", "mesh")
    h = h0; hist = []
    for it in range(5):
        hs = build_case_to(casedir, mesh, Iout, h, Tamb, fin, rad, thscale)
        rc, log = solve(casedir)
        if rc != 0:
            print("SOLVE FAIL", log[-1500:]); raise SystemExit("solve fail " + casedir)
        Ts = surf_mean_T(casedir, mdir, bcm["outer"])
        dT = Ts - (Tamb + 273.15)
        hnew = h_natural(dT, Tamb)
        hist.append({"it": it, "h": round(h, 3), "Tsurf_avg_C": round(Ts-273.15, 2),
                     "dT": round(dT, 2), "h_next": round(hnew, 3)})
        print("  h-it", it, "h=%.2f" % h, "Tsurf=%.1f" % (Ts-273.15), "dT=%.1f" % dT,
              "h_new=%.2f" % hnew, flush=True)
        if abs(hnew - h) / h < 0.02 or it >= 3:
            h = hnew; break
        h = 0.5*h + 0.5*hnew
    build_case_to(casedir, mesh, Iout, h, Tamb, fin, rad, thscale)
    rc, log = solve(casedir)
    if rc != 0:
        raise SystemExit("final solve fail " + casedir)
    return h, hist


def metrics(casedir, comps):
    P, T = read_vtu(casedir)
    T = T - 273.15
    x, y, z = P[:, 0]*1e3, P[:, 1]*1e3, P[:, 2]*1e3
    rec = {"Tmax_C": round(float(T.max()), 1),
           "Tmax_xyz_mm": [round(float(x[np.argmax(T)]), 2), round(float(y[np.argmax(T)]), 2),
                           round(float(z[np.argmax(T)]), 2)],
           "Tmin_C": round(float(T.min()), 1)}
    for nm, (w, l, cx, cy, hh) in comps.items():
        s = (np.abs(x-cx) <= w/2) & (np.abs(y-cy) <= l/2) & (z > 0.745 + 0.001)
        rec[nm+"_body_max_C"] = round(float(T[s].max()), 1) if s.sum() else None
        rec[nm+"_body_avg_C"] = round(float(T[s].mean()), 1) if s.sum() else None
    # board-interface (z = +0.745) mean under each device
    for nm, (w, l, cx, cy, hh) in comps.items():
        s = (np.abs(x-cx) <= w/2) & (np.abs(y-cy) <= l/2) & (np.abs(z - 0.745) < 0.002)
        rec[nm+"_iface_C"] = round(float(T[s].mean()), 1) if s.sum() else None
    return rec, P, T


COMPS = {"U21": (2.10, 1.60, 8.25, 26.88, 0.60), "D1": (4.30, 2.60, 14.10, 6.56, 2.20),
         "L2": (6.50, 6.50, 6.22, 33.74, 3.00), "L1": (4.20, 4.20, 16.38, 18.75, 2.00)}


def conservation(casedir, Tamb, eps=0.9):
    """Q_out by convection+radiation vs Sum P (from heat_sources)."""
    hs = json.load(open(os.path.join(casedir, "heat_sources.json")))
    Ptot = hs["total_W"]
    h = hs["h_W_per_m2K"]; eps = hs["emissivity"]
    # rough from scalar: use Tmax of loaded case; exact integral needs surface;
    # approximate Q_out = h*A*(Tsurf-Tamb)+eps*sig*A*(Tsurf^4-Tamb^4) using board+comp outer area
    return {"Ptot_W": Ptot}


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("mesh", "all"):
        build_meshes()
    if stage in ("setup", "all"):
        for n in MESHES:
            bcm = setup_elmer(n)
            json.dump(bcm, open(f"bc_{n}.json", "w"), indent=1)
            print("BC", n, bcm, flush=True)


def run_one(name, mesh, Iout, Tamb, fin=True, rad=True, hmode="conv", hval=5.0,
            thscale=1.0, do_conv=False):
    cd = os.path.join(ROOT, "cases", name)
    os.makedirs(cd, exist_ok=True)
    if hmode == "conv":
        h, hist = h_iterate(cd, mesh, Iout, Tamb, fin, rad, thscale=thscale)
    else:
        h = hval
        build_case_to(cd, mesh, Iout, h, Tamb, fin, rad, thscale)
        rc, log = solve(cd)
        if rc != 0:
            print("FAIL", name, log[-1200:]); return None
        hist = []
    rec, P, T = metrics(cd, COMPS)
    hs = json.load(open(os.path.join(cd, "heat_sources.json")))
    rec.update({"mesh": mesh, "Iout_A": Iout, "Tamb_C": Tamb, "fin": fin, "rad": rad,
                "h_final": h, "h_hist": hist, "total_W": hs["total_W"],
                "powers_W": hs["powers_W"], "theta_scale": thscale,
                "comp_k": hs["comp_effective_k_W_per_mK"]})
    # conservation: radiative+convective flux integral over outer boundary
    try:
        bcm = json.load(open(f"bc_{mesh}.json"))
        mdir = os.path.join(f"elmer_{mesh}", "mesh")
        rec["Tsurf_avg_C"] = round(surf_mean_T(cd, mdir, bcm["outer"]) - 273.15, 2)
    except Exception:
        pass
    json.dump(rec, open(os.path.join(cd, "result.json"), "w"), indent=1)
    print(name, rec, flush=True)
    return rec


def run_all():
    res = {}
    # convergence: 3 mesh levels, full load
    for lvl in ["L1", "L2", "L3"]:
        r = run_one(f"CONV_{lvl}", lvl, 3.0, 25, fin=True, rad=True, hmode="conv")
        res[f"CONV_{lvl}"] = r
        json.dump(res, open(os.path.join(ROOT, "numbers_thermal_v2.json"), "w"), indent=1)
    # main map case (reuse CONV_L2)
    res["MAIN_full_hconv_T25"] = res["CONV_L2"]
    # sensitivity on L2
    sens = [
        ("S_h5fix",  3.0, 25, True, True,  "fix", 5.0, 1.0),
        ("S_h10fix", 3.0, 25, True, True,  "fix", 10.0, 1.0),
        ("S_T40",    3.0, 40, True, True,  "conv", 5.0, 1.0),
        ("S_fin_off",3.0, 25, False, True, "conv", 5.0, 1.0),
        ("S_half",   1.5, 25, True, True,  "conv", 5.0, 1.0),
        ("S_light",  0.6, 25, True, True,  "conv", 5.0, 1.0),
        ("S_noRad",  3.0, 25, True, False, "conv", 5.0, 1.0),
        ("S_theta15",3.0, 25, True, True,  "conv", 5.0, 1.5),
    ]
    for (nm, I, T, fin, rad, hm, hv, ts) in sens:
        r = run_one(nm, "L2", I, T, fin=fin, rad=rad, hmode=hm, hval=hv, thscale=ts)
        res[nm] = r
        json.dump(res, open(os.path.join(ROOT, "numbers_thermal_v2.json"), "w"), indent=1)
    # convergence json
    conv = {}
    for lvl in ["L1", "L2", "L3"]:
        if res.get(f"CONV_{lvl}"):
            conv[lvl] = {k: res[f"CONV_{lvl}"][k] for k in
                         ("Tmax_C", "Tmax_xyz_mm", "U21_body_max_C", "U21_iface_C",
                          "D1_body_max_C", "L2_body_max_C", "L1_body_max_C", "h_final")}
    meta = {}
    for lvl in ["L1", "L2", "L3"]:
        m = json.load(open(f"mesh_{lvl}_meta.json"))
        meta[lvl] = {"tets": m["tets"], "lcmin_mm": m["lcmin_mm"], "lcmax_mm": m["lcmax_mm"]}
    json.dump({"levels": conv, "mesh": meta}, open(os.path.join(ROOT, "convergence_v2.json"), "w"), indent=1)
    print("ALL DONE", flush=True)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "run":
    run_all()

