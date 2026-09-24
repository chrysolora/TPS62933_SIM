#!/usr/bin/env python3
"""run_3x3_real.py — TPS62933 board thermal 3x3 driver, REAL cutout polygons.

Identical physics/params/boundary/h-chain/algorithm as thermal_3x3/run_3x3.py.
ONLY the keep-out GEOMETRY differs: true source polygons (topcut != dualcut)
instead of the bounding rectangles.
Usage: run_3x3_real.py [mesh|all]
"""
import os, sys, json, glob, subprocess, shutil, re
import numpy as np

ROOT = "/mnt/raid10/sim-work/tps62933/thermal_3x3_real"
os.chdir(ROOT)
Lx, Ly, Lz = 26.50e-3, 50.50e-3, 1.49e-3
Lc = (26.50 * 50.50 * 1e-6) / (2 * (26.50 + 50.50) * 1e-3)
KZ_VIA = 1.3612
GEOS = ["fullcu", "topcut", "dualcut"]
LOADS = [3.0, 1.5, 0.6]
H0 = {3.0: 11.8, 1.5: 9.79, 0.6: 7.78}
COMPS = {"U21": (2.10, 1.60, 8.25, 26.88), "D1": (4.30, 2.60, 14.10, 6.56),
         "L2": (6.50, 6.50, 6.22, 33.74), "L1": (4.20, 4.20, 16.38, 18.75)}
MESHES = {"main": (1.6, 0.35), "coarse": (2.5, 0.60)}

_pj = json.load(open(os.path.join(ROOT, "cutout_polygons.json")))
POLY = {"topcut": {}, "dualcut": {}}
for g in ("topcut", "dualcut"):
    for lbl, rec in _pj["polygons"][g].items():
        POLY[g][lbl] = rec["vertices_mesh_mm"]


def point_in_poly(x, y, poly):
    inside = False; n = len(poly); j = n - 1
    for i in range(n):
        xi, yi = poly[i]; xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
            inside = not inside
        j = i
    return inside


_T = np.array([250, 300, 350, 400, 450.0])
_NU = np.array([11.44e-6, 15.89e-6, 20.92e-6, 26.41e-6, 32.39e-6])
_K = np.array([0.0223, 0.0263, 0.0300, 0.0338, 0.0373])
_PR = np.array([0.720, 0.707, 0.700, 0.690, 0.686])
G = 9.81


def h_natural(dT, Tamb_C):
    if dT <= 0.5:
        return 5.0
    Tf = min(max(Tamb_C + dT / 2.0 + 273.15, 250), 450)
    k = np.interp(Tf, _T, _K); nu = np.interp(Tf, _T, _NU); pr = np.interp(Tf, _T, _PR)
    beta = 1.0 / Tf
    Ra = G * beta * dT * Lc ** 3 / nu ** 2 * pr
    Nu = 0.54 * Ra ** 0.25 if Ra < 1e7 else 0.15 * Ra ** (1.0 / 3.0)
    return float(Nu * k / Lc)


def build_mesh(geo, level):
    lmax, lmin = MESHES[level]
    out = f"mesh_{geo}_{level}"
    if os.path.exists(out + ".msh"):
        return
    print("building mesh", out, flush=True)
    r = subprocess.run(["python3", "build_mesh_3x3_real.py", str(lmax), str(lmin), geo, out],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:]); raise SystemExit("mesh fail " + out)


def setup_elmer(geo, level):
    out = f"elmer_{geo}_{level}"
    src = os.path.abspath(f"mesh_{geo}_{level}.msh")
    if not os.path.isdir(os.path.join(out, "mesh")):
        os.makedirs(out, exist_ok=True)
        r = subprocess.run(["ElmerGrid", "14", "2", src, "-out", "mesh", "-autoclean"],
                           cwd=out, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            print(r.stdout[-1500:], r.stderr[-1500:]); raise SystemExit("elmgrid fail")
    return detect_bc(os.path.join(out, "mesh"), geo)


def detect_bc(mdir, geo):
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
            groups.setdefault(int(p[1]), []).extend([int(t) for t in p[5:]])
        except ValueError:
            pass
    bel = {}
    for ln in open(os.path.join(mdir, "mesh.elements")):
        p = ln.split()
        if len(p) < 6:
            continue
        try:
            b = int(p[1])
            if int(p[2]) == 504:
                bel.setdefault(b, []).extend([int(t) for t in p[3:]])
        except ValueError:
            pass
    bcen = {}
    for b, nds in bel.items():
        pts = np.array([N[n] for n in set(nds) if n in N])
        if len(pts):
            bcen[b] = pts.mean(axis=0)
    board_main, cutout, comps = [], [], {}
    kp = list(POLY.get(geo, {}).values())
    for b, c in bcen.items():
        X, Y = c[0] * 1e3, c[1] * 1e3
        if abs(c[2]) < Lz / 2 * 0.9:
            if geo != "fullcu" and any(point_in_poly(X, Y, p) for p in kp):
                cutout.append(int(b))
            else:
                board_main.append(int(b))
        else:
            comps[int(b)] = (X, Y)
    bmap = {}
    for nm, (w, l, cx, cy) in COMPS.items():
        best, bd = None, 1e18
        for b, (X, Y) in comps.items():
            d = (X - cx) ** 2 + (Y - cy) ** 2
            if d < bd:
                bd, best = d, b
        bmap[nm] = best
    info = {}
    for k, ns in groups.items():
        pts = np.array([N[n] for n in set(ns) if n in N])
        info[k] = {"n": len(set(ns)), "cen": pts.mean(axis=0) if len(pts) else None}
    outer = max(info, key=lambda k: info[k]["n"])
    terms = {}
    for nm, (tx, ty) in [("TERM_IN", (4.95e-3, 46.69e-3)), ("TERM_OUT", (5.00e-3, 4.76e-3))]:
        best, bd = None, 1e18
        for k, v in info.items():
            if k == outer or v["cen"] is None:
                continue
            d = (v["cen"][0] - tx) ** 2 + (v["cen"][1] - ty) ** 2
            if d < bd:
                bd, best = d, k
        terms[nm] = best
    return {"outer": outer, "term_in": int(terms["TERM_IN"]), "term_out": int(terms["TERM_OUT"]),
            "board_main": [int(b) for b in board_main], "cutout": [int(b) for b in cutout],
            "comps": {k: int(v) for k, v in bmap.items()}}


def read_vtu(casedir):
    fs = glob.glob(os.path.join(casedir, "mesh", "*.vtu")) or glob.glob(os.path.join(casedir, "*.vtu"))
    def num(f):
        m = re.search(r"_t(\d+)\.vtu$", f); return int(m.group(1)) if m else -1
    f = sorted(fs, key=num, reverse=True)[0]
    raw = open(f, "rb").read(); txt = raw.decode("latin1")
    npts = int(re.search(r'NumberOfPoints="(\d+)"', txt).group(1))
    base = txt.index("_", txt.index("<AppendedData")) + 1
    def arr(off, dt, n):
        p = base + off
        ln = int(np.frombuffer(raw[p:p + 4], dtype="<u4")[0])
        return np.frombuffer(raw[p + 4:p + 4 + ln], dtype=dt)[:n]
    toff = int(re.search(r'Name="temperature"[^>]*offset="(\d+)"', txt).group(1))
    poff = int(re.search(r'<Points>\s*<DataArray[^>]*offset="(\d+)"', txt, re.S).group(1))
    T = arr(toff, "<f8", npts); P = arr(poff, "<f8", npts * 3).reshape(-1, 3)
    return P, T


def solve(casedir):
    r = subprocess.run(["ElmerSolver", "case.sif"], cwd=casedir, capture_output=True, text=True,
                       timeout=2400, env={**os.environ, "OMP_NUM_THREADS": "2"})
    return r.returncode, r.stdout[-2500:]


def link_mesh(geo, level, casedir):
    os.makedirs(casedir, exist_ok=True)
    m = os.path.join(casedir, "mesh")
    src = os.path.abspath(os.path.join(f"elmer_{geo}_{level}", "mesh"))
    if os.path.isdir(m):
        shutil.rmtree(m)
    shutil.copytree(src, m, copy_function=os.link)
    for f in glob.glob(os.path.join(m, "*.vtu")):
        os.remove(f)


def build_case(geo, level, Iout, h, Tamb, bc, casedir, kz=KZ_VIA):
    env = {**os.environ, "KZ": str(kz), "CUT_MODE": geo}
    json.dump({"board_main": bc["board_main"], "cutout": bc["cutout"], "comps": bc["comps"]},
              open("bodies.json", "w"))
    args = ["python3", "make_case_3x3.py", str(Iout), str(h), str(Tamb), casedir,
            str(bc["outer"]), str(bc["term_in"]), str(bc["term_out"]), "bodies.json"]
    r = subprocess.run(args, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(r.stdout, r.stderr); raise SystemExit("make_case fail " + casedir)
    link_mesh(geo, level, casedir)
    return json.load(open(os.path.join(casedir, "heat_sources.json")))


def surf_mean_T(casedir, mdir, bcidx):
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
            if int(p[1]) == bcidx:
                ns += [int(t) for t in p[5:]]
        except ValueError:
            continue
    P, T = read_vtu(casedir)
    cmap = {(round(P[j, 0], 9), round(P[j, 1], 9), round(P[j, 2], 9)): j for j in range(len(P))}
    vals = [T[cmap[N[n]]] for n in ns if n in N and N[n] in cmap]
    return float(np.mean(vals)) if vals else None


def metrics(casedir):
    P, T = read_vtu(casedir); T = T - 273.15
    x, y, z = P[:, 0] * 1e3, P[:, 1] * 1e3, P[:, 2] * 1e3
    board = np.abs(z) <= 0.7455
    rec = {"Tmax_body_C": round(float(T.max()), 2), "Tmin_C": round(float(T.min()), 2),
           "board_mean_C": round(float(T[board].mean()), 2),
           "board_max_C": round(float(T[board].max()), 2)}
    for nm, (w, l, cx, cy) in COMPS.items():
        s = (np.abs(x - cx) <= w / 2) & (np.abs(y - cy) <= l / 2) & (z > 0.7455)
        rec[nm + "_body_max_C"] = round(float(T[s].max()), 2) if s.sum() else None
        rec[nm + "_body_avg_C"] = round(float(T[s].mean()), 2) if s.sum() else None
        s2 = (np.abs(x - cx) <= w / 2) & (np.abs(y - cy) <= l / 2) & (np.abs(z - 0.745) < 0.003)
        rec[nm + "_iface_C"] = round(float(T[s2].mean()), 2) if s2.sum() else None
    r = np.hypot(x - 8.25, y - 26.88)
    refs = {}
    for rp in (1.0, 2.0, 2.5, 3.0, 4.0):
        m = (np.abs(z - 0.745) < 0.02) & (r <= rp)
        refs[str(rp)] = round(float(T[m].mean()), 2) if m.sum() else None
    rec["U21_board_ref_C"] = refs
    rec["U21_board_ref_r2p5mm_C"] = refs["2.5"]
    return rec


def conservation(casedir, mdir, hs):
    N = {}
    for ln in open(os.path.join(mdir, "mesh.nodes")):
        p = ln.split()
        if len(p) < 5:
            continue
        try:
            N[int(p[0])] = (round(float(p[2]), 9), round(float(p[3]), 9), round(float(p[4]), 9))
        except ValueError:
            pass
    grp = {}
    for ln in open(os.path.join(mdir, "mesh.boundary")):
        p = ln.split()
        if len(p) < 6:
            continue
        try:
            grp.setdefault(int(p[1]), []).append([int(t) for t in p[5:]])
        except ValueError:
            pass
    P, T = read_vtu(casedir)
    cmap = {(round(P[j, 0], 9), round(P[j, 1], 9), round(P[j, 2], 9)): j for j in range(len(P))}
    Tamb = hs["Tamb_C"] + 273.15; h = hs["h_W_per_m2K"]; eps = hs["emissivity"]
    hfin = hs["G_fin_per_terminal_W_per_K"] / (3.81e-3 * 3.81e-3) if hs.get("fin_on") else 0.0
    SB = 5.67e-8
    res = {}
    for tag, isterm in [(1, False), (2, True), (3, True)]:
        if tag not in grp:
            continue
        A = Qc = Qr = 0.0
        for tri in grp[tag]:
            try:
                pts = [np.array(N[i]) for i in tri]; Tt = [T[cmap[N[i]]] for i in tri]
            except KeyError:
                continue
            a, b, c = pts
            ar = 0.5 * np.linalg.norm(np.cross(b - a, c - a)); A += ar
            Tf = np.mean(Tt); hh = h + (hfin if isterm else 0)
            Qc += hh * (Tf - Tamb) * ar
            if eps is not None:
                Qr += eps * SB * (Tf ** 4 - Tamb ** 4) * ar
        res[tag] = {"area_mm2": round(A * 1e6, 1), "Qconv_W": round(Qc, 3), "Qrad_W": round(Qr, 3)}
    Qt = sum(v["Qconv_W"] + v["Qrad_W"] for v in res.values())
    return {"Ptot_W": hs["total_W"], "Qout_W": round(Qt, 3),
            "imbalance_pct": round(100 * (Qt - hs["total_W"]) / hs["total_W"], 2), "detail": res}


def run_case(geo, level, Iout, Tamb=25.0, hmax=4, tag=None):
    if tag is None:
        tag = f"{geo}_L{int(round(Iout * 10)):02d}"
    cd = os.path.join(ROOT, "cases", tag)
    os.makedirs(cd, exist_ok=True)
    bc = setup_elmer(geo, level)
    mdir = os.path.join(f"elmer_{geo}_{level}", "mesh")
    h = H0[Iout]; hist = []
    for it in range(hmax):
        hs = build_case(geo, level, Iout, h, Tamb, bc, cd)
        rc, log = solve(cd)
        if rc != 0:
            print("SOLVE FAIL", tag, log[-1200:], flush=True); return None
        Ts = surf_mean_T(cd, mdir, bc["outer"])
        dT = Ts - (Tamb + 273.15)
        hnew = h_natural(dT, Tamb)
        hist.append({"it": it, "h": round(h, 3), "Tsurf_avg_C": round(Ts - 273.15, 2),
                     "dT": round(dT, 2), "h_next": round(hnew, 3)})
        print(f"  {tag} h-it{it} h={h:.2f} Tsurf={Ts-273.15:.1f} h_new={hnew:.2f}", flush=True)
        if abs(hnew - h) / h < 0.02:
            h = hnew; break
        h = 0.5 * h + 0.5 * hnew
    hs = build_case(geo, level, Iout, h, Tamb, bc, cd)
    rc, log = solve(cd)
    if rc != 0:
        print("FINAL SOLVE FAIL", tag, log[-1200:], flush=True); return None
    rec = metrics(cd)
    rec["U21_P_W"] = hs["powers_W"]["U21"]
    rec["U21_thetaJB_C_W"] = hs["comp_Rtheta_C_per_W"]["U21"]
    rec["U21_Tj_best_C"] = round(rec["U21_board_ref_C"]["2.5"] + rec["U21_P_W"] * 19.3, 2)
    rec["U21_Tj_conservative_C"] = round(rec["U21_iface_C"] + rec["U21_P_W"] * 19.3, 2)
    rec.update({"geo": geo, "mesh_level": level, "Iout_A": Iout, "Tamb_C": Tamb,
                "h_final": h, "h_hist": hist, "total_W": hs["total_W"],
                "powers_W": hs["powers_W"], "efficiency_est": hs["efficiency_est"],
                "board_k_xy_bg_W_mK": hs["board_k_xy_bg_W_mK"],
                "board_k_xy_cut_W_mK": hs["board_k_xy_cut_W_mK"],
                "board_k_z_W_mK": hs["board_k_z_W_mK"]})
    rec["Tsurf_avg_C"] = round(surf_mean_T(cd, mdir, bc["outer"]) - 273.15, 2)
    rec["conservation"] = conservation(cd, mdir, hs)
    rec["gate_Tj_lt_150"] = bool(rec["U21_Tj_conservative_C"] < 150.0)
    json.dump(rec, open(os.path.join(cd, "result.json"), "w"), indent=1)
    print(tag, "KXYcut=%.2f" % hs["board_k_xy_cut_W_mK"],
          {k: rec[k] for k in ("board_mean_C", "U21_iface_C", "U21_Tj_best_C",
                                "U21_Tj_conservative_C", "D1_body_max_C", "L2_body_max_C",
                                "L1_body_max_C", "Tmax_body_C")}, flush=True)
    return rec


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    numbers = {}
    npath = os.path.join(ROOT, "numbers_thermal_3x3_real.json")
    if os.path.exists(npath):
        numbers = json.load(open(npath))
    for geo in GEOS:
        build_mesh(geo, "main")
        build_mesh(geo, "coarse")
    if stage == "mesh":
        print("MESH DONE", flush=True); return
    conv = {}
    for geo in GEOS:
        cm = run_case(geo, "coarse", 3.0)
        mm = run_case(geo, "main", 3.0)
        for lvl, rec in [("coarse", cm), ("main", mm)]:
            if rec:
                numbers[f"CONV_{geo}_{lvl}"] = rec
        json.dump(numbers, open(npath, "w"), indent=1)
        if cm and mm:
            conv[geo] = {"rel_change_pct": {
                "board_mean_C": round(abs(mm["board_mean_C"] - cm["board_mean_C"]) / abs(mm["board_mean_C"]) * 100, 3),
                "U21_board_ref_r2p5mm_C": round(abs(mm["U21_board_ref_r2p5mm_C"] - cm["U21_board_ref_r2p5mm_C"]) / abs(mm["U21_board_ref_r2p5mm_C"]) * 100, 3),
                "U21_Tj_best_C": round(abs(mm["U21_Tj_best_C"] - cm["U21_Tj_best_C"]) / abs(mm["U21_Tj_best_C"]) * 100, 3)}}
    json.dump(conv, open(os.path.join(ROOT, "convergence_3x3_real.json"), "w"), indent=1)
    for geo in GEOS:
        for Iout in LOADS:
            rec = run_case(geo, "main", Iout)
            if rec:
                numbers[f"M_{geo}_L{int(round(Iout * 10)):02d}"] = rec
            json.dump(numbers, open(npath, "w"), indent=1)
    json.dump(numbers, open(npath, "w"), indent=1)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
