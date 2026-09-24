import vtk, numpy as np, sys, os, json
d = sys.argv[1]; out = sys.argv[2]
BINS = [(0, 3), (3, 5), (5, 7), (7, 10), (10, 15), (15, 25), (25, 40)]
RHO_IN = 3.25  # inductor half-size mm

def load(p):
    r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(p); r.Update()
    g = r.GetOutput(); n = g.GetNumberOfPoints()
    pts = np.array([g.GetPoint(i) for i in range(n)])
    J = np.array(g.GetPointData().GetArray("Jmag"))
    rho = np.sqrt(pts[:, 0]**2 + pts[:, 1]**2) * 1e3
    return pts, J, rho

def stats(p):
    if not os.path.exists(p): return None
    pts, J, rho = load(p)
    if len(J) == 0: return None
    s = dict(n=int(len(J)), Jmax=float(J.max()), r_at_max=float(rho[np.argmax(J)]),
             frac_gt1e7=float((J > 1e7).mean()), frac_gt5e6=float((J > 5e6).mean()),
             mean_in=float(J[rho < RHO_IN].mean()) if (rho < RHO_IN).sum() else 0.0,
             n_in=int((rho < RHO_IN).sum()))
    s["bins"] = []
    for lo, hi in BINS:
        m = (rho >= lo) & (rho < hi)
        if m.sum():
            s["bins"].append([lo, hi, int(m.sum()), float(J[m].mean()), float(J[m].max())])
        else:
            s["bins"].append([lo, hi, 0, 0.0, 0.0])
    return s

def comb(p1, p2):
    pd = []
    for p in (p1, p2):
        if os.path.exists(p):
            _, J, rho = load(p)
            pd.append((rho, J))
    rho = np.concatenate([a for a, _ in pd]) if pd else np.array([])
    J = np.concatenate([b for _, b in pd]) if pd else np.array([])
    s = dict(n=int(len(J)), Jmax=float(J.max()) if len(J) else 0.0,
             mean_in=float(J[rho < RHO_IN].mean()) if len(J) and (rho < RHO_IN).sum() else 0.0)
    s["bins"] = []
    for lo, hi in BINS:
        m = (rho >= lo) & (rho < hi)
        s["bins"].append([lo, hi, int(m.sum()), float(J[m].mean()) if m.sum() else 0.0,
                          float(J[m].max()) if m.sum() else 0.0])
    return s

res = dict(top=stats(d + "/mesh/copper_top.vtu"),
           bot=stats(d + "/mesh/copper_bot.vtu"),
           comb=comb(d + "/mesh/copper_top.vtu", d + "/mesh/copper_bot.vtu"))
json.dump(res, open(out, "w"), indent=1)
print(json.dumps(res, indent=1))
