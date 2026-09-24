#!/usr/bin/env python3
"""poly_from_epru.py — re-parse the SOURCE cutout polygons directly from
epro/pourSim.epru (EasyEDA Pro project), for the thermal_3x3_real rerun.

Read-only. Confirms the polygon vertices against geom_verify/copper_areas.json
and writes cutout_polygons.json (with an explicit coordinate conversion).

Coordinate conventions (documented for the mesh builder):
  epru REGION vertices are in mil, y is DOWN-positive (y_epru).
  Thermal model frame (same as thermal_3x3 build_mesh): x in [0,26.50] mm,
  y in [0,50.50] mm, y UP-positive, board centred so that a component at
  epru (cx_mil, cy_mil) sits at mesh (cx_mil*0.0254, 50.50 + cy_mil*0.0254).
  Verified against the L2/L1 keep-out rectangles used by thermal_3x3
  (L2 bbox 6.600x7.000 mm at centre (6.223,33.736) mm).
  => CONVERSION:  x_mm = x_mil * 0.0254 ;  y_mm = 50.50 + y_mil * 0.0254
"""
import json, math, os

EPRU = "/mnt/raid10/sim-work/tps62933/epro/pourSim.epru"
OUT = "/mnt/raid10/sim-work/tps62933/thermal_3x3_real"
MIL2MM = 0.0254
BOARD_W_MM, BOARD_H_MM = 26.50, 50.50
os.makedirs(OUT, exist_ok=True)


def load_docs(f):
    lines = open(f, encoding="utf-8", errors="replace").read().split("\n")
    docs = []; cur = None
    for l in lines:
        if not l.strip() or "||" not in l:
            continue
        h, p = l.split("||", 1)
        if p.endswith("|"):
            p = p[:-1]
        try:
            hj = json.loads(h)
        except Exception:
            continue
        t = hj.get("type")
        if t == "DOCHEAD":
            try:
                pj = json.loads(p)
            except Exception:
                pj = {}
            cur = {"type": pj.get("docType"), "uuid": pj.get("uuid"), "items": []}
            docs.append(cur)
        elif cur is not None:
            try:
                pj = json.loads(p)
            except Exception:
                pj = {}
            cur["items"].append((t, pj, hj))
    return docs


def r_rect(x, y, w, h, rot):
    corners = [(0, 0), (w, 0), (w, -h), (0, -h)]
    th = math.radians(rot); ct, st = math.cos(th), math.sin(th)
    return [(x + dx * ct - dy * st, y + dx * st + dy * ct) for dx, dy in corners]


def parse_path(path):
    if not isinstance(path, list):
        return []
    if len(path) == 1 and isinstance(path[0], list):
        path = path[0]
    if path and path[0] == "R":
        x, y, w, h = path[1], path[2], path[3], path[4]
        rot = path[5] if len(path) > 5 else 0
        return r_rect(x, y, w, h, rot)
    pts = []; i = 0; n = len(path)
    while i < n:
        tok = path[i]
        if isinstance(tok, str):
            if tok == "ARC":
                i += 1
                if i < n and isinstance(path[i], (int, float)):
                    i += 1
                continue
            i += 1; continue
        else:
            if i + 1 < n and isinstance(path[i + 1], (int, float)):
                pts.append((float(path[i]), float(path[i + 1]))); i += 2
            else:
                i += 1
    return pts


def poly_area(pts):
    a = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]; x2, y2 = pts[(i + 1) % len(pts)]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def bbox(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return [min(xs), min(ys), max(xs), max(ys)]


def to_mesh(pts):
    """epru mil (y down) -> thermal mesh mm (y up, 0..50.5)."""
    return [(round(x * MIL2MM, 5), round(BOARD_H_MM + y * MIL2MM, 5)) for (x, y) in pts]


def main():
    docs = load_docs(EPRU)
    pcbs = [d for d in docs if d["type"] == "PCB"]
    titles = {}
    for i, d in enumerate(pcbs):
        titles[i] = [p for t, p, _ in d["items"] if t == "META"][0].get("title")

    def region_ids(doc):
        out = {}
        for t, p, hj in doc["items"]:
            if t == "REGION" and p.get("regionType") == "PROHIBIT":
                out[hj.get("id")] = p
        return out

    ids = [region_ids(d) for d in pcbs]
    new1 = set(ids[1].keys()) - set(ids[0].keys())
    new2 = set(ids[2].keys()) - set(ids[0].keys())

    def recs(i, newids):
        out = []
        for t, p, hj in pcbs[i]["items"]:
            if t != "REGION" or p.get("regionType") != "PROHIBIT":
                continue
            if hj.get("id") not in newids:
                continue
            pts = parse_path(p.get("path"))
            out.append({"id": hj.get("id"), "layer": p.get("layerId"),
                        "vertices_mil": [[round(x, 3), round(y, 3)] for x, y in pts],
                        "bbox_mil": [round(v, 3) for v in bbox(pts)],
                        "poly_area_mm2": round(poly_area(pts) * MIL2MM ** 2, 3)})
        return out

    top = recs(1, new1)
    dual = recs(2, new2)

    # cross-check vs geom_verify/copper_areas.json
    ref_path = "/mnt/raid10/sim-work/tps62933/geom_verify/copper_areas.json"
    checks = []
    if os.path.exists(ref_path):
        ref = json.load(open(ref_path))
        for key, got in [("cutouts_topcut", top), ("cutouts_dualcut", dual)]:
            refmap = {r["id"]: r for r in ref[key]}
            for g in got:
                r = refmap.get(g["id"])
                if r is None:
                    checks.append({"id": g["id"], "status": "MISSING_IN_REF"}); continue
                dv = max(abs(a - b) for va, vb in zip(g["vertices_mil"], r["vertices_mil"])
                         for a, b in zip(va, vb))
                da = abs(g["poly_area_mm2"] - r["poly_area_mm2"])
                checks.append({"id": g["id"], "n_v": len(g["vertices_mil"]),
                               "max_vertex_diff_mil": round(dv, 6),
                               "area_diff_mm2": round(da, 6),
                               "status": "MATCH" if (dv < 1e-6 and da < 1e-6) else "DIFF"})

    # map regions to labels by bbox centre: L2 region has centre ~ (245,-660),
    #   L1 region centre ~ (645,-1250)  (mil, epru frame)
    def label(r):
        c = ((r["bbox_mil"][0] + r["bbox_mil"][2]) / 2.0,
             (r["bbox_mil"][1] + r["bbox_mil"][3]) / 2.0)
        return "L2" if c[1] > -950 else "L1"

    poly = {"topcut": {}, "dualcut": {}}
    for key, lst in [("topcut", top), ("dualcut", dual)]:
        for r in lst:
            poly[key][label(r)] = {"id": r["id"], "layer": r["layer"],
                                   "vertices_mil": r["vertices_mil"],
                                   "vertices_mesh_mm": to_mesh(r["vertices_mil"]),
                                   "bbox_mil": r["bbox_mil"],
                                   "poly_area_mm2": r["poly_area_mm2"]}

    out = {
        "note": "cutout polygons re-parsed from epro/pourSim.epru; "
                "conversion x_mm=x_mil*0.0254, y_mm=50.50+y_mil*0.0254 (y flipped to up-positive)",
        "source": EPRU,
        "pcb_titles": {str(i): titles[i] for i in titles},
        "new_region_ids_topcut": sorted(new1),
        "new_region_ids_dualcut": sorted(new2),
        "cross_check_vs_geom_verify": checks,
        "polygons": poly,
    }
    json.dump(out, open(os.path.join(OUT, "cutout_polygons.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
