lines = open("/mnt/raid10/sim-work/tps62933/fea/case1/mesh.msh", errors="replace").read().splitlines()
# nodes
ni = [i for i, l in enumerate(lines) if l.strip() == "$Nodes"][0]
n = int(lines[ni + 1].split()[0])
nodes = {}
for k in range(n):
    p = lines[ni + 2 + k].split()
    nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
ei = [i for i, l in enumerate(lines) if l.strip() == "$Elements"][0]
ne = int(lines[ei + 1].split()[0])
from collections import defaultdict
agg = defaultdict(lambda: [0, 0.0, 0.0, 0.0])
for k in range(ne):
    p = lines[ei + 2 + k].split()
    et = int(p[1]); nt = int(p[2]); tags = p[3:3 + nt]; nids = [int(x) for x in p[3 + nt:]]
    if et != 4:
        continue
    tag = tags[-1]
    cx = sum(nodes[i][0] for i in nids) / 4
    cy = sum(nodes[i][1] for i in nids) / 4
    cz = sum(nodes[i][2] for i in nids) / 4
    a = agg[tag]
    a[0] += 1; a[1] += cx; a[2] += cy; a[3] += cz
for tag, a in sorted(agg.items(), key=lambda x: int(x[0])):
    print("tag %s : n=%d centroid=(%.2f, %.2f, %.2f) mm" % (tag, a[0], 1000*a[1]/a[0], 1000*a[2]/a[0], 1000*a[3]/a[0]))
