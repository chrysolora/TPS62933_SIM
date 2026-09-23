from collections import Counter
lines = open("/mnt/raid10/sim-work/tps62933/fea/case1/mesh.msh", errors="replace").read().splitlines()
# locate Elements
idx = None
for i, l in enumerate(lines):
    if l.strip() == "$Elements":
        idx = i
        break
i = idx + 1
n = int(lines[i].split()[0])
i += 1
c = Counter()
for k in range(n):
    parts = lines[i + k].split()
    et = int(parts[1]); ntags = int(parts[2])
    tags = parts[3:3 + ntags]
    if et == 4:
        c[tags[-1]] += 1
print("tets per physical volume tag:", dict(sorted(c.items(), key=lambda x: int(x[0]))))
