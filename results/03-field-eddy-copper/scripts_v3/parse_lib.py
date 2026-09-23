import json

def load_pcbs(f):
    lines=open(f,encoding="utf-8",errors="replace").read().split("\n")
    docs=[];cur=None
    for l in lines:
        if not l.strip() or "||" not in l: continue
        h,p=l.split("||",1)
        if p.endswith("|"): p=p[:-1]
        try: hj=json.loads(h)
        except: continue
        t=hj.get("type")
        if t=="DOCHEAD":
            try: pj=json.loads(p)
            except: pj={}
            cur={"type":pj.get("docType"),"uuid":pj.get("uuid"),"items":[]};docs.append(cur)
        elif cur is not None:
            try: pj=json.loads(p)
            except: pj={}
            cur["items"].append((t,pj))
    return [d for d in docs if d["type"]=="PCB"]

def r_rect(x,y,w,h,rot):
    """R = corner (x,y), extends +x and -y; rotate about corner by rot deg (screen, y down)."""
    import math
    corners=[(0,0),(w,0),(w,-h),(0,-h)]
    th=math.radians(rot)
    ct,st=math.cos(th),math.sin(th)
    pts=[]
    for dx,dy in corners:
        rx=dx*ct - dy*st
        ry=dx*st + dy*ct
        pts.append((x+rx, y+ry))
    return pts

def parse_path(path):
    """Return list of (x,y) points; handles R and L/ARC."""
    if not isinstance(path,list): return []
    # unwrap single-element nested
    if len(path)==1 and isinstance(path[0],list): path=path[0]
    if path and path[0]=="R":
        x,y,w,h=path[1],path[2],path[3],path[4]
        rot=path[5] if len(path)>5 else 0
        return r_rect(x,y,w,h,rot)
    pts=[];i=0;n=len(path)
    while i<n:
        tok=path[i]
        if isinstance(tok,str):
            if tok=="ARC":
                i+=1
                if i<n and isinstance(path[i],(int,float)): i+=1  # angle
                continue
            i+=1;continue
        else:
            if i+1<n and isinstance(path[i+1],(int,float)):
                pts.append((float(path[i]),float(path[i+1])));i+=2
            else: i+=1
    return pts
