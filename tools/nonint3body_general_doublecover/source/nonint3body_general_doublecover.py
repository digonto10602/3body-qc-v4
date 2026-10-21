#!/usr/bin/env python3
import argparse, math, os, sys, itertools, time
from dataclasses import dataclass
from typing import List, Tuple, Dict
import numpy as np

# ----------------------------- config helpers -----------------------------
@dataclass(frozen=True)
class Vec3:
    x:int; y:int; z:int
    def __iter__(self): return iter((self.x,self.y,self.z))
    def __add__(self,o): return Vec3(self.x+o.x,self.y+o.y,self.z+o.z)
    def __sub__(self,o): return Vec3(self.x-o.x,self.y-o.y,self.z-o.z)
    def norm2(self): return self.x*self.x+self.y*self.y+self.z*self.z
    def tag(self): return f"{self.x}{self.y}{self.z}"
    def s(self): return f"[{self.x},{self.y},{self.z}]"

@dataclass
class Particle:
    name:str; mass:float; J:float; parity:int; C:str

@dataclass
class Config:
    Emin:float=0.25; Emax:float=0.38; L:float=20.0; xi:float=3.444
    nmax:int=5; group_tol:float=1e-10; irrep_tol:float=1e-8
    P_list:List[Vec3]=None
    particles:List[Particle]=None
    outdir:str="output"; tag:str="nonint3body_general_dc_debug"
    double_cover:str="auto"   # auto,y,n
    print_empty_irreps:bool=False

def trim(s): return s.strip()
def read_kv(path):
    kv={}
    with open(path) as f:
        for line in f:
            line=line.split('#',1)[0].strip()
            if not line or '=' not in line: continue
            k,v=line.split('=',1); kv[k.strip()]=v.strip()
    return kv
def gs(kv,k,d): return kv.get(k,d)
def gd(kv,k,d): return float(kv.get(k,d))
def gi(kv,k,d): return int(kv.get(k,d))
def gb(kv,k,d):
    v=kv.get(k, 'y' if d else 'n').strip().lower()
    return v in ('y','yes','true','1','on')
def parse_parity(s): return -1 if '-' in s else +1
def parse_P(s):
    s=s.strip().replace('[','').replace(']','').replace(',',' ').replace('_',' ')
    parts=s.split()
    if len(parts)==3: return Vec3(int(parts[0]),int(parts[1]),int(parts[2]))
    s=''.join(ch for ch in s if not ch.isspace())
    if len(s)==3 and all(ch.isdigit() or ch=='-' for ch in s): return Vec3(int(s[0]),int(s[1]),int(s[2]))
    raise ValueError(f"cannot parse P token '{s}'")
def read_config(path):
    kv=read_kv(path)
    c=Config()
    c.Emin=gd(kv,'Ecm_min',c.Emin); c.Emax=gd(kv,'Ecm_max',c.Emax)
    c.L=gd(kv,'L',c.L); c.xi=gd(kv,'xi',c.xi); c.nmax=gi(kv,'nmax',c.nmax)
    c.group_tol=gd(kv,'energy_group_tol',c.group_tol); c.irrep_tol=gd(kv,'irrep_tol',c.irrep_tol)
    c.outdir=gs(kv,'output_dir',c.outdir); c.tag=gs(kv,'output_tag',c.tag)
    c.double_cover=gs(kv,'double_cover',c.double_cover).lower()
    c.print_empty_irreps=gb(kv,'print_empty_irreps',False)
    c.P_list=[parse_P(x) for x in gs(kv,'P_list','000 001 110 111 200').split()]
    parts=[]
    for i in range(1,4):
        pre=f'particle{i}_'
        parts.append(Particle(gs(kv,pre+'name',f'p{i}'), gd(kv,pre+'mass',0.1), gd(kv,pre+'J',0.0), parse_parity(gs(kv,pre+'P','+')), gs(kv,pre+'C','none')))
    c.particles=parts
    return c

def is_half_integer(J): return abs(2*J-round(2*J))<1e-9 and int(round(2*J))%2==1

# ----------------------------- integer matrices -----------------------------
def mat_tuple(M): return tuple(int(x) for row in M for x in row)
def mat_mul(A,B):
    return tuple(sum(A[3*i+k]*B[3*k+j] for k in range(3)) for i in range(3) for j in range(3))
def mat_apply(M,v:Vec3):
    a=list(M); x=[v.x,v.y,v.z]
    r=[sum(a[3*i+j]*x[j] for j in range(3)) for i in range(3)]
    return Vec3(*r)
def mat_det(M):
    a=M
    return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6])
def mat_trace(M): return M[0]+M[4]+M[8]
def mat_neg(M): return tuple(-x for x in M)
def signed_perm_mats():
    mats=[]
    for p in itertools.permutations(range(3)):
        for signs in itertools.product([-1,1], repeat=3):
            A=[0]*9
            for i in range(3): A[3*i+p[i]]=signs[i]
            mats.append(tuple(A))
    return mats

def little_group_mats(P):
    return [M for M in signed_perm_mats() if mat_apply(M,P)==P]

# ----------------------------- quaternions/SU2 lifts -----------------------------
def q_norm(q):
    n=math.sqrt(sum(x*x for x in q)); return tuple(x/n for x in q)
def q_neg(q): return tuple(-x for x in q)
def q_mul(a,b):
    w1,x1,y1,z1=a; w2,x2,y2,z2=b
    return q_norm((w1*w2-x1*x2-y1*y2-z1*z2,
                   w1*x2+x1*w2+y1*z2-z1*y2,
                   w1*y2-x1*z2+y1*w2+z1*x2,
                   w1*z2+x1*y2-y1*x2+z1*w2))
def q_key(q, inv): return (tuple(round(x,10) for x in q_norm(q)), int(inv))

def quat_to_mat(q):
    w,x,y,z=q_norm(q)
    R=np.array([[1-2*y*y-2*z*z, 2*x*y-2*z*w, 2*x*z+2*y*w],
                [2*x*y+2*z*w, 1-2*x*x-2*z*z, 2*y*z-2*x*w],
                [2*x*z-2*y*w, 2*y*z+2*x*w, 1-2*x*x-2*y*y]])
    return tuple(int(round(v)) for v in R.reshape(9))

def mat_to_quat(R):
    # R is a proper signed permutation rotation. Return one SU2 lift.
    m=np.array(R,dtype=float).reshape(3,3)
    tr=np.trace(m)
    if tr>0:
        S=math.sqrt(tr+1.0)*2.0
        w=0.25*S; x=(m[2,1]-m[1,2])/S; y=(m[0,2]-m[2,0])/S; z=(m[1,0]-m[0,1])/S
    elif (m[0,0]>m[1,1]) and (m[0,0]>m[2,2]):
        S=math.sqrt(1.0+m[0,0]-m[1,1]-m[2,2])*2.0
        w=(m[2,1]-m[1,2])/S; x=0.25*S; y=(m[0,1]+m[1,0])/S; z=(m[0,2]+m[2,0])/S
    elif m[1,1]>m[2,2]:
        S=math.sqrt(1.0+m[1,1]-m[0,0]-m[2,2])*2.0
        w=(m[0,2]-m[2,0])/S; x=(m[0,1]+m[1,0])/S; y=0.25*S; z=(m[1,2]+m[2,1])/S
    else:
        S=math.sqrt(1.0+m[2,2]-m[0,0]-m[1,1])*2.0
        w=(m[1,0]-m[0,1])/S; x=(m[0,2]+m[2,0])/S; y=(m[1,2]+m[2,1])/S; z=0.25*S
    return q_norm((w,x,y,z))

@dataclass(frozen=True)
class DElem:
    q:Tuple[float,float,float,float]
    inv:int
    G:Tuple[int,...]

def build_double_group(P):
    elems=[]; seen={}
    for G in little_group_mats(P):
        d=mat_det(G)
        inv=0 if d==1 else 1
        R=G if d==1 else mat_neg(G)  # G = inversion^inv * R
        q=mat_to_quat(R)
        for qq in [q, q_neg(q)]:
            e=DElem(q_norm(qq), inv, G)
            seen[q_key(e.q,e.inv)]=e
    elems=list(seen.values())
    # close under multiplication to guard against rounding/lift choices
    changed=True
    while changed:
        changed=False
        vals=list(seen.values())
        for a in vals:
            for b in vals:
                q=q_mul(a.q,b.q); inv=a.inv^b.inv
                R=quat_to_mat(q); G=R if inv==0 else mat_neg(R)
                if mat_apply(G,P)!=P: continue
                k=q_key(q,inv)
                if k not in seen:
                    seen[k]=DElem(q_norm(q),inv,G); changed=True
    elems=list(seen.values())
    key_to_idx={q_key(e.q,e.inv):i for i,e in enumerate(elems)}
    mult=np.empty((len(elems),len(elems)),dtype=int)
    for i,a in enumerate(elems):
        for j,b in enumerate(elems):
            q=q_mul(a.q,b.q); inv=a.inv^b.inv
            mult[i,j]=key_to_idx[q_key(q,inv)]
    # identity and inverses
    ididx=key_to_idx[q_key((1,0,0,0),0)]
    invs=[]
    for i in range(len(elems)):
        found=None
        for j in range(len(elems)):
            if mult[i,j]==ididx and mult[j,i]==ididx: found=j; break
        invs.append(found)
    return elems,mult,ididx,invs

def conjugacy_classes(elems,mult,invs):
    n=len(elems); unused=set(range(n)); classes=[]; elem_to_class={}
    while unused:
        a=next(iter(unused)); cls=set()
        for g in range(n): cls.add(mult[mult[g,a],invs[g]])
        cls=sorted(cls); ci=len(classes)
        for x in cls: elem_to_class[x]=ci
        classes.append(cls); unused-=set(cls)
    return classes, elem_to_class

# ----------------------------- automatic irreducible characters -----------------------------
def class_algebra_irreps(elems,mult,classes,elem_to_class,ididx):
    nc=len(classes); sizes=np.array([len(c) for c in classes],dtype=float); order=len(elems)
    # class multiplication matrices A_i: C_i*C_j=sum_k N_ij^k C_k, columns are j, rows are k.
    A=[]
    for i,ci in enumerate(classes):
        Ai=np.zeros((nc,nc),dtype=complex)
        for j,cj in enumerate(classes):
            counts=np.zeros(nc,dtype=float)
            for a in ci:
                for b in cj:
                    counts[elem_to_class[mult[a,b]]]+=1
            for k in range(nc): Ai[k,j]=counts[k]/sizes[k]
        A.append(Ai)
    rng=np.random.default_rng(12345)
    for attempt in range(20):
        coeff=rng.normal(size=nc)
        B=sum(coeff[i]*A[i] for i in range(nc))
        vals,vecs=np.linalg.eig(B)
        # if degenerate, retry
        if min(abs(vals[i]-vals[j]) for i in range(nc) for j in range(i)) if nc>1 else 1:
            pass
        irreps=[]; seen=[]
        ok=True
        for r in range(nc):
            v=vecs[:,r]
            # normalize arbitrary; compute lambda_i by Rayleigh quotient
            den=np.vdot(v,v)
            lambdas=[]
            resid=0.0
            for Ai in A:
                Av=Ai@v; lam=np.vdot(v,Av)/den; lambdas.append(lam)
                resid=max(resid, np.linalg.norm(Av-lam*v)/(np.linalg.norm(v)+1e-30))
            if resid>1e-6:
                ok=False; break
            denom=sum((abs(lambdas[i])**2)/sizes[i] for i in range(nc))
            d=math.sqrt(order/denom.real)
            chars=np.array([lambdas[i]*d/sizes[i] for i in range(nc)],dtype=complex)
            # fix roundoff
            chars=np.where(abs(chars.real)<1e-9, 0, chars.real) + 1j*np.where(abs(chars.imag)<1e-9,0,chars.imag)
            # ensure identity dim positive
            eid_cls=elem_to_class[ididx]
            dim=chars[eid_cls].real
            if dim<0: chars=-chars; dim=-dim
            # duplicate check
            if any(np.linalg.norm(chars-x)<1e-6 for x in seen): continue
            seen.append(chars)
            irreps.append({'dim':int(round(dim)), 'chars':chars})
        if ok and len(irreps)==nc and abs(sum(ir['dim']**2 for ir in irreps)-order)<1e-6:
            return irreps
    raise RuntimeError('failed to extract irreducible characters from class algebra')

# ----------------------------- labeling -----------------------------
def old_single_group_name(P, group_order):
    # In double-cover mode group_order = 2 * |single-cover little group|.
    # The common momentum classes are:
    # P=000: 2O_h order 96; [00n]: Dic4/C4v^D order 16;
    # [nn0]: Dic2/C2v^D order 8; [nnn]: Dic3/C3v^D order 12.
    if P.norm2()==0: return 'Oh'
    if group_order==16: return 'C4v'
    if group_order==12: return 'C3v'
    if group_order==8: return 'C2v'
    if group_order==4: return 'C2/Cs'
    if group_order==2: return 'C1'
    return 'C1'

def class_Oh(G):
    if mat_det(G)==1:
        t=mat_trace(G)
        if t==3: return 'E'
        if t==0: return 'C3'
        if t==1: return 'C4'
        if t==-1:
            diag=all((i==j or G[3*i+j]==0) for i in range(3) for j in range(3))
            return 'C2coord' if diag else 'C2edge'
    else:
        return 'i*'+class_Oh(mat_neg(G))
    return 'unknown'
def class_C4v(G,P):
    if mat_det(G)==1:
        if mat_trace(G)==3: return 'E'
        if mat_trace(G)==1: return 'C4'
        return 'C2'
    ax=0 if P.x!=0 and P.y==0 and P.z==0 else (1 if P.y!=0 and P.x==0 and P.z==0 else 2)
    perp=[i for i in range(3) if i!=ax]
    off=sum(1 for i in perp for j in perp if i!=j and G[3*i+j]!=0)
    return 'sigmad' if off else 'sigmav'
def class_C3v(G,P):
    if mat_det(G)==1: return 'E' if mat_trace(G)==3 else 'C3'
    return 'sigma'
def class_C2v(G,P):
    if mat_det(G)==1: return 'E' if mat_trace(G)==3 else 'C2'
    zaxis=0 if P.x==0 else (1 if P.y==0 else (2 if P.z==0 else -1))
    return 'sigma_zeroaxis' if zaxis>=0 and G[3*zaxis+zaxis]==-1 else 'sigma_swapaxis'
def class_C2(G,P): return 'E' if mat_trace(G)==3 else 'g'
def single_class(group,G,P):
    if group=='Oh': return class_Oh(G)
    if group=='C4v': return class_C4v(G,P)
    if group=='C3v': return class_C3v(G,P)
    if group=='C2v': return class_C2v(G,P)
    if group=='C2/Cs': return class_C2(G,P)
    return 'E'

single_tables={
 'Oh': {
  'classes':['E','C3','C2coord','C4','C2edge','i*E','i*C3','i*C2coord','i*C4','i*C2edge'],
  'irreps':{
   **{name+'p':vals+vals for name,vals in {'A1':[1,1,1,1,1], 'A2':[1,1,1,-1,-1], 'E':[2,-1,2,0,0], 'T1':[3,0,-1,1,-1], 'T2':[3,0,-1,-1,1]}.items()},
   **{name+'m':vals+[-x for x in vals] for name,vals in {'A1':[1,1,1,1,1], 'A2':[1,1,1,-1,-1], 'E':[2,-1,2,0,0], 'T1':[3,0,-1,1,-1], 'T2':[3,0,-1,-1,1]}.items()}
  }
 },
 'C4v': {'classes':['E','C4','C2','sigmav','sigmad'], 'irreps':{'A1':[1,1,1,1,1], 'A2':[1,1,1,-1,-1], 'B1':[1,-1,1,1,-1], 'B2':[1,-1,1,-1,1], 'E':[2,0,-2,0,0]}},
 'C3v': {'classes':['E','C3','sigma'], 'irreps':{'A1':[1,1,1], 'A2':[1,1,-1], 'E':[2,-1,0]}},
 'C2v': {'classes':['E','C2','sigma_zeroaxis','sigma_swapaxis'], 'irreps':{'A1':[1,1,1,1], 'A2':[1,1,-1,-1], 'B1':[1,-1,1,-1], 'B2':[1,-1,-1,1]}},
 'C2/Cs': {'classes':['E','g'], 'irreps':{'A':[1,1], 'B':[1,-1]}},
 'C1': {'classes':['E'], 'irreps':{'A':[1]}}
}

def single_char_for_label(label, group, elem, P):
    tab=single_tables[group]
    cls=single_class(group, elem.G, P)
    ci=tab['classes'].index(cls)
    return tab['irreps'][label][ci]

def find_central_indices(elems, ididx):
    zidx=None; inv_idx=None
    for i,e in enumerate(elems):
        if e.inv==0 and all(abs(e.q[j]-(-1,0,0,0)[j])<1e-8 for j in range(4)): zidx=i
        if e.inv==1 and all(abs(e.q[j]-(1,0,0,0)[j])<1e-8 for j in range(4)): inv_idx=i
    return zidx, inv_idx

def elem_char_from_class(ir, elem_idx, elem_to_class): return ir['chars'][elem_to_class[elem_idx]].real

def label_irreps(irreps, elems, classes, elem_to_class, ididx, P):
    group=old_single_group_name(P, len(elems))
    zidx, inv_idx=find_central_indices(elems,ididx)
    labels=[]; used={}
    # choose C4/C3 representative for spinorial labels
    c4idx=next((i for i,e in enumerate(elems) if e.inv==0 and mat_trace(e.G)==1 and e.q[0]>0), None)
    c3idx=next((i for i,e in enumerate(elems) if e.inv==0 and mat_trace(e.G)==0 and e.q[0]>0), None)
    for ir in irreps:
        dim=ir['dim']
        zratio=elem_char_from_class(ir,zidx,elem_to_class)/dim if zidx is not None else +1
        invratio=elem_char_from_class(ir,inv_idx,elem_to_class)/dim if inv_idx is not None else None
        label=None
        if zratio>0.5:
            # match to single-cover conventional label by character inner product over double group.
            best=None
            for lab in single_tables[group]['irreps']:
                diff=0.0
                for ei,e in enumerate(elems):
                    diff += abs(elem_char_from_class(ir,ei,elem_to_class)-single_char_for_label(lab,group,e,P))**2
                if best is None or diff<best[0]: best=(diff,lab)
            if best and best[0]<1e-6: label=best[1]
            else: label=f"single_dim{dim}"
        else:
            par='' if invratio is None else ('p' if invratio>0 else 'm')
            if group=='Oh':
                if dim==4: base='H'
                elif dim==2:
                    val=elem_char_from_class(ir,c4idx,elem_to_class) if c4idx is not None else 0
                    base='G1' if val>0 else 'G2'
                else: base=f'Spin{dim}'
                label=base+par
            elif group=='C4v':
                val=elem_char_from_class(ir,c4idx,elem_to_class) if c4idx is not None else 0
                base='E1' if val>0 else 'E3'
                label=base
            elif group=='C3v':
                val=elem_char_from_class(ir,c3idx,elem_to_class) if c3idx is not None else 0
                base='G1' if val>0 else 'G2'
                label=base
            elif group=='C2v': label='G'
            elif group=='C2/Cs': label='G'
            else: label=f'Spin{dim}'
        # avoid collisions
        if label in used:
            used[label]+=1; label=f"{label}_{used[label]}"
        else: used[label]=0
        labels.append(label)
    # sort: single labels first conventional, then spinorial; preserve mapping
    out=[]
    for ir,label in zip(irreps,labels):
        out.append({'name':label,'dim':ir['dim'],'chars':ir['chars']})
    return out, group

# ----------------------------- physics enumeration -----------------------------
@dataclass
class Row:
    n:Tuple[Vec3,Vec3,Vec3]; Ecm:float
@dataclass
class Shell:
    E:float; rows:List[Row]

def Ecm_of(cfg,P,ns):
    q=2*math.pi/(cfg.xi*cfg.L)
    Elab=sum(math.sqrt(cfg.particles[i].mass**2 + q*q*ns[i].norm2()) for i in range(3))
    s=Elab*Elab - q*q*P.norm2()
    return math.sqrt(s) if s>0 else float('nan')
def enumerate_shells(cfg,P):
    rows=[]; nmax=cfg.nmax
    rng=range(-nmax,nmax+1)
    for ax,ay,az in itertools.product(rng,rng,rng):
        n1=Vec3(ax,ay,az)
        for bx,by,bz in itertools.product(rng,rng,rng):
            n2=Vec3(bx,by,bz); n3=Vec3(P.x-n1.x-n2.x,P.y-n1.y-n2.y,P.z-n1.z-n2.z)
            if abs(n3.x)>nmax or abs(n3.y)>nmax or abs(n3.z)>nmax: continue
            ns=(n1,n2,n3); E=Ecm_of(cfg,P,ns)
            if cfg.Emin-1e-12 <= E <= cfg.Emax+1e-12: rows.append(Row(ns,E))
    rows.sort(key=lambda r:r.Ecm)
    shells=[]
    for r in rows:
        if not shells or abs(r.Ecm-shells[-1].E)>cfg.group_tol: shells.append(Shell(r.Ecm,[r]))
        else: shells[-1].rows.append(r)
    return shells

def su2_char(J,q):
    q=q_norm(q); w=max(-1.0,min(1.0,q[0])); theta=2*math.acos(w)
    den=math.sin(theta/2)
    if abs(den)<1e-12:
        # q=+1 or q=-1
        return (2*J+1) if w>0 else ((-1)**int(round(2*J)))*(2*J+1)
    return math.sin((2*J+1)*theta/2)/den

def particle_char(p,e): return (p.parity if e.inv else 1.0)*su2_char(p.J,e.q)
def spin_product_char(cfg,e):
    v=1.0
    for p in cfg.particles: v*=particle_char(p,e)
    return v

def shell_character(cfg,P,elems,sh):
    chi=np.zeros(len(elems),dtype=float)
    for ei,e in enumerate(elems):
        fixed=0
        for r in sh.rows:
            if all(mat_apply(e.G,r.n[i])==r.n[i] for i in range(3)): fixed+=1
        chi[ei]=fixed*spin_product_char(cfg,e)
    return chi

# ----------------------------- main -----------------------------
def main():
    ap=argparse.ArgumentParser(description='General three-body non-interacting energies with automatic single/double-cover irrep decomposition.')
    ap.add_argument('config')
    args=ap.parse_args()
    cfg=read_config(args.config); os.makedirs(cfg.outdir,exist_ok=True)
    use_dc = (cfg.double_cover in ('y','yes','true','1','on')) or (cfg.double_cover=='auto' and any(is_half_integer(p.J) for p in cfg.particles))
    print(f"[setup] config={args.config} Ecm=[{cfg.Emin},{cfg.Emax}] L={cfg.L} xi={cfg.xi} nmax={cfg.nmax}", flush=True)
    print(f"[setup] double_cover={cfg.double_cover} -> use_double_cover={'y' if use_dc else 'n'}", flush=True)
    for i,p in enumerate(cfg.particles,1): print(f"[setup] particle{i} name={p.name} mass={p.mass:.17g} J={p.J:.17g} P={'+' if p.parity>0 else '-'} C={p.C}", flush=True)
    allpath=os.path.join(cfg.outdir,cfg.tag+'_all_irreps_nonint3body.dat')
    with open(allpath,'w') as allf:
        allf.write('# General three-body non-interacting energies decomposed by little-group irreps.\n')
        allf.write(f"# double_cover_used {'y' if use_dc else 'n'}\n")
        allf.write('# columns: Ptag group irrep level Ecm group_degeneracy multiplicity irrep_dim raw_energy_multiplicity example_n1 example_n2 example_n3 central_z class_chars\n')
        for ip,P in enumerate(cfg.P_list,1):
            t0=time.time(); elems,mult,ididx,invs=build_double_group(P)
            if not use_dc:
                # keep only + lift to emulate old single cover, but automatic code works either way
                elems=[e for e in elems if e.q[0]>=-1e-10]  # not exact; keep dc mode recommended
                elems,mult,ididx,invs=build_double_group(P) # use dc and label single irreps; spinorial won't appear for integer spins
            classes, elem_to_class=conjugacy_classes(elems,mult,invs)
            raw_irreps=class_algebra_irreps(elems,mult,classes,elem_to_class,ididx)
            irreps, gname=label_irreps(raw_irreps, elems, classes, elem_to_class, ididx, P)
            zidx,_=find_central_indices(elems,ididx)
            shells=enumerate_shells(cfg,P)
            print(f"[P] [{ip}/{len(cfg.P_list)}] P={P.s()} Ptag={P.tag()} group={gname}{'^D' if use_dc else ''} |G|={len(elems)} classes={len(classes)} irreps={len(irreps)} shells={len(shells)}", flush=True)
            print('[P] irreps: '+' '.join(f"{ir['name']}(dim={ir['dim']})" for ir in sorted(irreps,key=lambda x:x['name'])), flush=True)
            files={}
            for ir in irreps:
                path=os.path.join(cfg.outdir,f"{cfg.tag}_P{P.tag()}_{ir['name']}_nonint3body.dat")
                f=open(path,'w'); files[ir['name']]=f
                f.write(f"# P={P.s()} group={gname}{'^D' if use_dc else ''} irrep={ir['name']} |G|={len(elems)}\n")
                f.write('# particles:' + ''.join(f" {p.name}(M={p.mass},J={p.J},P={'+' if p.parity>0 else '-'},C={p.C})" for p in cfg.particles)+'\n')
                f.write('# columns: level Ecm group_degeneracy multiplicity irrep_dim raw_energy_multiplicity example_n1 example_n2 example_n3 central_z class_chars\n')
            levels={ir['name']:0 for ir in irreps}; wrote={ir['name']:0 for ir in irreps}
            for sh in shells:
                chi=shell_character(cfg,P,elems,sh)
                for ir in irreps:
                    irchi=np.array([ir['chars'][elem_to_class[i]].real for i in range(len(elems))])
                    mult_ir=float(np.dot(irchi,chi)/len(elems))
                    # clean numerical roundoff
                    if abs(mult_ir-round(mult_ir))<5e-8: mult_ir=round(mult_ir)
                    gdeg=mult_ir*ir['dim']
                    if gdeg>cfg.irrep_tol or cfg.print_empty_irreps:
                        if gdeg<=cfg.irrep_tol and not cfg.print_empty_irreps: continue
                        ex=sh.rows[0].n; lev=levels[ir['name']]; levels[ir['name']]+=1; wrote[ir['name']]+=1
                        zratio='NA'
                        if zidx is not None: zratio=f"{(ir['chars'][elem_to_class[zidx]].real/ir['dim']):.8g}"
                        # compact class chars by conjugacy class
                        class_ch=[]
                        for ci,cls in enumerate(classes):
                            class_ch.append(f"c{ci}:size{len(cls)}:chi{chi[cls[0]]:.8g}")
                        line=f"{lev} {sh.E:.17g} {gdeg:.17g} {mult_ir:.17g} {ir['dim']} {len(sh.rows)} {ex[0].s()} {ex[1].s()} {ex[2].s()} {zratio} [{' '.join(class_ch)}]\n"
                        files[ir['name']].write(line)
                        allf.write(f"{P.tag()} {gname}{'^D' if use_dc else ''} {ir['name']} "+line)
                        if gdeg>cfg.irrep_tol:
                            print(f"[nonint-irrep] P={P.s()} group={gname}{'^D' if use_dc else ''} irrep={ir['name']} level={lev} Ecm={sh.E:.17g} group_deg={gdeg:.8g} mult={mult_ir:.8g} irrep_dim={ir['dim']} raw_mult={len(sh.rows)} example n1={ex[0].s()} n2={ex[1].s()} n3={ex[2].s()}", flush=True)
            for f in files.values(): f.close()
            print(f"[P] complete P={P.s()} elapsed_sec={time.time()-t0:.3f}", flush=True)
    print(f"[done] wrote combined file {allpath}", flush=True)

if __name__=='__main__':
    main()
