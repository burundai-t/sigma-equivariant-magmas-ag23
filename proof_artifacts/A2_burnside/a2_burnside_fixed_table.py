#!/usr/bin/env python3
"""A2 proof artifact: Burnside fixed-point table for Omega_elem.

Scope:
    Omega_elem = {g1,...,g6} x Delta_sigma, where Delta_sigma has 27
    sigma-equivariant diagonal maps.

Claim verified:
    |Omega_elem / S3| = 90 under diagonal relabelling by S3.

Important non-claim:
    This is not a proof that the 90 orbits are full magma isomorphism
    classes under all 9! bijections of the carrier M.
"""
from __future__ import annotations
import csv, itertools
from pathlib import Path

S=(0,1,2)
M=tuple((r,c) for r in S for c in S)
RULES=('g1','g2','g3','g4','g5','g6')
DELTAS=tuple(itertools.product(S, repeat=3))
OMEGA=tuple((g,d) for g in RULES for d in DELTAS)
OUT=Path(__file__).resolve().parent

def comp(a,b):
    if a==b: raise ValueError('comp on equal')
    return next(x for x in S if x!=a and x!=b)

def delta(d,r,c): return (r + d[(c-r)%3])%3

def gv(g,r1,c1,r2,c2):
    if g=='g1': return r2
    if g=='g2': return c2
    if g=='g3': return comp(c1,c2) if c1!=c2 else c1
    if g=='g4': return c1
    if g=='g5': return r1
    if g=='g6': return comp(r2,c2) if r2!=c2 else r2
    raise ValueError(g)

def prod(m,x,y):
    g,d=m; r1,c1=x; r2,c2=y
    if r1!=r2: return (r1,gv(g,r1,c1,r2,c2))
    if c1!=c2: return (r1,comp(c1,c2))
    return (r1,delta(d,r1,c1))

def table(m): return tuple(prod(m,x,y) for x in M for y in M)

def inv(p):
    q=[0,0,0]
    for i,pi in enumerate(p): q[pi]=i
    return tuple(q)

def ap(p,x): return (p[x[0]],p[x[1]])

def conj_table(m,p):
    pinv=inv(p); out=[]
    for x in M:
        for y in M:
            out.append(ap(p, prod(m, ap(pinv,x), ap(pinv,y))))
    return tuple(out)

T2M={table(m):m for m in OMEGA}

def action(m,p):
    t=conj_table(m,p)
    if t not in T2M:
        raise RuntimeError(f'outside Omega_elem: p={p}, m={m}')
    return T2M[t]

def cyc(p):
    seen=set(); parts=[]
    for i in S:
        if i in seen: continue
        cur=[]; x=i
        while x not in seen:
            seen.add(x); cur.append(x); x=p[x]
        if len(cur)>1: parts.append('('+' '.join(map(str,cur))+')')
    return ''.join(parts) or '()'

def ctype(p):
    c=cyc(p)
    if c=='()': return 'identity'
    if c.count('(')==1 and len(c.replace('(','').replace(')','').split())==3: return '3-cycle'
    return 'transposition'

def label(p):
    return {(0,1,2):'id',(1,2,0):'sigma',(2,0,1):'sigma^2',
            (1,0,2):'tau_01',(2,1,0):'tau_02',(0,2,1):'tau_12'}[p]

def order(p):
    cur=(0,1,2)
    for k in range(1,7):
        cur=tuple(p[cur[i]] for i in S)
        if cur==(0,1,2): return k
    raise AssertionError

def ds(d): return ''.join(map(str,d))
def ms(m): return f'{m[0]}:{ds(m[1])}'

def main():
    perms=list(itertools.permutations(S))
    assert len(OMEGA)==162 and len(T2M)==162

    fixed={}
    for p in perms:
        fixed[p]=[]
        for m in OMEGA:
            if action(m,p)==m:
                fixed[p].append(m)

    seen=set(); orbits=[]
    for m in sorted(OMEGA, key=ms):
        if m in seen: continue
        orb=sorted({action(m,p) for p in perms}, key=ms)
        for x in orb: seen.add(x)
        stab=[label(p) for p in perms if action(m,p)==m]
        orbits.append((orb,stab))

    burnside_sum=sum(len(v) for v in fixed.values())
    assert burnside_sum==540
    assert burnside_sum//6==90
    assert len(orbits)==90
    assert all(len(fixed[p])==162 for p in perms if ctype(p) in ('identity','3-cycle'))
    assert all(len(fixed[p])==18 for p in perms if ctype(p)=='transposition')

    with (OUT/'A2_burnside_fixed_table.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['permutation_label','permutation_tuple','cycle_notation','cycle_type','order','fixed_count','fixed_count_by_rule','fixed_diagonals','burnside_numerator_term']
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for p in sorted(perms,key=lambda q:(ctype(q),label(q))):
            by={g:0 for g in RULES}
            for g,d in fixed[p]: by[g]+=1
            w.writerow({'permutation_label':label(p),'permutation_tuple':''.join(map(str,p)),'cycle_notation':cyc(p),'cycle_type':ctype(p),'order':order(p),'fixed_count':len(fixed[p]),'fixed_count_by_rule':';'.join(f'{g}={by[g]}' for g in RULES),'fixed_diagonals':';'.join(sorted({ds(d) for _,d in fixed[p]})),'burnside_numerator_term':len(fixed[p])})

    with (OUT/'A2_burnside_fixed_magmas.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['permutation_label','cycle_type','magma','rule','delta']
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for p in sorted(perms,key=lambda q:(ctype(q),label(q))):
            for m in sorted(fixed[p],key=ms):
                w.writerow({'permutation_label':label(p),'cycle_type':ctype(p),'magma':ms(m),'rule':m[0],'delta':ds(m[1])})

    with (OUT/'A2_burnside_orbits.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['orbit_id','orbit_size','representative','members','stabilizer_size','stabilizer']
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for i,(orb,stab) in enumerate(orbits,1):
            w.writerow({'orbit_id':i,'orbit_size':len(orb),'representative':ms(orb[0]),'members':';'.join(ms(x) for x in orb),'stabilizer_size':len(stab),'stabilizer':';'.join(stab)})

    hist={s:sum(1 for o,_ in orbits if len(o)==s) for s in sorted({len(o) for o,_ in orbits})}
    lines=['A2 Burnside fixed-point verification for Omega_elem',
           '===================================================',
           f'|Omega_elem| = {len(OMEGA)}','|S3| = 6','Fixed counts:']
    for p in sorted(perms,key=lambda q:(ctype(q),label(q))):
        lines.append(f'  {label(p):7s} {"".join(map(str,p))} {cyc(p):7s} {ctype(p):13s} fixed={len(fixed[p])}')
    lines += [f'Burnside numerator = {burnside_sum}',
              f'Orbit count = {burnside_sum}/6 = {burnside_sum//6}',
              f'Explicit orbit traversal count = {len(orbits)}',
              f'Orbit-size histogram = {hist}',
              'Scope: diagonal S3-orbits inside Omega_elem; not full isomorphism classes.',
              '']
    output='\n'.join(lines)
    (OUT/'A2_burnside_count_output.txt').write_text(output,encoding='utf-8')
    print(output)

if __name__ == '__main__':
    main()
