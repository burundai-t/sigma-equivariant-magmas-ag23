#!/usr/bin/env python3
from __future__ import annotations
import csv
from collections import Counter
from itertools import product
from pathlib import Path
from typing import Callable, Dict, List, Tuple

S=(0,1,2)
M=tuple((r,c) for r in S for c in S)
PATTERNS=("RRR","SRR","RRS","RSR","DIST")
OUTDIR=Path(__file__).resolve().parent
Element=Tuple[int,int]
Diagonal=Tuple[int,int,int]
RuleFunc=Callable[[int,int,int,int],int]

def comp(a:int,b:int)->int:
    if a==b: raise ValueError('comp equal')
    return 3-a-b  # valid for distinct elements of {0,1,2}

def delta_value(d:Diagonal,r:int,c:int)->int:
    return (d[(c-r)%3]+r)%3

def g1(r1,c1,r2,c2): return r2
def g2(r1,c1,r2,c2): return c2
def g3(r1,c1,r2,c2): return comp(c1,c2) if c1!=c2 else c1
def g4(r1,c1,r2,c2): return c1
def g5(r1,c1,r2,c2): return r1
def g6(r1,c1,r2,c2): return comp(r2,c2) if r2!=c2 else r2
RULES={"g1":g1,"g2":g2,"g3":g3,"g4":g4,"g5":g5,"g6":g6}

def mul(x:Element,y:Element,rule:RuleFunc,d:Diagonal)->Element:
    r1,c1=x; r2,c2=y
    if r1==r2:
        if c1==c2: return (r1,delta_value(d,r1,c1))
        return (r1,comp(c1,c2))
    return (r1,rule(r1,c1,r2,c2))

def pattern(x,y,z):
    r1,r2,r3=x[0],y[0],z[0]
    if r1==r2==r3: return "RRR"
    if r1!=r2 and r2==r3: return "SRR"
    if r1==r2 and r2!=r3: return "RRS"
    if r1==r3 and r1!=r2: return "RSR"
    return "DIST"

def assoc_counts(rule:RuleFunc,d:Diagonal)->Dict[str,int]:
    counts={p:0 for p in PATTERNS}
    for x,y,z in product(M,repeat=3):
        if mul(mul(x,y,rule,d),z,rule,d)==mul(x,mul(y,z,rule,d),rule,d):
            counts[pattern(x,y,z)]+=1
    counts["Assoc"]=sum(counts[p] for p in PATTERNS)
    return counts

def diag_params(d:Diagonal)->Dict[str,int]:
    pdist=len(set(d)); ds=sum(1 for k,dk in enumerate(d) if dk==k)
    fix=fixdiag=fixoff=0
    for r,c in M:
        if delta_value(d,r,c)==c:
            fix+=1
            if r==c: fixdiag+=1
            else: fixoff+=1
    return {"pdist":pdist,"ds":ds,"fix":fix,"fixdiag":fixdiag,"fixoff":fixoff}

def formula_rrr(params):
    if params["pdist"]==1: return 57
    if params["pdist"]==2: return 33 if params["ds"]==0 else 39
    if params["pdist"]==3: return 27
    raise AssertionError(params)

def formula_counts(rule_name,params):
    fix=params["fix"]; fixdiag=params["fixdiag"]; fixoff=params["fixoff"]
    out={"RRR":formula_rrr(params)}
    if rule_name=="g1": out.update({"RRS":9*fixoff,"RSR":9*fixoff,"SRR":162,"DIST":0})
    elif rule_name=="g2": out.update({"RRS":6*fix,"RSR":6*fix,"SRR":6*fix,"DIST":162})
    elif rule_name=="g3": out.update({"RRS":72-2*fix,"RSR":54,"SRR":54,"DIST":54})
    elif rule_name=="g4": out.update({"RRS":162,"RSR":6*fix,"SRR":162,"DIST":162})
    elif rule_name=="g5": out.update({"RRS":18*fixdiag,"RSR":18*fixdiag,"SRR":162,"DIST":162})
    elif rule_name=="g6": out.update({"RRS":6*fix,"RSR":81-3*fix,"SRR":6*fix,"DIST":54})
    else: raise AssertionError(rule_name)
    out["Assoc_formula"]=sum(out[p] for p in PATTERNS)
    return out

def write_csv(path, rows, fieldnames):
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)

def main():
    validations=[]; mismatches=[]; counter=Counter(); by_rule={r:Counter() for r in RULES}; witnesses={}
    for rule_name,rule in RULES.items():
        for d in product(S,repeat=3):
            d=tuple(d); params=diag_params(d); enum=assoc_counts(rule,d); formula=formula_counts(rule_name,params)
            assoc=enum["Assoc"]
            match=(assoc==formula["Assoc_formula"] and all(enum[p]==formula[p] for p in PATTERNS))
            row={"rule":rule_name,"d0":d[0],"d1":d[1],"d2":d[2],**params,
                 "RRR":enum["RRR"],"SRR":enum["SRR"],"RRS":enum["RRS"],"RSR":enum["RSR"],"DIST":enum["DIST"],"Assoc":assoc,
                 "RRR_formula":formula["RRR"],"SRR_formula":formula["SRR"],"RRS_formula":formula["RRS"],"RSR_formula":formula["RSR"],"DIST_formula":formula["DIST"],"Assoc_formula":formula["Assoc_formula"],"formula_match":int(match),
                 "is_standard_diagonal":int(d==(0,0,0)),"is_anchor_27_multiple":int(assoc%27==0),"anchor_multiplier":assoc//27 if assoc%27==0 else ""}
            validations.append(row)
            if not match: mismatches.append(row)
            counter[assoc]+=1; by_rule[rule_name][assoc]+=1; witnesses.setdefault(assoc,row.copy())
    values=sorted(counter); anchors=[v for v in values if v%27==0]
    minv,maxv=min(values),max(values)
    spectrum_rows=[]
    for v in values:
        w=witnesses[v]
        spectrum_rows.append({"assoc_value":v,"normalized_assoc":f"{v}/729","decimal_assoc":f"{v/729:.9f}","multiplicity_in_Omega_elem":counter[v],"is_minimum":int(v==minv),"is_maximum":int(v==maxv),"is_multiple_of_3":int(v%3==0),"is_multiple_of_27_anchor":int(v%27==0),"anchor_multiplier":v//27 if v%27==0 else "","witness_rule":w["rule"],"witness_diagonal":f"({w['d0']},{w['d1']},{w['d2']})","witness_pdist":w["pdist"],"witness_ds":w["ds"],"witness_fix":w["fix"],"witness_fixdiag":w["fixdiag"],"witness_fixoff":w["fixoff"],"witness_RRR":w["RRR"],"witness_SRR":w["SRR"],"witness_RRS":w["RRS"],"witness_RSR":w["RSR"],"witness_DIST":w["DIST"]})
    by_rule_rows=[]
    for rn in RULES:
        vals=sorted(by_rule[rn])
        by_rule_rows.append({"rule":rn,"distinct_spectrum_size":len(vals),"min_assoc":min(vals),"max_assoc":max(vals),"values":" ".join(map(str,vals)),"multiplicities":" ".join(f"{v}:{by_rule[rn][v]}" for v in vals)})
    witness_rows=[]
    for v in values:
        w=witnesses[v]
        witness_rows.append({"assoc_value":v,"rule":w["rule"],"diagonal":f"({w['d0']},{w['d1']},{w['d2']})","pdist":w["pdist"],"ds":w["ds"],"fix":w["fix"],"fixdiag":w["fixdiag"],"fixoff":w["fixoff"],"RRR":w["RRR"],"SRR":w["SRR"],"RRS":w["RRS"],"RSR":w["RSR"],"DIST":w["DIST"]})
    anchor_rows=[{"assoc_value":v,"multiplier_of_27":v//27,"multiplicity_in_Omega_elem":counter[v],"is_minimum":int(v==minv),"is_maximum":int(v==maxv),"witness_rule":witnesses[v]["rule"],"witness_diagonal":f"({witnesses[v]['d0']},{witnesses[v]['d1']},{witnesses[v]['d2']})"} for v in anchors]
    write_csv(OUTDIR/"spectrum_values.csv",spectrum_rows,["assoc_value","normalized_assoc","decimal_assoc","multiplicity_in_Omega_elem","is_minimum","is_maximum","is_multiple_of_3","is_multiple_of_27_anchor","anchor_multiplier","witness_rule","witness_diagonal","witness_pdist","witness_ds","witness_fix","witness_fixdiag","witness_fixoff","witness_RRR","witness_SRR","witness_RRS","witness_RSR","witness_DIST"])
    write_csv(OUTDIR/"spectrum_by_rule.csv",by_rule_rows,["rule","distinct_spectrum_size","min_assoc","max_assoc","values","multiplicities"])
    write_csv(OUTDIR/"spectrum_witnesses.csv",witness_rows,["assoc_value","rule","diagonal","pdist","ds","fix","fixdiag","fixoff","RRR","SRR","RRS","RSR","DIST"])
    write_csv(OUTDIR/"spectrum_full_validation.csv",validations,["rule","d0","d1","d2","pdist","ds","fix","fixdiag","fixoff","RRR","SRR","RRS","RSR","DIST","Assoc","RRR_formula","SRR_formula","RRS_formula","RSR_formula","DIST_formula","Assoc_formula","formula_match","is_standard_diagonal","is_anchor_27_multiple","anchor_multiplier"])
    write_csv(OUTDIR/"spectrum_anchor_values.csv",anchor_rows,["assoc_value","multiplier_of_27","multiplicity_in_Omega_elem","is_minimum","is_maximum","witness_rule","witness_diagonal"])
    write_csv(OUTDIR/"spectrum_mismatches.csv",mismatches,["rule","d0","d1","d2","pdist","ds","fix","fixdiag","fixoff","RRR","SRR","RRS","RSR","DIST","Assoc","RRR_formula","SRR_formula","RRS_formula","RSR_formula","DIST_formula","Assoc_formula","formula_match","is_standard_diagonal","is_anchor_27_multiple","anchor_multiplier"])
    out=["A4 spectrum artifact: Omega_elem",f"Total magmas checked: {len(validations)}",f"Formula mismatches: {len(mismatches)}",f"Spectrum size: {len(values)}",f"Spectrum min: {minv}",f"Spectrum max: {maxv}","Spectrum values: "+" ".join(map(str,values)),"Anchor values (multiples of 27): "+" ".join(map(str,anchors)),"Anchor multipliers: "+" ".join(map(str,[v//27 for v in anchors])),f"All values divisible by 3: {all(v%3==0 for v in values)}",f"Min witness: {witnesses[minv]['rule']} diagonal=({witnesses[minv]['d0']},{witnesses[minv]['d1']},{witnesses[minv]['d2']})",f"Max witness: {witnesses[maxv]['rule']} diagonal=({witnesses[maxv]['d0']},{witnesses[maxv]['d1']},{witnesses[maxv]['d2']})","By-rule ranges:"]
    for row in by_rule_rows:
        out.append(f"  {row['rule']}: size={row['distinct_spectrum_size']} range=[{row['min_assoc']},{row['max_assoc']}]")
    (OUTDIR/"spectrum_output.txt").write_text("\n".join(out)+"\n",encoding="utf-8")
    print("\n".join(out))
if __name__=="__main__": main()
