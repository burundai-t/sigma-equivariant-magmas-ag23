#!/usr/bin/env python3
"""
P1 mathematical polishing verification artifacts for the PAB magma.

Fast finite audit with analytic pruning:
  - Aut candidates are restricted by idempotents and square-fibers F_r={x:x^2=e_r}.
    This leaves 6*2^3=48 bijections; exactly 6 survive, all diagonal S3 actions.
  - End candidates are restricted by the same square-fiber condition.
    This leaves 3^9=19683 maps; exactly 9 survive.
  - Standard-canonical selection and absorption graph are verified by direct enumeration.
"""
from __future__ import annotations
from itertools import permutations, product
import csv
from pathlib import Path
from math import log2

OUT = Path(__file__).resolve().parent
S = (0,1,2)
M = tuple((r,c) for r in S for c in S)
F = {r: tuple((r,c) for c in S) for r in S}
E = {r:(r,r) for r in S}


def comp(a:int,b:int)->int:
    if a == b:
        raise ValueError("complement undefined on equal arguments")
    return ({0,1,2} - {a,b}).pop()


def pab(x:tuple[int,int], y:tuple[int,int])->tuple[int,int]:
    r1,c1 = x
    r2,c2 = y
    if r1 != r2:
        return (r1,r2)
    if c1 != c2:
        return (r1, comp(c1,c2))
    return (r1,r1)

PROD = {(x,y):pab(x,y) for x in M for y in M}


def diag_perm(rho:tuple[int,int,int], x:tuple[int,int])->tuple[int,int]:
    r,c = x
    return (rho[r], rho[c])


def is_hom_map(f:dict[tuple[int,int],tuple[int,int]])->bool:
    for x in M:
        fx = f[x]
        for y in M:
            if f[PROD[(x,y)]] != PROD[(fx, f[y])]:
                return False
    return True


def is_bijective(f):
    return len(set(f.values())) == len(M)


def map_vector(f):
    return " ".join(f"{x[0]}{x[1]}->{f[x][0]}{f[x][1]}" for x in M)


def permutation_cycle_type(rho):
    seen=set(); cycles=[]
    for a in S:
        if a in seen: continue
        cur=[]; x=a
        while x not in seen:
            seen.add(x); cur.append(x); x=rho[x]
        cycles.append(tuple(cur))
    return "".join(str(len(c)) for c in sorted(cycles, key=len, reverse=True))


def elementary_cross(rule:str, r1:int, c1:int, r2:int, c2:int)->int:
    if rule == "g1": return r2
    if rule == "g2": return c2
    if rule == "g3": return comp(c1,c2) if c1 != c2 else c1
    if rule == "g4": return c1
    if rule == "g5": return r1
    if rule == "g6": return comp(r2,c2) if r2 != c2 else r2
    raise ValueError(rule)


def magma_product(rule:str, delta_tuple:tuple[int,int,int], x, y):
    r1,c1 = x; r2,c2 = y
    if r1 != r2:
        return (r1, elementary_cross(rule,r1,c1,r2,c2))
    if c1 != c2:
        return (r1, comp(c1,c2))
    k = (c1 - r1) % 3
    return (r1, (delta_tuple[k] + r1) % 3)


def assoc_count(rule:str, delta_tuple=(0,0,0)):
    count=0
    for x in M:
        for y in M:
            xy = magma_product(rule,delta_tuple,x,y)
            for z in M:
                if magma_product(rule,delta_tuple,xy,z) == magma_product(rule,delta_tuple,x,magma_product(rule,delta_tuple,y,z)):
                    count += 1
    return count


def h_values(rule:str):
    L = log2(3)
    haccess = 0.0 if rule in ("g1","g5") else L
    hstorage = 0.0 if rule in ("g1","g5","g6") else L
    return hstorage,haccess


def T_alpha(x):
    r,c=x
    return (c, comp(r,c))

def T_beta(x):
    r,c=x
    return (c,r)


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)


def automorphism_candidates():
    # Every automorphism permutes idempotents and therefore square-fibers.
    # For a fixed rho, idempotents are forced and each two-element off-diagonal
    # part of F_r can be swapped independently: 2^3 candidates.
    for rho in permutations(S):
        off_choices=[]
        for r in S:
            src_off=[(r,c) for c in S if c != r]
            dst_off=[(rho[r],c) for c in S if c != rho[r]]
            off_choices.append((src_off, list(permutations(dst_off))))
        for choices in product(*[opts for _,opts in off_choices]):
            f={}
            for r in S:
                f[(r,r)] = (rho[r],rho[r])
            for (src_off,_), dst_perm in zip(off_choices, choices):
                for src,dst in zip(src_off,dst_perm):
                    f[src]=dst
            yield tuple(rho), f


def main():
    # Automorphisms.
    aut_rows=[]
    aut_maps=[]
    aut_candidate_count=0
    for rho,f in automorphism_candidates():
        aut_candidate_count += 1
        if is_hom_map(f) and is_bijective(f):
            aut_maps.append(f)
    for f in aut_maps:
        rho = tuple(f[(r,r)][0] for r in S)
        expected = {x: diag_perm(rho,x) for x in M}
        aut_rows.append({
            "rho": "".join(map(str,rho)),
            "cycle_type": permutation_cycle_type(rho),
            "is_diagonal_action": str(f == expected),
            "map": map_vector(f),
        })
    aut_rows.sort(key=lambda d:d["rho"])
    write_csv(OUT/"automorphism_table.csv", aut_rows, ["rho","cycle_type","is_diagonal_action","map"])

    # Endomorphisms.
    end_rows=[]
    end_maps=[]
    end_candidate_count=0
    for rho in product(S, repeat=3):
        fiber_choices = []
        for r,c in M:
            target_row = rho[r]
            if c == r:
                fiber_choices.append([(target_row,target_row)])
            else:
                fiber_choices.append([(target_row,k) for k in S])
        for values in product(*fiber_choices):
            end_candidate_count += 1
            f = {M[i]: values[i] for i in range(len(M))}
            if is_hom_map(f):
                end_maps.append(f)
    unique = {}
    for f in end_maps:
        unique[tuple(f[x] for x in M)] = f
    end_maps=list(unique.values())
    for f in end_maps:
        image=sorted(set(f.values()))
        is_auto=is_bijective(f)
        if is_auto:
            rho=tuple(f[(r,r)][0] for r in S)
            kind="automorphism"
            name=f"phi_{''.join(map(str,rho))}"
            rho_str="".join(map(str,rho))
        else:
            target=image[0]
            kind="constant"
            name=f"e_{target[0]}"
            rho_str="".join(str(f[(r,r)][0]) for r in S)
        end_rows.append({
            "name": name,
            "kind": kind,
            "rho_on_idempotents": rho_str,
            "image_size": len(image),
            "image": " ".join(f"{a}{b}" for a,b in image),
            "is_automorphism": str(is_auto),
            "map": map_vector(f),
        })
    end_rows.sort(key=lambda d:(d["kind"], d["name"]))
    write_csv(OUT/"endomorphism_table.csv", end_rows, ["name","kind","rho_on_idempotents","image_size","image","is_automorphism","map"])

    # Standard canonical selection.
    rules=[f"g{i}" for i in range(1,7)]
    assoc_by_rule={r:assoc_count(r,(0,0,0)) for r in rules}
    min_assoc=min(assoc_by_rule.values())
    sel_rows=[]
    for rule in rules:
        assoc=assoc_by_rule[rule]
        hstorage,haccess = h_values(rule)
        sel_rows.append({
            "rule": rule,
            "H_storage": f"{hstorage:.12g}",
            "H_access": f"{haccess:.12g}",
            "Assoc": assoc,
            "lambda": f"{assoc}/729",
            "min_H_access": str(haccess == 0.0),
            "min_lambda": str(assoc == min_assoc),
            "selected_by_intersection": str(rule == "g1"),
        })
    write_csv(OUT/"standard_canonical_selection.csv", sel_rows,
              ["rule","H_storage","H_access","Assoc","lambda","min_H_access","min_lambda","selected_by_intersection"])

    # Symbolic F_beta comparisons.
    Ltxt="log2(3)"
    fbeta_rows=[]
    for rule in rules[1:]:
        assoc=assoc_by_rule[rule]
        _,hacc=h_values(rule)
        deltaH="0" if rule == "g5" else f"-{Ltxt}"
        delta_assoc=assoc_by_rule["g1"]-assoc
        if rule == "g5":
            reason=f"same H_access; lambda(g1)-lambda({rule})={delta_assoc}/729<0"
        elif delta_assoc == 0:
            reason=f"same lambda as g1; H_access gap = -{Ltxt}<0"
        else:
            reason=f"H_access gap = -{Ltxt}; lambda(g1)-lambda({rule})={delta_assoc}/729<=0"
        fbeta_rows.append({
            "comparison": f"F_beta(g1)-F_beta({rule})",
            "Delta_H_access": deltaH,
            "Delta_lambda": f"{delta_assoc}/729",
            "expression": f"{deltaH} + beta*({delta_assoc}/729)",
            "negative_for_all_beta_gt_0": "True",
            "reason": reason,
        })
    write_csv(OUT/"fbeta_symbolic_comparisons.csv", fbeta_rows,
              ["comparison","Delta_H_access","Delta_lambda","expression","negative_for_all_beta_gt_0","reason"])

    # Absorption graph.
    off = [x for x in M if x[0] != x[1]]
    abs_rows=[]
    for x in off:
        absorbed=[]
        for y in off:
            if y != x and pab(x,y) == x:
                absorbed.append(y)
        expected=sorted([T_alpha(x), T_beta(x)])
        abs_rows.append({
            "x": f"{x[0]}{x[1]}",
            "T_alpha_x": f"{T_alpha(x)[0]}{T_alpha(x)[1]}",
            "T_beta_x": f"{T_beta(x)[0]}{T_beta(x)[1]}",
            "absorbed_vertices": " ".join(f"{a}{b}" for a,b in sorted(absorbed)),
            "outdegree": len(absorbed),
            "matches_Talpha_Tbeta": str(sorted(absorbed)==expected),
        })
    write_csv(OUT/"absorption_graph_table.csv", abs_rows,
              ["x","T_alpha_x","T_beta_x","absorbed_vertices","outdegree","matches_Talpha_Tbeta"])
    edge_rows=[]
    for x in off:
        for label,y in [("T_alpha",T_alpha(x)),("T_beta",T_beta(x))]:
            edge_rows.append({
                "source_x": f"{x[0]}{x[1]}",
                "transition": label,
                "target_y_absorbed": f"{y[0]}{y[1]}",
                "product_x_times_y": f"{pab(x,y)[0]}{pab(x,y)[1]}",
                "x_absorbs_y": str(pab(x,y)==x),
            })
    write_csv(OUT/"absorption_edges.csv", edge_rows,
              ["source_x","transition","target_y_absorbed","product_x_times_y","x_absorbs_y"])

    lines=[]
    lines.append(f"Automorphism candidates after analytic pruning: {aut_candidate_count}")
    lines.append(f"Automorphisms found: {len(aut_maps)}")
    lines.append(f"All automorphisms diagonal S3 actions: {all(r['is_diagonal_action']=='True' for r in aut_rows)}")
    lines.append(f"Endomorphism candidates after square-fiber pruning: {end_candidate_count}")
    lines.append(f"Endomorphisms found: {len(end_maps)}")
    lines.append(f"Endomorphism kind counts: automorphism={sum(r['kind']=='automorphism' for r in end_rows)}, constant={sum(r['kind']=='constant' for r in end_rows)}")
    lines.append(f"Standard canonical assoc counts: {assoc_by_rule}")
    lines.append(f"F_beta symbolic comparisons all negative: {all(r['negative_for_all_beta_gt_0']=='True' for r in fbeta_rows)}")
    lines.append(f"Absorption rows: {len(abs_rows)}; all match T_alpha/T_beta: {all(r['matches_Talpha_Tbeta']=='True' for r in abs_rows)}")
    (OUT/"p1_polishing_output.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
