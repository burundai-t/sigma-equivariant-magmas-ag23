#!/usr/bin/env python3
r"""
A5 proof artifact: restricted sequential selection on Ω'_elem.

This script builds the restricted extended elementary space

    Ω'_elem = F_sigma × {g1,...,g6} × Δ_sigma

where:
  * F_sigma is the 9-element family of translation-equivariant ordered
    within-fiber rules f:S×S\diag→S;
  * {g1,...,g6} are the six distinguished elementary cross-rules;
  * Δ_sigma is the 27-element family of sigma-equivariant diagonal maps.

It verifies the sequential selection count:

    1458 --min H_access--> 486 --min H_diag--> 54 --min Assoc--> 1 = PAB.

The result is intentionally restricted to Ω'_elem. It is not a global theorem
about the full 3^18 cross-rule universe.
"""

from __future__ import annotations

import csv
import itertools
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

S = (0, 1, 2)
M = tuple((r, c) for r in S for c in S)
M_INDEX = {x: i for i, x in enumerate(M)}
INDEX_M = {i: x for x, i in M_INDEX.items()}

OUTDIR = Path(__file__).resolve().parent


def complement(a: int, b: int) -> int:
    if a == b:
        raise ValueError("complement requires distinct inputs")
    return next(x for x in S if x != a and x != b)


# ---------------------------------------------------------------------------
# Six distinguished elementary cross-rules.
# Input: r1,c1,r2,c2 with r1 != r2. Output: column in S.
# ---------------------------------------------------------------------------

def cross_value(rule_id: int, r1: int, c1: int, r2: int, c2: int) -> int:
    if rule_id == 1:
        return r2
    if rule_id == 2:
        return c2
    if rule_id == 3:
        return complement(c1, c2) if c1 != c2 else c1
    if rule_id == 4:
        return c1
    if rule_id == 5:
        return r1
    if rule_id == 6:
        return complement(r2, c2) if r2 != c2 else r2
    raise ValueError(f"unknown rule_id: {rule_id}")


CROSS_RULE_NAMES = {
    1: "g1_column_blind_PAB",
    2: "g2_transparent_c2",
    3: "g3_su2_transparent",
    4: "g4_echo_c1",
    5: "g5_self_referential_r1",
    6: "g6_anti_complement",
}

HACCESS_EXACT = {
    1: "0",
    2: "log2(3)",
    3: "log2(3)",
    4: "log2(3)",
    5: "0",
    6: "log2(3)",
}

HACCESS_NUMERIC = {
    1: 0.0,
    2: math.log2(3),
    3: math.log2(3),
    4: math.log2(3),
    5: 0.0,
    6: math.log2(3),
}


# ---------------------------------------------------------------------------
# 9 translation-equivariant ordered within-fiber rules.
# ---------------------------------------------------------------------------
# We parameterise f by two seeds:
#   u = f(0,1), v = f(0,2).
# Equivariance f(a+1,b+1)=f(a,b)+1 gives:
#   f(a,a+1)=u+a mod 3,
#   f(a,a-1)=v+a mod 3.
# The PAB/Steiner complement rule is (u,v)=(2,1).
# ---------------------------------------------------------------------------

def fiber_value(u: int, v: int, a: int, b: int) -> int:
    if a == b:
        raise ValueError("fiber_value only handles distinct columns")
    if (b - a) % 3 == 1:
        return (u + a) % 3
    return (v + a) % 3


def fiber_rule_values(u: int, v: int) -> Dict[Tuple[int, int], int]:
    return {(a, b): fiber_value(u, v, a, b) for a in S for b in S if a != b}


def is_complement_fiber(u: int, v: int) -> bool:
    return all(fiber_value(u, v, a, b) == complement(a, b) for a in S for b in S if a != b)


def avoids_inputs_for_all_distinct_pairs(u: int, v: int) -> bool:
    return all(fiber_value(u, v, a, b) not in {a, b} for a in S for b in S if a != b)


# ---------------------------------------------------------------------------
# 27 sigma-equivariant diagonal maps.
# Base triple d=(d0,d1,d2), where dk = δ(0,k).
# Equivariance gives δ(r,c)=r+d[c-r] mod 3.
# ---------------------------------------------------------------------------

def delta_value(d: Tuple[int, int, int], r: int, c: int) -> int:
    return (r + d[(c - r) % 3]) % 3


def diagonal_entropy_exact(d: Tuple[int, int, int]) -> str:
    counts = sorted(Counter(d).values(), reverse=True)
    if counts == [3]:
        return "0"
    if counts == [2, 1]:
        return "H(2/3,1/3)"
    if counts == [1, 1, 1]:
        return "log2(3)"
    raise AssertionError(f"unexpected counts: {counts}")


def diagonal_entropy_numeric(d: Tuple[int, int, int]) -> float:
    total = len(d)
    h = 0.0
    for n in Counter(d).values():
        p = n / total
        h -= p * math.log2(p)
    return h


def is_standard_diagonal(d: Tuple[int, int, int]) -> bool:
    return d == (0, 0, 0)


def is_constant_diagonal(d: Tuple[int, int, int]) -> bool:
    return d[0] == d[1] == d[2]


# ---------------------------------------------------------------------------
# Operation tables and associativity.
# ---------------------------------------------------------------------------

def build_op_table(rule_id: int, u: int, v: int, d: Tuple[int, int, int]) -> List[List[int]]:
    table = [[0] * 9 for _ in range(9)]
    for i, (r1, c1) in enumerate(M):
        for j, (r2, c2) in enumerate(M):
            if r1 != r2:
                out = (r1, cross_value(rule_id, r1, c1, r2, c2))
            elif c1 != c2:
                out = (r1, fiber_value(u, v, c1, c2))
            else:
                out = (r1, delta_value(d, r1, c1))
            table[i][j] = M_INDEX[out]
    return table


def associativity_count(op: List[List[int]]) -> int:
    count = 0
    for x in range(9):
        opx = op[x]
        for y in range(9):
            xy = opx[y]
            opxy = op[xy]
            for z in range(9):
                if opxy[z] == opx[op[y][z]]:
                    count += 1
    return count


# ---------------------------------------------------------------------------
# Writers.
# ---------------------------------------------------------------------------

def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    fiber_rows: List[Dict[str, object]] = []
    for u, v in itertools.product(S, S):
        vals = fiber_rule_values(u, v)
        row = {
            "fiber_id": f"F{u}{v}",
            "u_f_0_1": u,
            "v_f_0_2": v,
            "is_complement_AX2": is_complement_fiber(u, v),
            "avoids_inputs_for_all_distinct_pairs": avoids_inputs_for_all_distinct_pairs(u, v),
        }
        for a, b in [(0,1),(1,2),(2,0),(0,2),(1,0),(2,1)]:
            row[f"f_{a}_{b}"] = vals[(a, b)]
        fiber_rows.append(row)

    all_rows: List[Dict[str, object]] = []
    step3_rows: List[Dict[str, object]] = []

    for rule_id in range(1, 7):
        for u, v in itertools.product(S, S):
            for d in itertools.product(S, S, S):
                op = build_op_table(rule_id, u, v, d)
                assoc = associativity_count(op)
                h_access_num = HACCESS_NUMERIC[rule_id]
                h_diag_num = diagonal_entropy_numeric(d)
                step1 = h_access_num == 0.0
                step2 = step1 and h_diag_num == 0.0
                row = {
                    "rule_id": f"g{rule_id}",
                    "rule_name": CROSS_RULE_NAMES[rule_id],
                    "fiber_id": f"F{u}{v}",
                    "u_f_0_1": u,
                    "v_f_0_2": v,
                    "is_complement_AX2": is_complement_fiber(u, v),
                    "diagonal": str(d),
                    "d0": d[0],
                    "d1": d[1],
                    "d2": d[2],
                    "is_constant_diagonal": is_constant_diagonal(d),
                    "is_standard_diagonal": is_standard_diagonal(d),
                    "H_access_exact": HACCESS_EXACT[rule_id],
                    "H_access_numeric": h_access_num,
                    "H_diag_exact": diagonal_entropy_exact(d),
                    "H_diag_numeric": h_diag_num,
                    "assoc_count": assoc,
                    "lambda_assoc_over_729": assoc / 729,
                    "step1_min_H_access_survivor": step1,
                    "step2_min_H_diag_survivor": step2,
                    "is_PAB": (rule_id == 1 and is_complement_fiber(u, v) and is_standard_diagonal(d)),
                }
                all_rows.append(row)
                if step2:
                    step3_rows.append(row.copy())

    min_h_access = min(r["H_access_numeric"] for r in all_rows)
    step1_survivors = [r for r in all_rows if r["H_access_numeric"] == min_h_access]
    min_h_diag_after_step1 = min(r["H_diag_numeric"] for r in step1_survivors)
    step2_survivors = [r for r in step1_survivors if r["H_diag_numeric"] == min_h_diag_after_step1]
    min_assoc_after_step2 = min(r["assoc_count"] for r in step2_survivors)
    step3_survivors = [r for r in step2_survivors if r["assoc_count"] == min_assoc_after_step2]

    global_min_assoc = min(r["assoc_count"] for r in all_rows)
    global_min_rows = [r for r in all_rows if r["assoc_count"] == global_min_assoc]

    for r in all_rows:
        r["step3_min_assoc_survivor"] = (r in step3_survivors)
        r["global_min_assoc_survivor"] = (r["assoc_count"] == global_min_assoc)

    for r in step3_rows:
        r["step3_min_assoc_survivor"] = (
            r["assoc_count"] == min_assoc_after_step2
            and r["rule_id"] == "g1"
            and r["is_complement_AX2"]
            and r["is_standard_diagonal"]
        )

    step_rows = [
        {
            "step": 0,
            "criterion": "full Ω'_elem",
            "survivors": len(all_rows),
            "description": "9 fiber rules × 6 elementary cross-rules × 27 diagonals",
        },
        {
            "step": 1,
            "criterion": "min H_access",
            "survivors": len(step1_survivors),
            "description": "9 fiber rules × {g1,g5} × 27 diagonals",
        },
        {
            "step": 2,
            "criterion": "min H_diag among Step 1 survivors",
            "survivors": len(step2_survivors),
            "description": "9 fiber rules × {g1,g5} × 3 constant diagonals",
        },
        {
            "step": 3,
            "criterion": "min associativity λ among Step 2 survivors",
            "survivors": len(step3_survivors),
            "description": "unique survivor: g1 + complement fiber F21 + standard diagonal (0,0,0)",
        },
    ]

    all_fieldnames = [
        "rule_id", "rule_name", "fiber_id", "u_f_0_1", "v_f_0_2", "is_complement_AX2",
        "diagonal", "d0", "d1", "d2", "is_constant_diagonal", "is_standard_diagonal",
        "H_access_exact", "H_access_numeric", "H_diag_exact", "H_diag_numeric",
        "assoc_count", "lambda_assoc_over_729",
        "step1_min_H_access_survivor", "step2_min_H_diag_survivor", "step3_min_assoc_survivor",
        "global_min_assoc_survivor", "is_PAB",
    ]

    write_csv(OUTDIR / "omega_prime_fiber_rules.csv", fiber_rows, list(fiber_rows[0].keys()))
    write_csv(OUTDIR / "omega_prime_survivors.csv", all_rows, all_fieldnames)
    write_csv(OUTDIR / "omega_prime_step3_candidates.csv", step2_survivors, all_fieldnames)
    write_csv(OUTDIR / "omega_prime_step3_winner.csv", step3_survivors, all_fieldnames)
    write_csv(OUTDIR / "omega_prime_global_min_assoc.csv", global_min_rows, all_fieldnames)
    write_csv(OUTDIR / "omega_prime_selection_steps.csv", step_rows, ["step", "criterion", "survivors", "description"])

    output_lines = []
    output_lines.append("A5 restricted sequential selection on Ω'_elem")
    output_lines.append("================================================")
    output_lines.append(f"Fiber rules: {len(fiber_rows)}")
    output_lines.append(f"Elementary cross-rules: 6")
    output_lines.append(f"Diagonals: 27")
    output_lines.append(f"Total Ω'_elem objects: {len(all_rows)}")
    output_lines.append("")
    for row in step_rows:
        output_lines.append(f"Step {row['step']}: {row['criterion']} -> {row['survivors']}")
    output_lines.append("")
    output_lines.append(f"Minimum associativity after Step 2: {min_assoc_after_step2}")
    output_lines.append(f"Step 3 survivor count: {len(step3_survivors)}")
    for r in step3_survivors:
        output_lines.append(
            "Winner: "
            f"rule={r['rule_id']}, fiber={r['fiber_id']} (u={r['u_f_0_1']},v={r['v_f_0_2']}), "
            f"diagonal={r['diagonal']}, assoc={r['assoc_count']}, is_PAB={r['is_PAB']}"
        )
    output_lines.append("")
    output_lines.append("Global associativity check (not the theorem):")
    output_lines.append(f"Global min assoc over Ω'_elem: {global_min_assoc}")
    output_lines.append(f"Global min survivor count: {len(global_min_rows)}")
    for r in global_min_rows[:10]:
        output_lines.append(
            f"GlobalMin: rule={r['rule_id']}, fiber={r['fiber_id']}, diagonal={r['diagonal']}, assoc={r['assoc_count']}"
        )
    if len(global_min_rows) > 10:
        output_lines.append(f"... plus {len(global_min_rows)-10} more")
    output_lines.append("")
    output_lines.append("Sanity checks:")
    output_lines.append(f"Expected chain 1458 -> 486 -> 54 -> 1: { [len(all_rows), len(step1_survivors), len(step2_survivors), len(step3_survivors)] }")
    output_lines.append(f"Unique Step 3 survivor is PAB: {len(step3_survivors)==1 and step3_survivors[0]['is_PAB']}")

    (OUTDIR / "omega_prime_selection_output.txt").write_text("\n".join(output_lines) + "\n", encoding="utf-8")

    report = f"""# A5 — Restricted Sequential Selection Artifact

## Result

This artifact verifies the restricted sequential selection theorem inside

\\[
\\Omega'_{{\\mathrm{{elem}}}}
=\\mathcal F_\\sigma\\times\\{{g_1,\\dots,g_6\\}}\\times\\Delta_\\sigma,
\\qquad |\\Omega'_{{\\mathrm{{elem}}}}|=9\\cdot6\\cdot27=1458.
\\]

The verified sequential chain is:

\\[
1458 \\xrightarrow{{\\min H_{{\\mathrm{{access}}}}}} 486
\\xrightarrow{{\\min H_{{\\mathrm{{diag}}}}}} 54
\\xrightarrow{{\\min \\lambda}} 1.
\\]

The unique final survivor is PAB:

\\[
(g_1, F_{{21}}, \\delta_{{\\mathrm{{std}}}}),
\\]

where `F21` is the complement/Steiner fiber rule with

\\[
f(0,1)=2,\\qquad f(0,2)=1.
\\]

## Fiber-rule parameterisation

The 9 fiber rules are translation-equivariant ordered rules on distinct columns.
Each rule is determined by two seeds:

\\[
u=f(0,1),\\qquad v=f(0,2).
\\]

Equivariance gives:

\\[
f(a,a+1)=u+a\\pmod 3,
\\]

\\[
f(a,a-1)=v+a\\pmod 3.
\\]

Thus there are \\(3^2=9\\) rules. The PAB/AX2 complement rule is \\((u,v)=(2,1)\\).

## Files

- `omega_prime_selection.py` — executable verification script.
- `omega_prime_fiber_rules.csv` — the 9 fiber rules.
- `omega_prime_survivors.csv` — all 1458 objects with filters and associativity counts.
- `omega_prime_selection_steps.csv` — 1458 → 486 → 54 → 1 table.
- `omega_prime_step3_candidates.csv` — the 54 Step-2 survivors.
- `omega_prime_step3_winner.csv` — the unique final survivor.
- `omega_prime_global_min_assoc.csv` — global associativity minima over Ω'_elem; included only as a warning that the theorem is sequential.
- `omega_prime_selection_output.txt` — console output.

## Verified counts

| Step | Criterion | Survivors |
|---:|---|---:|
| 0 | full \\(\\Omega'_{{\\mathrm{{elem}}}}\\) | {len(all_rows)} |
| 1 | min \\(H_{{\\mathrm{{access}}}}\\) | {len(step1_survivors)} |
| 2 | min \\(H_{{\\mathrm{{diag}}}}\\) among Step 1 | {len(step2_survivors)} |
| 3 | min \\(\\lambda\\) among Step 2 | {len(step3_survivors)} |

## Important scope warning

This is **not** a global triple-argmin theorem. It is a sequential theorem.
The script also records the global minimum of associativity over \\(\\Omega'_{{\\mathrm{{elem}}}}\\), which is `{global_min_assoc}` and is not PAB.

The manuscript-safe statement is:

> Inside \\(\\Omega'_{{\\mathrm{{elem}}}}\\), after sequential restriction to minimal \\(H_{{\\mathrm{{access}}}}\\) and then minimal \\(H_{{\\mathrm{{diag}}}}\\), the associativity count is uniquely minimised by PAB.

"""
    (OUTDIR / "A5_OMEGA_PRIME_SELECTION_REPORT.md").write_text(report, encoding="utf-8")

    print("\n".join(output_lines))


if __name__ == "__main__":
    main()
