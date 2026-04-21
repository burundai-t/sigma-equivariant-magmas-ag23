#!/usr/bin/env python3
"""
A3 proof artifact: mixed-pattern formulas for Ω_elem.

Generates:
  - A3_mixed_pattern_formula_table.csv
  - A3_rrr_formula_table.csv
  - A3_mixed_pattern_validation.csv
  - A3_mixed_pattern_summary.csv
  - A3_mixed_pattern_count_output.txt

The script verifies, by exhaustive enumeration over Ω_elem = {g1,...,g6} × Δσ,
that the row-pattern associativity counts match the formula table used in the manuscript.
"""
from __future__ import annotations

import csv
import itertools
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple

S = (0, 1, 2)
M = tuple((r, c) for r in S for c in S)
PATTERNS = ("RRR", "SRR", "RRS", "RSR", "DIST")
OUT = Path(__file__).resolve().parent


def comp(a: int, b: int) -> int:
    """Complement in S={0,1,2}. Requires a != b."""
    if a == b:
        raise ValueError(f"comp({a},{b}) undefined for equal inputs")
    return next(x for x in S if x not in (a, b))


def g1(r1: int, c1: int, r2: int, c2: int) -> int:
    return r2


def g2(r1: int, c1: int, r2: int, c2: int) -> int:
    return c2


def g3(r1: int, c1: int, r2: int, c2: int) -> int:
    return comp(c1, c2) if c1 != c2 else c1


def g4(r1: int, c1: int, r2: int, c2: int) -> int:
    return c1


def g5(r1: int, c1: int, r2: int, c2: int) -> int:
    return r1


def g6(r1: int, c1: int, r2: int, c2: int) -> int:
    return comp(r2, c2) if r2 != c2 else r2


RULES: Dict[str, Callable[[int, int, int, int], int]] = {
    "g1": g1,
    "g2": g2,
    "g3": g3,
    "g4": g4,
    "g5": g5,
    "g6": g6,
}

RULE_NAMES = {
    "g1": "column-blind/PAB",
    "g2": "transparent",
    "g3": "su2-transparent",
    "g4": "echo",
    "g5": "self-referential",
    "g6": "anti-complement",
}


def sigma(x: int, k: int = 1) -> int:
    return (x + k) % 3


def delta_value(d: Tuple[int, int, int], r: int, c: int) -> int:
    """σ-equivariant diagonal map determined by d_k = δ(0,k)."""
    return (d[(c - r) % 3] + r) % 3


def op(rule: Callable[[int, int, int, int], int], d: Tuple[int, int, int], x: Tuple[int, int], y: Tuple[int, int]) -> Tuple[int, int]:
    r1, c1 = x
    r2, c2 = y
    if r1 != r2:
        return (r1, rule(r1, c1, r2, c2))
    if c1 != c2:
        return (r1, comp(c1, c2))
    return (r1, delta_value(d, r1, c1))


def pattern(x: Tuple[int, int], y: Tuple[int, int], z: Tuple[int, int]) -> str:
    r1, r2, r3 = x[0], y[0], z[0]
    if r1 == r2 == r3:
        return "RRR"
    if r1 != r2 and r2 == r3:
        return "SRR"
    if r1 == r2 and r2 != r3:
        return "RRS"
    if r1 == r3 and r1 != r2:
        return "RSR"
    if len({r1, r2, r3}) == 3:
        return "DIST"
    raise AssertionError("unreachable row pattern")


def actual_pattern_counts(rule: Callable[[int, int, int, int], int], d: Tuple[int, int, int]) -> Counter:
    counts = Counter({p: 0 for p in PATTERNS})
    totals = Counter({p: 0 for p in PATTERNS})
    for x, y, z in itertools.product(M, repeat=3):
        p = pattern(x, y, z)
        totals[p] += 1
        lhs = op(rule, d, op(rule, d, x, y), z)
        rhs = op(rule, d, x, op(rule, d, y, z))
        if lhs == rhs:
            counts[p] += 1
    assert totals == Counter({"RRR": 81, "SRR": 162, "RRS": 162, "RSR": 162, "DIST": 162})
    return counts


def diag_params(d: Tuple[int, int, int]) -> Dict[str, int | str]:
    pdist = len(set(d))
    ds = sum(1 for k, dk in enumerate(d) if dk == k)
    fixdiag = 3 if d[0] == 0 else 0
    fixoff = (3 if d[1] == 1 else 0) + (3 if d[2] == 2 else 0)
    fix = fixdiag + fixoff
    return {
        "d": "".join(map(str, d)),
        "d0": d[0],
        "d1": d[1],
        "d2": d[2],
        "pdist": pdist,
        "ds": ds,
        "fix": fix,
        "fixdiag": fixdiag,
        "fixoff": fixoff,
    }


def predicted_rrr(d: Tuple[int, int, int]) -> int:
    p = diag_params(d)
    pdist = int(p["pdist"])
    ds = int(p["ds"])
    if pdist == 1:
        return 57
    if pdist == 2 and ds == 0:
        return 33
    if pdist == 2 and ds >= 1:
        return 39
    if pdist == 3:
        return 27
    raise AssertionError(f"unhandled RRR class d={d}")


def predicted_counts(rule_name: str, d: Tuple[int, int, int]) -> Dict[str, int]:
    p = diag_params(d)
    fix = int(p["fix"])
    fixdiag = int(p["fixdiag"])
    fixoff = int(p["fixoff"])
    pred = {"RRR": predicted_rrr(d)}
    if rule_name == "g1":
        pred.update({"RRS": 9 * fixoff, "RSR": 9 * fixoff, "SRR": 162, "DIST": 0})
    elif rule_name == "g2":
        pred.update({"RRS": 6 * fix, "RSR": 6 * fix, "SRR": 6 * fix, "DIST": 162})
    elif rule_name == "g3":
        pred.update({"RRS": 72 - 2 * fix, "RSR": 54, "SRR": 54, "DIST": 54})
    elif rule_name == "g4":
        pred.update({"RRS": 162, "RSR": 6 * fix, "SRR": 162, "DIST": 162})
    elif rule_name == "g5":
        pred.update({"RRS": 18 * fixdiag, "RSR": 18 * fixdiag, "SRR": 162, "DIST": 162})
    elif rule_name == "g6":
        pred.update({"RRS": 6 * fix, "RSR": 81 - 3 * fix, "SRR": 6 * fix, "DIST": 54})
    else:
        raise AssertionError(rule_name)
    return pred


def formula_string(rule_name: str, patt: str) -> str:
    formulas = {
        "g1": {"RRS": "9*fixoff", "RSR": "9*fixoff", "SRR": "162", "DIST": "0"},
        "g2": {"RRS": "6*fix", "RSR": "6*fix", "SRR": "6*fix", "DIST": "162"},
        "g3": {"RRS": "72-2*fix", "RSR": "54", "SRR": "54", "DIST": "54"},
        "g4": {"RRS": "162", "RSR": "6*fix", "SRR": "162", "DIST": "162"},
        "g5": {"RRS": "18*fixdiag", "RSR": "18*fixdiag", "SRR": "162", "DIST": "162"},
        "g6": {"RRS": "6*fix", "RSR": "81-3*fix", "SRR": "6*fix", "DIST": "54"},
    }
    return formulas[rule_name][patt]


def write_csv(path: Path, fieldnames: List[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> None:
    diagnostics: List[str] = []
    diagonals = list(itertools.product(S, repeat=3))

    # Main mixed-pattern formula table.
    formula_rows = []
    for rn in RULES:
        row = {
            "rule": rn,
            "rule_name": RULE_NAMES[rn],
            "RRR_formula": "RRR(pdist,ds) (see A3_rrr_formula_table.csv)",
            "RRS_formula": formula_string(rn, "RRS"),
            "RSR_formula": formula_string(rn, "RSR"),
            "SRR_formula": formula_string(rn, "SRR"),
            "DIST_formula": formula_string(rn, "DIST"),
            "scope": "Ω_elem only: six distinguished elementary cross-rules × 27 σ-equivariant diagonals",
        }
        formula_rows.append(row)
    write_csv(
        OUT / "A3_mixed_pattern_formula_table.csv",
        ["rule", "rule_name", "RRR_formula", "RRS_formula", "RSR_formula", "SRR_formula", "DIST_formula", "scope"],
        formula_rows,
    )

    # RRR formula table.
    rrr_rows = [
        {"pdist": 1, "ds_condition": "any", "RRR": 57, "fiber_associativity": "19/27", "notes": "constant base diagonal"},
        {"pdist": 2, "ds_condition": "0", "RRR": 33, "fiber_associativity": "11/27", "notes": "two values, no base fixed position"},
        {"pdist": 2, "ds_condition": ">=1", "RRR": 39, "fiber_associativity": "13/27", "notes": "two values, at least one base fixed position"},
        {"pdist": 3, "ds_condition": "any possible", "RRR": 27, "fiber_associativity": "9/27", "notes": "base diagonal is a permutation of S"},
    ]
    write_csv(
        OUT / "A3_rrr_formula_table.csv",
        ["pdist", "ds_condition", "RRR", "fiber_associativity", "notes"],
        rrr_rows,
    )

    # Full validation over Ω_elem.
    validation_rows = []
    mismatches = []
    spectrum_counter = Counter()
    by_rule_summary = defaultdict(lambda: Counter())

    for rn, gf in RULES.items():
        for d in diagonals:
            actual = actual_pattern_counts(gf, d)
            pred = predicted_counts(rn, d)
            total_actual = sum(actual[p] for p in PATTERNS)
            total_pred = sum(pred[p] for p in PATTERNS)
            params = diag_params(d)
            row = {
                "rule": rn,
                "rule_name": RULE_NAMES[rn],
                **params,
                "RRR_actual": actual["RRR"],
                "SRR_actual": actual["SRR"],
                "RRS_actual": actual["RRS"],
                "RSR_actual": actual["RSR"],
                "DIST_actual": actual["DIST"],
                "Assoc_actual": total_actual,
                "RRR_predicted": pred["RRR"],
                "SRR_predicted": pred["SRR"],
                "RRS_predicted": pred["RRS"],
                "RSR_predicted": pred["RSR"],
                "DIST_predicted": pred["DIST"],
                "Assoc_predicted": total_pred,
                "matches": "yes" if all(actual[p] == pred[p] for p in PATTERNS) else "no",
            }
            validation_rows.append(row)
            if row["matches"] != "yes":
                mismatches.append(row)
            spectrum_counter[total_actual] += 1
            by_rule_summary[rn][total_actual] += 1

    validation_fields = [
        "rule", "rule_name", "d", "d0", "d1", "d2", "pdist", "ds", "fix", "fixdiag", "fixoff",
        "RRR_actual", "SRR_actual", "RRS_actual", "RSR_actual", "DIST_actual", "Assoc_actual",
        "RRR_predicted", "SRR_predicted", "RRS_predicted", "RSR_predicted", "DIST_predicted", "Assoc_predicted", "matches",
    ]
    write_csv(OUT / "A3_mixed_pattern_validation.csv", validation_fields, validation_rows)

    # Summary rows by rule.
    summary_rows = []
    for rn in RULES:
        rule_rows = [r for r in validation_rows if r["rule"] == rn]
        for patt in PATTERNS:
            vals = sorted({int(r[f"{patt}_actual"]) for r in rule_rows})
            summary_rows.append({
                "rule": rn,
                "rule_name": RULE_NAMES[rn],
                "pattern": patt,
                "distinct_actual_values": " ".join(map(str, vals)),
                "min": min(vals),
                "max": max(vals),
                "formula": "RRR(pdist,ds)" if patt == "RRR" else formula_string(rn, patt),
            })
        totals = sorted({int(r["Assoc_actual"]) for r in rule_rows})
        summary_rows.append({
            "rule": rn,
            "rule_name": RULE_NAMES[rn],
            "pattern": "Assoc",
            "distinct_actual_values": " ".join(map(str, totals)),
            "min": min(totals),
            "max": max(totals),
            "formula": "sum of pattern formulas",
        })
    write_csv(
        OUT / "A3_mixed_pattern_summary.csv",
        ["rule", "rule_name", "pattern", "distinct_actual_values", "min", "max", "formula"],
        summary_rows,
    )

    diagnostics.append("A3 mixed-pattern formula verification")
    diagnostics.append(f"Ω_elem magmas checked: {len(validation_rows)}")
    diagnostics.append(f"Rules checked: {len(RULES)}")
    diagnostics.append(f"Diagonals checked: {len(diagonals)}")
    diagnostics.append(f"Row-pattern triples per magma: RRR=81 SRR=162 RRS=162 RSR=162 DIST=162 total=729")
    diagnostics.append(f"Formula mismatches: {len(mismatches)}")
    diagnostics.append(f"Associativity spectrum size from validation table: {len(spectrum_counter)}")
    diagnostics.append("Spectrum values: " + " ".join(map(str, sorted(spectrum_counter))))
    diagnostics.append("Formula table: A3_mixed_pattern_formula_table.csv")
    diagnostics.append("RRR table: A3_rrr_formula_table.csv")
    diagnostics.append("Validation table: A3_mixed_pattern_validation.csv")
    diagnostics.append("Summary table: A3_mixed_pattern_summary.csv")

    if mismatches:
        diagnostics.append("ERROR: mismatches detected")
        for mm in mismatches[:5]:
            diagnostics.append(str(mm))
        raise SystemExit("Formula mismatches detected; see output.")

    (OUT / "A3_mixed_pattern_count_output.txt").write_text("\n".join(diagnostics) + "\n", encoding="utf-8")
    print("\n".join(diagnostics))


if __name__ == "__main__":
    main()
