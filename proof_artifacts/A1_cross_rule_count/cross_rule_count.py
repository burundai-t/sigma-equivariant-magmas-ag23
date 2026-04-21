#!/usr/bin/env python3
"""
A1 proof artifact: count sigma-equivariant cross-rules for CTA/PAB.

A cross-rule is a function
    g(r1, c1, r2, c2) -> c_out in S
for cross-row inputs r1 != r2, where S = Z/3Z.

Sigma-equivariance means
    g(r1+1, c1+1, r2+1, c2+1) = g(r1,c1,r2,c2)+1 mod 3.

The 54 cross-row inputs split into 18 free Z/3Z-orbits; each orbit has one
free output value in S. Therefore the full space has 3^18 cross-rules.

Running this script writes:
    cross_rule_orbits.csv
    elementary_rule_vectors.csv
and prints a compact verification summary.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, List, Sequence, Tuple

S = (0, 1, 2)
Input = Tuple[int, int, int, int]
Rule = Callable[[int, int, int, int], int]


def sigma_value(x: int, k: int = 1) -> int:
    return (x + k) % 3


def sigma_input(t: Input, k: int = 1) -> Input:
    return tuple(sigma_value(v, k) for v in t)  # type: ignore[return-value]


def complement(a: int, b: int) -> int:
    """Return the unique element of S minus {a,b}; requires a != b."""
    if a == b:
        raise ValueError(f"complement undefined for equal inputs: {a}, {b}")
    return ({0, 1, 2} - {a, b}).pop()


def all_cross_inputs() -> List[Input]:
    return [(r1, c1, r2, c2) for r1 in S for c1 in S for r2 in S for c2 in S if r1 != r2]


def orbit_of(t: Input) -> Tuple[Input, Input, Input]:
    return (t, sigma_input(t, 1), sigma_input(t, 2))


def canonical_orbit(t: Input) -> Tuple[Input, Input, Input]:
    """Return orbit ordered from its lexicographically minimal representative."""
    reps = orbit_of(t)
    rep = min(reps)
    return orbit_of(rep)


def compute_orbits() -> List[Tuple[Input, Input, Input]]:
    seen = set()
    orbits: List[Tuple[Input, Input, Input]] = []
    for t in all_cross_inputs():
        if t in seen:
            continue
        orb = canonical_orbit(t)
        for u in orb:
            seen.add(u)
        orbits.append(orb)
    orbits.sort(key=lambda o: o[0])
    return orbits


# The six distinguished elementary rules used in Omega_elem.
def g1(r1: int, c1: int, r2: int, c2: int) -> int:
    return r2


def g2(r1: int, c1: int, r2: int, c2: int) -> int:
    return c2


def g3(r1: int, c1: int, r2: int, c2: int) -> int:
    return complement(c1, c2) if c1 != c2 else c1


def g4(r1: int, c1: int, r2: int, c2: int) -> int:
    return c1


def g5(r1: int, c1: int, r2: int, c2: int) -> int:
    return r1


def g6(r1: int, c1: int, r2: int, c2: int) -> int:
    return complement(r2, c2) if r2 != c2 else r2


ELEMENTARY_RULES: Sequence[Tuple[str, Rule, str]] = (
    ("g1", g1, "r2 / PAB column-blind"),
    ("g2", g2, "c2 / transparent"),
    ("g3", g3, "c1c2 if c1 != c2 else c1"),
    ("g4", g4, "c1 / echo"),
    ("g5", g5, "r1 / self-referential"),
    ("g6", g6, "r2c2 if r2 != c2 else r2"),
)


def is_equivariant_rule(rule: Rule) -> bool:
    for t in all_cross_inputs():
        lhs = rule(*sigma_input(t, 1))
        rhs = sigma_value(rule(*t), 1)
        if lhs != rhs:
            return False
    return True


def tuple_str(t: Input) -> str:
    return f"({t[0]},{t[1]},{t[2]},{t[3]})"


def orbit_str(orb: Sequence[Input]) -> str:
    return " -> ".join(tuple_str(t) for t in orb)


def equivariant_values_for_seed(seed: int) -> str:
    return f"{seed},{(seed+1)%3},{(seed+2)%3}"


def rule_vector(rule: Rule, orbits: Sequence[Sequence[Input]]) -> List[int]:
    """Values on canonical orbit representatives; this is the seed vector in S^18."""
    return [rule(*orb[0]) for orb in orbits]


def write_orbits_csv(path: Path, orbits: Sequence[Sequence[Input]]) -> None:
    headers = [
        "orbit_id",
        "representative",
        "rep_r1",
        "rep_c1",
        "rep_r2",
        "rep_c2",
        "orbit_inputs_sigma_order",
        "orbit_size",
        "seed_0_outputs_on_orbit",
        "seed_1_outputs_on_orbit",
        "seed_2_outputs_on_orbit",
    ]
    for name, _rule, desc in ELEMENTARY_RULES:
        headers.append(f"{name}_rep_value")
    for name, _rule, desc in ELEMENTARY_RULES:
        headers.append(f"{name}_orbit_outputs")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for i, orb in enumerate(orbits, start=1):
            rep = orb[0]
            row = {
                "orbit_id": i,
                "representative": tuple_str(rep),
                "rep_r1": rep[0],
                "rep_c1": rep[1],
                "rep_r2": rep[2],
                "rep_c2": rep[3],
                "orbit_inputs_sigma_order": orbit_str(orb),
                "orbit_size": len(orb),
                "seed_0_outputs_on_orbit": equivariant_values_for_seed(0),
                "seed_1_outputs_on_orbit": equivariant_values_for_seed(1),
                "seed_2_outputs_on_orbit": equivariant_values_for_seed(2),
            }
            for name, rule, _desc in ELEMENTARY_RULES:
                row[f"{name}_rep_value"] = rule(*rep)
            for name, rule, _desc in ELEMENTARY_RULES:
                row[f"{name}_orbit_outputs"] = ",".join(str(rule(*t)) for t in orb)
            writer.writerow(row)


def write_rule_vectors_csv(path: Path, orbits: Sequence[Sequence[Input]]) -> None:
    headers = ["rule", "description", "is_sigma_equivariant"] + [f"orbit_{i}" for i in range(1, len(orbits) + 1)]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for name, rule, desc in ELEMENTARY_RULES:
            vec = rule_vector(rule, orbits)
            row = {"rule": name, "description": desc, "is_sigma_equivariant": is_equivariant_rule(rule)}
            row.update({f"orbit_{i}": vec[i - 1] for i in range(1, len(orbits) + 1)})
            writer.writerow(row)


def verify_orbits(orbits: Sequence[Sequence[Input]]) -> None:
    inputs = all_cross_inputs()
    flat = [t for orb in orbits for t in orb]
    if len(inputs) != 54:
        raise AssertionError(f"Expected 54 cross-row inputs, got {len(inputs)}")
    if len(orbits) != 18:
        raise AssertionError(f"Expected 18 orbits, got {len(orbits)}")
    if len(flat) != len(set(flat)):
        raise AssertionError("Orbit decomposition contains duplicate inputs")
    if set(flat) != set(inputs):
        raise AssertionError("Orbit decomposition does not cover all cross-row inputs")
    if any(len(orb) != 3 for orb in orbits):
        raise AssertionError("A sigma-orbit has size different from 3")
    for name, rule, _desc in ELEMENTARY_RULES:
        if not is_equivariant_rule(rule):
            raise AssertionError(f"Elementary rule {name} is not sigma-equivariant")


def main() -> int:
    parser = argparse.ArgumentParser(description="Count sigma-equivariant cross-rules for CTA/PAB.")
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent,
                        help="Directory for generated CSV artifacts. Default: script directory.")
    args = parser.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    orbits = compute_orbits()
    verify_orbits(orbits)

    orbit_count = len(orbits)
    cross_rule_count = 3 ** orbit_count

    orbits_csv = outdir / "cross_rule_orbits.csv"
    vectors_csv = outdir / "elementary_rule_vectors.csv"
    write_orbits_csv(orbits_csv, orbits)
    write_rule_vectors_csv(vectors_csv, orbits)

    print("A1 sigma-equivariant cross-rule count")
    print(f"S size: 3")
    print(f"Cross-row input count: {len(all_cross_inputs())}")
    print(f"Sigma-orbit count on inputs: {orbit_count}")
    print(f"Orbit sizes: {sorted({len(o) for o in orbits})}")
    print(f"Cross-rule count: 3^{orbit_count} = {cross_rule_count}")
    print("Elementary rules sigma-equivariant:")
    for name, rule, desc in ELEMENTARY_RULES:
        print(f"  {name}: {is_equivariant_rule(rule)}  ({desc})")
    print(f"Wrote: {orbits_csv}")
    print(f"Wrote: {vectors_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
