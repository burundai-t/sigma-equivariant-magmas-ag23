# Proof Artifacts

This directory contains the finite verification layer supporting the results stated in the main paper.

Each subdirectory corresponds to a specific theorem block and contains:
- a reproducible script,
- generated CSV/TXT outputs,
- (optionally) a short report describing the verification.

All computations are performed by exhaustive enumeration over finite state spaces.

---

## Structure

### A1 — Cross-rule count

Verifies the size of the σ-equivariant cross-rule space.

- Result:
  \[
  |\mathrm{CrossRules}_\sigma| = 3^{18}.
  \]

- Method:
  decomposition of 54 cross-row inputs into 18 σ-orbits.

---

### A2 — Burnside orbit count

Computes the number of diagonal \(S_3\)-orbits in the elementary family.

- Result:
  \[
  |\Omega_{\mathrm{elem}} / S_3| = 90.
  \]

- Method:
  Burnside’s lemma with explicit fixed-point enumeration.

---

### A3 — Associativity master formula

Validates the row-pattern decomposition formula for associativity.

- Scope:
  all \(162\) elementary magmas and all \(729\) triples per magma.

- Result:
  exact agreement between formula and enumeration (zero mismatches).

---

### A4 — Associativity spectrum

Computes the full associativity spectrum in the elementary family.

- Result:
  - exactly 28 distinct values,
  - range: \(162 \leq \mathrm{Assoc} \leq 567\),
  - all values divisible by 3.

---

### A5 — Restricted sequential selection

Verifies the sequential selection process in the extended elementary space

\[
\Omega'_{\mathrm{elem}} = \mathcal F_\sigma \times \{g_1,\dots,g_6\} \times \Delta_\sigma.
\]

- Result:
  \[
  1458 \to 486 \to 54 \to 1,
  \]
  with PAB as the unique final survivor.

- Important:
  this is a **sequential** selection theorem, not a global minimisation over the full space.

---

### P1 — Structural polishing

Provides finite verification tables for structural properties of the PAB magma.

Includes:
- automorphism group (\(\cong S_3\)),
- endomorphism monoid (6 automorphisms + 3 constants),
- standard canonical selection table,
- symbolic \(F_\beta\)-comparisons,
- absorption graph verification.

---

## Reproducibility notes

- All scripts are deterministic and use exact arithmetic.
- No floating-point comparisons are used for algebraic claims.
- Generated CSV files serve as stable verification artifacts.
- Re-running the scripts should reproduce all outputs exactly.

---

## Scope

All results in this directory concern one of the following spaces:

- \(\Omega_{\mathrm{elem}}\) — elementary family (162 objects),
- \(\Omega'_{\mathrm{elem}}\) — extended elementary family (1458 objects).

They do **not** claim results for the full cross-rule universe
\[
\Omega_{\mathrm{full}}.
\]