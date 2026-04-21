# σ-Equivariant Elementary Magmas on AG(2, 3)

> **σ-Equivariant Elementary Magmas on AG(2, 3): Classification, associativity spectrum, and canonical selection in a six-rule family**

## Main verified results

- The full space of σ-equivariant cross-row rules has cardinality \(3^{18} = 387,420,489\).
- The elementary family has 162 objects and decomposes into 90 diagonal \(S_3\)-orbits.
- The associativity spectrum inside the elementary family has exactly 28 values, ranging from 162 to 567.
- Inside the restricted extended space
  \[
  \Omega'_{\mathrm{elem}}=\mathcal F_\sigma\times\{g_1,\dots,g_6\}\times\Delta_\sigma,
  \]
  the verified sequential chain is
  \[
  1458 \to 486 \to 54 \to 1,
  \]
  with PAB as the unique final survivor.

## Repository layout

```text
.
├── README.md
├── paper/
│   ├── main.tex
│   ├── references.bib
│   └── main.pdf
└── proof_artifacts/
    ├── README.md
    ├── A1_cross_rule_count/
    ├── A2_burnside/
    ├── A3_mixed_pattern_formula/
    ├── A4_spectrum/
    ├── A5_omega_prime_selection/
    └── P1_structural_polish/
```

## Directory guide

### `paper/`

Current paper source and compiled PDF.

- `main.tex` — current preprint source.
- `references.bib` — bibliography.
- `main.pdf` — compiled preprint.

### `proof_artifacts/`

Finite verification scripts together with their generated CSV/TXT outputs.

- `A1_cross_rule_count/` — cross-rule orbit decomposition and the 318318 count.
- `A2_burnside/` — Burnside fixed-point table and the 90-orbit count.
- `A3_mixed_pattern_formula/` — row-pattern formula tables and full validation over the elementary family.
- `A4_spectrum/` — spectrum values, witnesses, and full validation.
- `A5_omega_prime_selection/` — restricted sequential selection on the elementary family.
- `P1_structural_polish/` — automorphisms, endomorphisms, standard selection table, symbolic comparisons, and absorption graph checks.

## Reproducibility

The repository includes both the verification scripts and the generated finite artifacts used in the current release. The paper source is in `paper/main.tex`, and `paper/main.pdf` is the corresponding compiled version.