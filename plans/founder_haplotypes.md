# Founder Haplotype Generation — Design Plan

## Core Principle

All founder haplotype methods write the same output:
- `founder_haplotypes` table in DuckDB (`hap_id`, `line_name`, `locus_1..locus_n`)
- `founder_allele_freq` column in `genome_meta`

The source of those haplotypes varies enormously, so functions are split by **source**, not by method.

---

## Status Summary

| Function | Status |
|---|---|
| `define_founder_haplotypes()` | ✅ Complete |
| `simulate_founder_haplotypes()` | ⬜ Not started |
| `import_founder_haplotypes()` | ⬜ Not started |

---

## 1. `define_founder_haplotypes()` — parametric, no external dependencies ✅ IMPLEMENTED

**File**: `R/define_founder_haplotypes.R`, `R/founder_haplotype_helpers.R`

Generates haplotypes internally — no external software needed.

### Implemented methods

| Method | Description | LD structure |
|---|---|---|
| `"fixed"` | Every locus gets the same `allele_freq` (default 0.5) | None |
| `"uniform"` | Per-locus frequency from Uniform(`min_allele_freq`, `max_allele_freq`) | None |
| `"beta"` | Per-locus frequency from Beta(`beta_shape1`, `beta_shape2`); default shape1=shape2=0.5 (Jeffreys prior) | None |
| `"balding_nichols"` | Beta parameterized by `fst` and `mean_freq`; models drift around a mean | None |
| `"mosaic"` | Builds haplotypes by copying from `n_templates` templates with Poisson switches at rate `switch_rate` per Mb; requires `chr`/`pos_Mb` in `genome_meta` | Yes (block LD) |
| `"gaussian_copula"` | AR(1) latent Gaussian thresholded at per-locus frequencies; LD decays as ρ = exp(−`decay_rate` × d_Mb); requires `chr`/`pos_Mb` in `genome_meta` | Yes (exponential decay) |

### Key design points

- Multi-line support via `line_name` — rows tagged, `hap_id` prefixed (`"LineA_hap_1"`)
- Wrong-method argument validation: passing `fst` to `method = "uniform"` errors cleanly
- `founder_allele_freq` written to `genome_meta` after every call
- Shared internal helper `.write_founder_haplotypes()` handles all DB writes

---

## 2. `simulate_founder_haplotypes()` — coalescent software ⬜ NOT STARTED

Shells out to a coalescent simulator, parses stdout, loads haplotypes into the DB. Produces realistic LD structure from population genetics models. **Priority: implement MaCS first.**

### Target software (all output MS format)

| Software | Notes | Priority |
|---|---|---|
| `MaCS` | Efficient for large genomes; low memory footprint | **First** |
| `scrm` | ms-compatible, fast, widely available | Second |
| `ms` (Hudson) | Original; reference implementation | Third |
| `msms` | ms + selection | Later |
| `fastsimcoal2` | Complex demographic models; different output format | Later |

SLiM is different enough (script-based, forward-time, completely different I/O) to warrant a separate function if supported — not part of this design.

### Proposed signature

```r
simulate_founder_haplotypes(pop, software = "macs", n_haplotypes, ...)
```

### Simulation behaviors

- Reads chromosome structure and positions from `genome_meta` (position in Mb → recombination distance)
- Constructs the simulator command string internally from `genome_meta` + user parameters
- Checks that the requested binary is on `PATH` (`Sys.which(software)`) and errors clearly
- Shells out via `system2()`, captures stdout, parses MS format
- Writes result to `founder_haplotypes` table via `.write_founder_haplotypes()`
- Supports `line_name` argument (same semantics as `define_founder_haplotypes()`)

### MaCS-specific notes

- MaCS command takes recombination rate (4Ner) and mutation rate (4Neu) per bp; derive from `genome_meta` positions
- Output: MS format (one header block, then `//` separated haplotype matrices per chromosome segment)
- Each locus position in MS output must be mapped back to `genome_meta` locus names by position

---

## 3. `import_founder_haplotypes()` — user-provided real data ⬜ NOT STARTED

Loads existing haplotypes from a reference panel or real population. The path for **realistic simulations** grounded in observed genetic variation.

### Proposed `from` argument values

| Value | Description | Dependency |
|---|---|---|
| `"matrix"` | User passes an R matrix (rows = haplotypes, cols = loci, values = 0/1) | None |
| `"vcf"` | Read a phased VCF file | `vcfR` or `VariantAnnotation` |
| `"plink"` | Read PLINK `.bed/.bim/.fam` | `BEDMatrix` or `genio` |

### Import function signature

```r
import_founder_haplotypes(pop, from = "matrix", haplotype_matrix = NULL,
                          file = NULL, line_name = NULL, ...)
```

### Import behaviors

- **Locus alignment** is the central challenge: imported loci must match `genome_meta` exactly
  (same locus names and/or positions). Function validates alignment and errors clearly
  when loci are missing or mismatched.
- Allele frequencies are derived from the imported data (column means), not sampled
- Writes result via `.write_founder_haplotypes()` — same output schema as other functions

### From-argument implementation order

1. `from = "matrix"` — no file I/O, no new dependencies; validates column names match `locus_name` in `genome_meta`
2. `from = "vcf"` — requires phased VCF; position-based locus matching
3. `from = "plink"` — `.bim` file provides locus IDs for matching

---

## Shared Internal Helper (implemented)

```r
.write_founder_haplotypes(pop, haplotype_matrix, allele_freqs, line_name)
```

Handles: table creation if absent, `line_name` column guard (old-format detection),
column naming (`locus_1 … locus_n`), `hap_id` prefixing, `pop$tables` registration,
`founder_allele_freq` write-back to `genome_meta`.

---

## Implementation Order (remaining work)

1. `simulate_founder_haplotypes()` with `software = "macs"` — **next priority**
2. `import_founder_haplotypes()` with `from = "matrix"` — no dependencies, high user value
3. `import_founder_haplotypes()` with `from = "vcf"` and `from = "plink"`
4. Additional coalescent simulators (`scrm`, `ms`) in `simulate_founder_haplotypes()`
