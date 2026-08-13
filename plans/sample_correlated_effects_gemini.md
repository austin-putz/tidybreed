# Plans: Sampling Correlated Stochastic Effects Across Simulated Time

**Status**: Draft for Review (Gemini Proposal)  
**Author**: Gemini CLI

---

## 1. The Core Challenge

In `tidybreed`, breeding simulations are highly dynamic. Phenotypes are often recorded at different points in simulated time (e.g., birth, weaning, test-entry, test-exit, slaughter), and culling or selection occurs between these events. 

If a user defines a residual covariance or correlated random effects between multiple phenotypes, those covariances must be honored even if the phenotypes are recorded:
* Days, weeks, or months apart.
* In separate `add_phenotype()` R calls.
* On different but overlapping subsets of individuals (due to culling).

Currently, `tidybreed` only draws correlated residuals/effects when they are recorded in the **same R call** on **identical subsets of individuals**. In any other scenario, it silently falls back to independent draws, completely ignoring the specified off-diagonal covariances.

---

## 2. Defects in the Current Implementation

1. **Residuals are stateless:** The helper `sample_residuals()` draws a joint matrix in memory. This is immediately consumed to compute liability/phenotype values and then discarded. Nothing about the realized residual is stored, making across-time conditioning mathematically impossible.
2. **Strict same-call/identical-subset gating:** Joint MVN residual draws are gated by a check requiring all requested phenotypes to have byte-identical, sorted lists of `id_ind`. If an animal is culled in between, the subsets differ, and the joint path is bypassed in favor of independent draws.
3. **Correlated named random effects lose correlation on partial re-draw:** While random effects are stored in `trait_random_effects`, drawing a new phenotype's correlated random effect for an existing level (e.g., pen `P1` has a stored draw for `ADG` but not `BW`) generates a fresh joint draw, writes the `BW` component, and discards the corresponding joint `ADG` component in favor of the previously stored, independent `ADG` value. This results in zero correlation on disk.

---

## 3. The Mathematical Solution: Conditional MVN

Rather than trying to eagerly pre-draw everything at birth (which fails for new animals or complex branching), we can use standard **conditional multivariate normal (MVN) conditioning**.

When drawing deviations for a set of target phenotypes $n$ (new), we check if any other phenotypes in the same correlated block $o$ (observed) already have stored deviations. We partition the full covariance block $R$:

$$
R = \begin{bmatrix}
R_{nn} & R_{no} \\
R_{on} & Roo
\end{bmatrix}
$$

The conditional distribution of the new deviations $e_n$ given the observed deviations $e_o$ is:

$$
e_n \mid e_o \sim \mathcal{N}\left( R_{no} R_{oo}^{-1} e_o, \,\, R_{nn} - R_{no} R_{oo}^{-1} R_{on} \right)
$$

For individuals with no prior observations in the block, $e_o$ is empty, the conditional mean is $0$, and the conditional variance is the full $R_{nn}$ (falling back to a standard unconditional joint draw).

---

## 4. Architectural Storage Options

### Option A: `residual_value` Column on `ind_phenotype`

We alter the `ind_phenotype` table on demand (or during initialization) to add a `residual_value DOUBLE` column. Realized residuals are written directly alongside the phenotype records they produced.

* **Pros:**
  * **Zero new tables** and zero new user-facing concepts.
  * **Automatic deletion consistency:** Deleting a phenotype record via `remove_rows()` automatically deletes its residual. No orphaned records or manual syncing is required.
  * **Highly diagnostic:** Users can easily verify the simulation's behavior using standard queries:
    ```r
    get_table("ind_phenotype") |> collect() |> summarise(var(residual_value))
    ```
  * Matches existing table alteration patterns (e.g., adding `liability_value` and `cat_name` on demand).
* **Cons:**
  * Solves residuals only. Correlated named random effects in `trait_random_effects` still require a separate (but mathematically identical) conditioning logic fix on their existing table.
  * If a phenotype record does not exist (e.g., censored/skipped), no residual is stored to condition on.

### Option B: Unified `random_draw` Table

We treat residuals as just another random effect where the "level" is the individual ID, merging `trait_random_effects` and residuals into a single unified table:

```sql
CREATE TABLE random_draw (
  effect_name    VARCHAR NOT NULL,   -- 'residual', 'pen', 'litter', 'pe'
  phenotype_name VARCHAR NOT NULL,
  level          VARCHAR NOT NULL,   -- id_ind for residuals; pen_id for pen
  record_number  INTEGER NOT NULL DEFAULT 1,
  draw_value     DOUBLE,
  PRIMARY KEY (effect_name, phenotype_name, level, record_number)
)
```

* **Pros:**
  * Unifies all stochastic simulation deviations into a single elegant concept.
  * Resolves both residuals and named effects through a single code path.
* **Cons:**
  * `record_number` is meaningless for all named effects (always `1`).
  * Shadows `ind_phenotype` row-for-row for residuals, introducing deletion-sync vulnerabilities (calling `remove_rows("ind_phenotype")` would leave orphaned rows in `random_draw` unless carefully synced).
  * `level VARCHAR` erases the explicit `id_ind` foreign-key typing.

### Option C: Dedicated `ind_residual` Table

We keep `ind_phenotype` strictly clean by storing residuals in a dedicated table:

```sql
CREATE TABLE ind_residual (
  id_ind          VARCHAR NOT NULL,
  phenotype_name  VARCHAR NOT NULL,
  pheno_number    INTEGER NOT NULL DEFAULT 1,
  residual_value  DOUBLE,
  PRIMARY KEY (id_ind, phenotype_name, pheno_number)
)
```

* **Pros:**
  * Keeps `ind_phenotype` clean.
  * Allows residuals to exist independently of a phenotype record (e.g., for pre-drawn lifetime trajectories at birth).
* **Cons:**
  * Adds an extra table to maintain.
  * Requires manual syncing on row deletion.
  * Over-engineers for a feature (pre-drawn trajectories) that is not currently planned.

---

## 5. Comparative Matrix

| Feature | Option A (`ind_phenotype` col) | Option B (`random_draw`) | Option C (`ind_residual` table) |
|---|---|---|---|
| **New Tables** | 0 | 1 (replaces 1) | 1 |
| **Fixes Residuals** | Yes | Yes | Yes |
| **Fixes Named Effects** | Yes (via shared helper) | Yes (same path) | Yes (via shared helper) |
| **Repeated Records** | Free (`pheno_number`) | Needs `record_number` | Free (`pheno_number`) |
| **Deletion Consistency** | Free / Automatic | Manual | Manual |
| **Typed `id_ind` FK** | Yes | No (`level VARCHAR`) | Yes |
| **Diagnostic Access** | High (same table) | Medium (separate table) | Medium (separate table) |
| **New Concepts** | 0 | 0 | 1 table |

---

## 6. Recommendation

**We strongly recommend Option A, combined with a conditioning upgrade to the existing `trait_random_effects` path.**

This approach is mathematically robust, maintains perfect referential integrity, and introduces zero schema bloat.

### Recommended Decisions:
1. **Repeated records are independent:** For a repeatable phenotype, records are independent across time. Across-time correlation of the *same* phenotype is the job of permanent-environment (PE) random effects (`define_effect_random(..., source_column = "id_ind")`), which already exists and persists. Therefore, we match residuals on `(id_ind, phenotype_name)` at the same `pheno_number` across *distinct* phenotypes.
2. **Always write residuals:** Populating `residual_value` for every record (even when uncorrelated) is extremely clean and provides massive diagnostic value for simulation builders checking their models.
3. **Zero new user concepts:** The user simply calls `define_residual_cov()` and `add_phenotype()` exactly as they do today. The system seamlessly handles the covariance behind the scenes, regardless of simulated time gaps or selection/culling.

---

## 7. Implementation Sketch

### Step 1: Add internal helper in `R/phenotype_helpers.R`
We write a robust conditional drawer:

```r
#' Draw new deviations conditional on already-realized ones
#'
#' @param R        Full covariance matrix over the correlated block, dimnamed.
#' @param new      Character vector of names to draw.
#' @param observed Numeric matrix (n x length(obs_names)) of realized deviations,
#'                 column-dimnamed; NULL or zero-column for an unconditional draw.
#' @return n x length(new) matrix, column-dimnamed by `new`.
draw_conditional_mvn <- function(R, new, observed = NULL) {
  n <- if (!is.null(observed)) nrow(observed) else 1L
  if (is.null(observed) || ncol(observed) == 0L) {
    return(MASS::mvrnorm(n, mu = rep(0, length(new)), Sigma = R[new, new, drop = FALSE]))
  }
  
  o    <- colnames(observed)
  Roo  <- R[o, o, drop = FALSE]
  Rno  <- R[new, o, drop = FALSE]
  
  # Solve Roo %*% B_t = Rno_t -> B = Rno %*% solve(Roo)
  # Using robust solve or Cholesky decomposition
  B    <- Rno %*% solve(Roo)
  Sc   <- R[new, new, drop = FALSE] - B %*% t(Rno) # Conditional covariance
  
  # Ensure conditional covariance remains positive semi-definite
  # (e.g., handle tiny numerical negatives on diagonal)
  diag(Sc)[diag(Sc) < 0] <- 0
  
  mu_c <- observed %*% t(B) # Conditional mean: n x length(new)
  
  draws <- mu_c + MASS::mvrnorm(n, rep(0, length(new)), Sc)
  colnames(draws) <- new
  draws
}
```

### Step 2: Update `add_phenotype()` Residual Path
We group the individuals in the current call by their **observation pattern** (which correlated residuals they already have stored on disk). For each group, we construct the `observed` matrix, call `draw_conditional_mvn()`, and store the newly drawn `residual_value` along with the record in `ind_phenotype`.

### Step 3: Update `add_phenotype()` Random Effects Path
We apply the exact same `draw_conditional_mvn()` helper to the named random effects drawn at Step 7.5. Instead of overwriting and discarding joint draws, we condition on any previously written levels in `trait_random_effects`.

---

## 8. Verification & Testing Strategy

Since sequential draws are order-dependent, a sequential draw path and a joint draw path will not produce identical numeric results under the same seed (as the number of RNG calls differs). Our test suite must focus on:

1. **Exact within-path reproducibility:** Same seed and same sequential call sequence must yield byte-identical results.
2. **Distributional parity:** At large $N$ (e.g., 5,000 animals), sequential draws of `on_test_wt` followed by `off_test_wt` (with culling) must yield a realized residual correlation matching the target covariance matrix within a small statistical tolerance.
3. **Culling edge case:** Phenotype A on 5,000 animals, Phenotype B on a 3,000-animal subset. The realized correlation on the 3,000 overlapping animals must match the target, while the 2,000 culled animals remain untouched.
4. **Random effect parity:** Stored level draws in `trait_random_effects` are correctly conditioned upon, and the realized covariance across many levels matches the target covariance matrix.
