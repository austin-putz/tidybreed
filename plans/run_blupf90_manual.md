Good news: the folder creation and BLUPF90 execution is already built

`add_ebv(software = "blupf90")` handles all of this automatically — it creates the timestamped eval folder, writes data/pedigree/genotype files, runs renumf90 + blupf90+, parses the solutions, and loads EBVs into ind_ebv. The eval_blup_20260619_212546/ you see in git status was created by a prior call.

Standard call (works for any trait already in trait_meta + phenotype_meta)

```r
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen >= 2L) |>
  add_ebv(
    trait_name     = c("ADG", "BF"),   # must exist in trait_meta
    software       = "blupf90",
    model          = "blup_v1",
    run_dir        = ".",              # eval_blup_v1_YYYYMMDD_HHMMSS/ created here
    n_gen_pedigree = 3L
  )
```

The folder will contain: data.txt, pedigree.txt, renum.par, meta.txt, renumf90.out/err, blupf90.out/err, solutions.orig.

---
For a complex trait the built-in write_renum_par() doesn't support yet

The current code generates a standard animal model (mu + fixed class/cov effects + random animal). If you need something more complex — maternal effects, permanent environment, multiple random effects, custom covariance structure — the workflow is:

Step 1: Let tidybreed write the data files, then customize renum.par

Run add_ebv() but have renumf90/blupf90+ fail (or not be on PATH temporarily) to get the data

# Or just use the internal helpers directly to dump files manually

Actually the simpler path: create the eval folder yourself, export the data, write your own rload results back.

Step 2: Export data from the database

```r
library(dplyr)

eval_dir <- "my_complex_eval"
dir.create(eval_dir, showWarnings = FALSE)

# --- Pedigree ---
ped <- get_table(pop, "ind_meta") |>
  collect() |>
  mutate(
    id_parent_1 = if_else(is.na(id_parent_1), "0", id_parent_1),
    id_parent_2 = if_else(is.na(id_parent_2), "0", id_parent_2)
  ) |>
  select(id_ind, id_parent_1, id_parent_2)

write.table(ped, file.path(eval_dir, "pedigree.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# --- Phenotype data ---
pheno <- get_table(pop, "ind_phenotype") |>
  filter(phenotype_name %in% c("MyTrait")) |>
  collect() |>
  select(id_ind, phenotype_name, pheno_value)

# Join covariates from ind_meta (sex, gen, etc.)
meta <- get_table(pop, "ind_meta") |> collect()
data_df <- pheno |>
  left_join(meta |> select(id_ind, sex, gen), by = "id_ind") |>
  mutate(mu = 1L)

write.table(data_df, file.path(eval_dir, "data.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)
```

Step 3: Write your own renum.par

Put the parameter file in `my_complex_eval/renum.par` with whatever model you need (maternal, P

Step 4: Run renumf90 + blupf90+ manually

```r
old_wd <- getwd()
setwd(eval_dir)
system2("renumf90", args = "renum.par", stdout = "renumf90.out", stderr = "renumf90.err")
system2("blupf90+", args = "renf90.par", stdout = "blupf90.out",  stderr = "blupf90.err")
setwd(old_wd)
```

Step 5: Load EBVs back into tidybreed

```r
# Parse solutions file (5 cols with OPTION origID: trait effect level orig_id solution)
sols_raw <- read.table(file.path(eval_dir, "solutions.orig"),
                       header = FALSE, stringsAsFactors = FALSE)
colnames(sols_raw) <- c("trait_num", "effect_num", "level", "id_ind", "ebv_value")

# Keep only animal effect rows (effect_num = your animal effect number)
animal_effect_num <- 3L   # adjust to match your model
ebv_df <- sols_raw |>
  filter(effect_num == animal_effect_num) |>
  mutate(
    trait_name  = "MyTrait",
    model       = "custom_maternal_v1",
    eval_number = 1L,
    acc         = NA_real_,
    se          = NA_real_
  ) |>
  select(id_ind, trait_name, model, ebv_value, acc, se, eval_number)

# Insert into ind_ebv via DBI
DBI::dbAppendTable(pop$db_conn, "ind_ebv",
  dplyr::mutate(ebv_df,
    id_ebv = seq.int(
      tidybreed:::next_int_id(pop$db_conn, "ind_ebv", "id_ebv"),
      length.out = nrow(ebv_df)
    )
  )
)
```

---
What trait complexity matters for choosing the path

┌───────────────────────────────────────────────────────────────────────────────────────────────────────────┬─────────────────────────────────────────────────────────────────────────────────┐
│                                                Complexity                                                            Use                                       │
├───────────────────────────────────────────────────────────────────────────────────────────────────────────┼─────────────────────────────────────────────────────────────────────────────────┤
│ Standard animal model (any single/multi-trait with fixed effects defined via define_effect_ware = "blupf90") — fully automated                                 │
├───────────────────────────────────────────────────────────────────────────────────────────────────────────┼─────────────────────────────────────────────────────────────────────────────────┤
│ Custom model structure (maternal, PE, dominance)                                           t → custom renum.par → manual run → DBI::dbAppendTable to load back │
└───────────────────────────────────────────────────────────────────────────────────────────────────────────┴─────────────────────────────────────────────────────────────────────────────────┘                      

The `meta.txt` file the automated path writes (`my_eval_dir/meta.txt`) is a good reference for the column layout of `data.txt` and the effect numbering — worth checking even if you're doing a custom run, since the data export logic is shared.






Yes — it's `.create_run_dir()`. Here's the full picture:

The folder system

`open_pop()` sets up a 4-layer folder structure. You declare tool names when opening the pop:

```r
pop <- open_pop(
  pop_name     = "swine_sim",
  output_dir   = "tidybreed_output",
  scenario_dir = "scenario_A",
  tools        = c("blupf90", "plink")   # <-- registers tool dirs at layer 4
)
```

This creates:
tidybreed_output/
  scenario_A/
    sim.duckdb
    blupf90/          # pop$run_dirs[["blupf90"]]
    plink/            # pop$run_dirs[["plink"]]

`.create_run_dir()` — the function you're thinking of

It creates a timestamped layer-5 subfolder inside any registered tool dir and returns the path:

```r
run_dir <- tidybreed:::.create_run_dir(pop, tool = "blupf90")
# e.g. tidybreed_output/scenario_A/blupf90/20260619_212546_x7m3qk/
```

The name format is `YYYYMMDD_HHMMSS_<6-char random suffix>` so every call gets a unique isolated folder.

Full manual workflow

# 1. Create the run folder

```r
run_dir <- tidybreed:::.create_run_dir(pop, tool = "blupf90")
```

# 2. Write your files there

```r
write.table(ped_df,   file.path(run_dir, "pedigree.txt"), ...)
write.table(data_df,  file.path(run_dir, "data.txt"), ...)
writeLines(par_lines, file.path(run_dir, "renum.par"))
```

# 3. Run BLUPF90 manually

```
old_wd <- getwd()
setwd(run_dir)
system2("renumf90", args = "renum.par", stdout = "renumf90.out")
system2("blupf90+", args = "renf90.par", stdout = "blupf90.out")
setwd(old_wd)
```

# 4. Load results back 

The key requirement: the tool name ("blupf90") must have been declared in tools = when you called open_pop(), otherwise .create_run_dir() errors listing what's registered. If your pop was opened without tools, you can check pop$run_dirs — if it's empty or missing "blupf90", you'd need to reopen with clean st create the folder manually with dir.create().



