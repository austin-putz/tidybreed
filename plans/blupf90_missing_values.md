# Plan: missing values in the BLUPF90 files

**Status: to do next, with the `add_ebv()` redesign. Not implemented.**

Found while fixing the empty `phenotype` filter in 0.76.15. Left alone on
purpose: the redesign will rework how `add_ebv()` writes its files, and missing
values need to be designed into that, not patched beforehand.

## Problem

`build_data_file()` (`R/blupf90_helpers.R`) writes one `data.txt` that mixes
three kinds of column, and they don't share one missing-value convention.

### Responses (trait columns)

- An animal with no value for a trait is written as `0` (the "Replace NA with
  0" step after the pivot), and `write_renum_par()` sets `OPTION missing 0`.
- So **an observed value of exactly 0 is also read as missing.** This affects:
  - binary traits (every 0 in a 0/1 trait such as mortality is dropped, so the
    trait is fitted on its 1s only);
  - count traits (a litter size of 0);
  - any continuous trait centered at or near 0, such as a residual-scale trait
    or a trait defined with `mean = 0`.
- Not yet confirmed against a BLUPF90 run; this is from reading the code.

### Fixed effects: no missing code exists

- The user found out the hard way that `OPTION missing` applies to **responses
  only**. BLUPF90 has no missing-value code for a fixed effect.
- **Fixed class:** a character `NA` is written as the string `"0"`
  (`col_vals[is.na(col_vals)] <- "0"`), so `0` becomes one more **level** of
  the effect. All animals missing the class are pooled into one group. This may
  be acceptable in some cases (an "unknown herd" level), but it is a modelling
  choice that is made silently today. A real level that happens to be called
  `"0"` would merge with the missing animals.
- **Fixed covariate:** a numeric `NA` is left as `NA`, and `write.table()`
  writes the literal text `NA` into `data.txt`. How renumf90 reads that is not
  known: it might error, read it as 0 (silently, with the covariate set to 0),
  or do something else. This needs checking.

### Animals whose records are all NULL

- An animal whose only records have a NULL `pheno_value` still gets a row: the
  `list_avg` mean is NULL, so the row is written as missing.
- `OPTION remove_all_missing` then drops rows where every response is missing,
  but this interacts with the response problem above: an animal whose real
  values are all 0 would be dropped as well.
- The per-trait count message added in 0.76.15 already excludes these animals,
  and every trait at 0 is an error, so an all-missing file never reaches the
  solver.

## Questions for the redesign

1. **Response sentinel.** Use a value that can't occur in real data (commonly
   `-999` or similar) and set `OPTION missing` to match. Should it be fixed, or
   chosen per evaluation from the observed range? Does BLUPF90 need it to be an
   integer?
2. **Missing fixed class.** For each effect, either drop the record, or give it
   an explicit, named level (`"__missing__"`, not `"0"`). Which should be the
   default, and should the user be able to choose per effect (a field on
   `phenotype_effects`)?
3. **Missing fixed covariate.** Dropping the record is probably the only safe
   option; imputing 0 or the mean is a modelling decision. Confirm what
   renumf90 does with `NA` today before choosing.
4. **Report what was dropped.** Extend the per-trait `Records [...]` message
   with the number of records lost to a missing covariate or class.
5. **Verify against BLUPF90 itself:**
   - a binary trait containing 0s;
   - a covariate with an `NA`;
   - a class with `NA` next to a real level named `"0"`.

   Check the BLUPF90 / renumf90 documentation for `OPTION missing`,
   `remove_all_missing` and how renumf90 handles non-numeric covariate fields.

## Tests to write then

- A binary trait with 0s keeps every record in `data.txt` (not coded as
  missing).
- A missing response and an observed 0 are written differently.
- A missing class level never shares a code with a real level.
- A missing covariate follows the chosen rule (dropped or an explicit level),
  never a silent 0 and never the text `NA`.
