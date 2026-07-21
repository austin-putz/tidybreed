# Upsert TBV rows — column-named INSERT with ON CONFLICT DO UPDATE SET

Only the columns present in tbv_df are updated on conflict; all other
columns (user-defined extras) are left untouched in existing rows.

## Usage

``` r
upsert_ind_tbv(pop, tbv_df)
```
