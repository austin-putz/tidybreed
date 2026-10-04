# Value components of `ind_tgv`

The closed set of `ind_tgv.component_name` values
(plans/import_qtl_effect_methods.md §6B): the `contrast_name` of a
one-locus term, or `interaction` for a term over two or more loci.
Validated in R, not by an SQL `CHECK`, so a name can be added in one
line. `"total"` is not a row; it is the `ind_tgv_total` view.

## Usage

``` r
TGV_COMPONENT_NAMES
```
