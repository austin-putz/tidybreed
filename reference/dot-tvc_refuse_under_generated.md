# Refuse a genetic target under terms a generator calibrated (Q21)

A `"generated"` term is calibrated to the target it was generated with,
and the prevalence threshold trusts the stored target for that reason.
Writing a target under such terms (even after removing the old one)
would break that. Only the exported
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
calls this; a generator's `G =` writes through
[`.tvc_write_block()`](https://austin-putz.github.io/tidybreed/reference/dot-tvc_write_block.md)
together with the terms it calibrates, so it is never refused here.

## Usage

``` r
.tvc_refuse_under_generated(conn, effect_name, traits, line_name)
```

## Details

Scope as in
[`.tvc_generated_traits()`](https://austin-putz.github.io/tidybreed/reference/dot-tvc_generated_traits.md).
