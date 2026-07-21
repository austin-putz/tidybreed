# Draw a Poisson count from the current dqrng stream (log-accumulation inversion)

Knuth's multiplicative inversion in log space: accumulate `log(u)` and
count steps until the sum drops to `-lambda`. Consumes ONLY dqrng
uniforms (no base `rpois`, no `std::poisson_distribution`), so the R
reference and the future C++ kernel match by construction from the
shared uniform stream. Log form (not the raw product) avoids underflow
when `exp(-lambda)` is denormal.

## Usage

``` r
.draw_poisson_dqrng(lambda)
```

## Arguments

- lambda:

  Non-negative Poisson rate. `lambda <= 0` returns `0L` and draws
  nothing. Errors above a documented ceiling (30 Morgans / 3000 cM)
  rather than spinning `O(lambda)` uniforms.

## Value

Integer crossover count.
