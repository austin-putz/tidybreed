# Write the BLUPF90 genotype file (fixed-width ID + concatenated 0/1/2 string)

Write the BLUPF90 genotype file (fixed-width ID + concatenated 0/1/2
string)

## Usage

``` r
write_geno_file(pop, all_ped_ids, chip_name, eval_dir)
```

## Arguments

- pop:

  tidybreed_pop

- all_ped_ids:

  character vector of all pedigree animals

- chip_name:

  character chip name (e.g. "HD")

- eval_dir:

  path to evaluation folder

## Value

list with geno_file path, n_loci, id_width, geno_ids
