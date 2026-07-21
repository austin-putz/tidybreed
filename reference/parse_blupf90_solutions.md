# Parse the blupf90+ solutions file and return EBV tibble

With OPTION origID the solutions file has 5 columns: trait_num
effect_num level original_id solution

## Usage

``` r
parse_blupf90_solutions(
  eval_dir,
  trait_name,
  animal_effect_num,
  all_ped_ids,
  model,
  eval_nums
)
```

## Arguments

- eval_dir:

  path to evaluation folder

- trait_name:

  character vector of trait names (in model order)

- animal_effect_num:

  integer effect number for the animal random effect

- all_ped_ids:

  character vector of all pedigree animal IDs (used to distinguish
  animal solutions from other random effect solutions)

- model:

  character model label

- eval_nums:

  named integer vector of eval_number per trait name

## Value

tibble: id_ind, trait_name, model, ebv_value, acc, se, eval_number
