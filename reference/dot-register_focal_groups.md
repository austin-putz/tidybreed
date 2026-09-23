# The group value of each focal individual, registered as `__ap_focal_groups`

Focals with a `NULL` group value are left out of the view. The caller
unregisters the view.

## Usage

``` r
.register_focal_groups(conn, focal_ids, group_column, group_table, what)
```

## Value

Logical: which focals have a group value.
