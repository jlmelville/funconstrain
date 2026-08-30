# Discover funconstrain test problems

`funconstrain_catalog()` lists the 35 Moré–Garbow–Hillstrom test
problems in paper order. `funconstrain_problem()` resolves one exact
factory name and an optional dimension request to a concrete,
solver-neutral problem.

## Usage

``` r
funconstrain_catalog()

funconstrain_problem(name, n = NULL, m = NULL)
```

## Arguments

- name:

  A single exact factory name from the `name` column of
  `funconstrain_catalog()`.

- n:

  Optional number of parameters. `NULL` uses the problem default.

- m:

  Optional number of summand functions. `NULL` uses the problem default
  or the value derived from `n`.

## Value

`funconstrain_catalog()` returns a data frame with one row per problem
and these columns:

- `number`, `name`, and `title` identify the problem;

- `n_kind` is either `"fixed"` or `"variable"`;

- `n_default` is the default parameter dimension;

- `m_configurable` reports whether `m` is a factory argument; and

- `m_default` is the effective `m` at `n_default`.

`funconstrain_problem()` returns a plain list with `fn`, `gr`, `he`, and
`fg` callbacks; a concrete numeric `x0`; and `problem`, `configuration`,
and `reference` metadata. The callbacks accept only vectors whose length
equals `configuration$effective$n`.

`configuration` records normalized explicit requests under `requested`
(omitted values remain `NULL`) and concrete integer values under
`effective`.

Each `fmin` and `xmin` reference record contains `value`, `status`,
`reason`, and `source_configuration`. Status is one of `"applicable"`,
`"not_applicable"`, `"unavailable"`, or `"unknown"`. A stored value is
returned even when it is not applicable so that its documented source
configuration remains inspectable; unavailable values are `NULL`.
Reasons currently distinguish documented applicability, a source
configuration mismatch, no single stored minimizer, and an unencoded
applicability rule. `source_configuration` is a partial rule: `n` and
`m` are exact values, `_min` and `_max` suffixes are inclusive bounds,
and an omitted component applies to every package-supported value.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41. <https://doi.org/10.1145/355934.355936>

## Examples

``` r
catalog <- funconstrain_catalog()
catalog[c(1, 20, 35), ]
#>    number      name               title   n_kind n_default m_configurable
#> 1       1     rosen Rosenbrock Function    fixed         2          FALSE
#> 20     20    watson     Watson Function variable         6          FALSE
#> 35     35 chebyquad  Chebyquad Function variable        50           TRUE
#>    m_default
#> 1          2
#> 20        31
#> 35        50

fixed <- funconstrain_problem("rosen")
fixed$fn(fixed$x0)
#> [1] 24.2

variable_n <- funconstrain_problem("chebyquad", n = 10)
length(variable_n$x0)
#> [1] 10

configurable_m <- funconstrain_problem("jenn_samp", m = 20)
configurable_m$configuration$effective
#> $n
#> [1] 2
#> 
#> $m
#> [1] 20
#> 

different_references <- funconstrain_problem("ex_rosen", n = 4)
different_references$reference$fmin$status
#> [1] "applicable"
different_references$reference$xmin$status
#> [1] "not_applicable"
```
