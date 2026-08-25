# Using funconstrain with mize

This article shows a small, repeatable integration between
`funconstrain` and [`mize`](https://jlmelville.github.io/mize/). It
resolves exact problem instances, adapts their callback names, and keeps
the optimizer’s native termination record separate from independent
final-quality checks.

The first version supports GitHub development `mize` version 0.2.5.9001
or later, rather than CRAN `mize` 0.2.5. Install the supported version
before running the examples:

``` r

pak::pak("jlmelville/mize")
```

## Load the shared harness

The sourceable companion used here is installed with `funconstrain`. It
remains an unexported example rather than a package API:

``` r

library(funconstrain)

harness_path <- system.file(
  "examples",
  "mize-harness.R",
  package = "funconstrain"
)
stopifnot(nzchar(harness_path))
source(harness_path)
```

At the callback boundary, `mize` uses `hs` for a Hessian while
`funconstrain` uses `he`. The adapter preserves the objective, gradient,
and combined callback and changes only that Hessian name:

``` r

problem <- funconstrain_problem("rosen")
names(as_mize_fg(problem))
#> [1] "fn" "gr" "fg" "hs"
```

The companion also wraps `fn`, `gr`, `fg`, and `hs` separately. These
physical call counts answer a different question from `mize`’s native
`nf` and `ng` evaluation totals, especially when one combined `fg` call
supplies both an objective and a gradient.

## Run one ordinary solve

`run_mize_problem()` owns the resolved starting point, callback list,
and method argument. Solver controls are supplied explicitly in a named
list, so the run record retains the choices that differ from `mize`
defaults.

``` r

ordinary <- run_mize_problem(
  "rosen",
  method = "BFGS",
  controls = list(
    max_iter = 200,
    abs_tol = NULL,
    rel_tol = NULL,
    ginf_tol = 1e-6,
    step_tol = NULL
  )
)

data.frame(
  status = ordinary$outcome$status,
  termination = ordinary$outcome$terminate$what,
  converged = ordinary$outcome$converged,
  final_objective = signif(ordinary$quality$final$objective, 3),
  gradient_inf_norm = signif(
    ordinary$quality$final$gradient_inf_norm,
    3
  )
)
#>      status termination converged final_objective gradient_inf_norm
#> 1 converged    ginf_tol      TRUE        5.59e-16          6.23e-07
```

`status`, `message`, `converged`, and `terminate` above are native
`mize` observations. The final objective and gradient infinity norm are
evaluated afterward through the original resolved callbacks. They are
independent audit observations, not a reinterpretation of the native
termination status.

The stored objective reference is kept as a structured record. Compare
with its value only when its status is `"applicable"`; a stored
minimizer is not a required equality target for a successful run.

``` r

ordinary$reference[c("status", "value", "source_configuration")]
#> $status
#> [1] "applicable"
#> 
#> $value
#> [1] 0
#> 
#> $source_configuration
#> list()
```

## Preserve a hard-budget outcome

A hard callback budget is a valid experimental outcome, not a malformed
row. This run permits only one combined objective-or-gradient
evaluation:

``` r

budget_limited <- run_mize_problem(
  "rosen",
  method = "BFGS",
  controls = list(
    max_iter = 3,
    max_fg = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
)

data.frame(
  status = budget_limited$outcome$status,
  termination = budget_limited$outcome$terminate$what,
  converged = budget_limited$outcome$converged,
  native_objective_available = !is.null(
    budget_limited$native_quality$objective
  ),
  audited_objective = signif(
    budget_limited$quality$final$objective,
    3
  )
)
#>             status termination converged native_objective_available
#> 1 budget_exhausted      max_fg     FALSE                      FALSE
#>   audited_objective
#> 1              24.2
```

Here `mize` legitimately omits its final objective because another
optimizer callback would violate the hard budget. The harness leaves
that native field absent. Its separately labelled audited objective is
evaluated only after the run, through the original callback and outside
the optimizer callback counts.

## Inspect work and provenance

Keep native evaluation totals and physical callback invocations as two
labelled views rather than trying to reconcile them into one number:

``` r

ordinary$work
#> $native
#> $native$nf
#> [1] 48
#> 
#> $native$ng
#> [1] 48
#> 
#> $native$iter
#> [1] 38
#> 
#> 
#> $physical
#> fn gr fg hs 
#>  1  1 47  0
```

The record also keeps the exact configurations and environment
provenance. The commit identifies the development `mize` source used for
this rendered run; it records provenance without turning the article
into a permanent commit pin.

``` r

ordinary$problem
#> $number
#> [1] 1
#> 
#> $name
#> [1] "rosen"
#> 
#> $title
#> [1] "Rosenbrock Function"
#> 
#> $configuration
#> $configuration$requested
#> $configuration$requested$n
#> NULL
#> 
#> $configuration$requested$m
#> NULL
#> 
#> 
#> $configuration$effective
#> $configuration$effective$n
#> [1] 2
#> 
#> $configuration$effective$m
#> [1] 2
ordinary$provenance
#> $r_version
#> [1] "R version 4.6.1 (2026-06-24)"
#> 
#> $platform
#> [1] "x86_64-pc-linux-gnu"
#> 
#> $funconstrain_version
#> [1] "0.1.1"
#> 
#> $mize_version
#> [1] "0.2.5.9001"
#> 
#> $mize_commit
#> [1] "6f975bc7897700ab0ba7c652c73c5ead8abd203e"
```

## Run a small explicit manifest

For a small integration check, declare every problem instance. Use `NA`
for a default dimension and explicit values for variable `n` or
configurable `m`:

``` r

manifest <- data.frame(
  name = c("rosen", "ex_rosen", "jenn_samp"),
  n = c(NA_integer_, 4L, NA_integer_),
  m = c(NA_integer_, NA_integer_, 20L)
)

suite <- run_mize_suite(
  manifest,
  method = "BFGS",
  controls = list(max_iter = 10)
)

suite_table <- do.call(rbind, lapply(suite$runs, function(run) {
  data.frame(
    problem = run$problem$name,
    effective_n = run$problem$configuration$effective$n,
    effective_m = run$problem$configuration$effective$m,
    status = run$outcome$status,
    final_objective = signif(run$quality$final$objective, 3)
  )
}))
rownames(suite_table) <- NULL
suite_table
#>     problem effective_n effective_m           status final_objective
#> 1     rosen           2           2 budget_exhausted            2.60
#> 2  ex_rosen           4           4 budget_exhausted            4.72
#> 3 jenn_samp           2          20 budget_exhausted        36200.00
```

This is intentionally a small integration and reproducibility harness.
It does not rank methods or prescribe line searches and stopping
controls. See the [`mize`
documentation](https://jlmelville.github.io/mize/) for optimizer
behavior and method selection; use this companion when exact
`funconstrain` problem resolution and an inspectable run record are the
task.
