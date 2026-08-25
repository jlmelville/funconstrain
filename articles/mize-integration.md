# Using funconstrain with mize

`funconstrain` supplies the test problem;
[`mize`](https://jlmelville.github.io/mize/) supplies the optimizer. We
will fit Rosenbrock with BFGS, look at the optimization as it
progresses, and then use a small harness to repeat the same setup
consistently.

> **Version requirement.** These examples use GitHub development `mize`
> 0.2.5.9002 or later. Install it with:
>
> ``` r
>
> pak::pak("jlmelville/mize")
> ```

## Fit Rosenbrock with BFGS

Start by resolving Rosenbrock and passing its callbacks to `mize()`:

``` r

library(funconstrain)

problem <- funconstrain_problem("rosen")
fg <- list(
  fn = problem$fn,
  gr = problem$gr,
  fg = problem$fg,
  hs = problem$he
)

fit <- mize::mize(
  problem$x0,
  fg,
  method = "BFGS",
  max_iter = 200,
  abs_tol = NULL,
  rel_tol = NULL,
  ginf_tol = 1e-6,
  step_tol = NULL,
  store_progress = TRUE
)

data.frame(
  status = fit$status,
  stopped_by = fit$terminate$what,
  iterations = fit$iter,
  objective_evaluations = fit$nf,
  gradient_evaluations = fit$ng,
  hessian_evaluations = fit$nh,
  inverse_hessian_evaluations = fit$nhi,
  final_objective = signif(fit$f, 3),
  gradient_inf_norm = signif(fit$ginfn, 3)
)
#>      status stopped_by iterations objective_evaluations gradient_evaluations
#> 1 converged   ginf_tol         38                    48                   48
#>   hessian_evaluations inverse_hessian_evaluations final_objective
#> 1                   0                           0        5.59e-16
#>   gradient_inf_norm
#> 1          6.23e-07
```

`mize` reports that the gradient infinity-norm tolerance stopped this
run. The small final objective and gradient norm provide a useful check
that the returned point is good for this problem.

The callback list is a direct name-for-name match except for the
Hessian: `funconstrain` calls it `he`, while `mize` calls it `hs`. BFGS
does not use the Hessian; `hs` is included so the same callback list can
also be passed to Hessian-based methods. Consequently `fit$nh` is zero.
`fit$nhi` is also zero because this adapter does not supply an
inverse-Hessian callback; BFGS updates its approximation internally
rather than calling one.

## Follow the optimization

With `store_progress = TRUE`, `mize` records the initial state and then
one row per outer iteration. The objective trace shows the rapid early
decrease and the smaller improvements near the solution. Progress
storage does not request an otherwise-unneeded objective value for the
starting row, so the article evaluates that value once for the plot.
This display-only call is outside the objective count in `fit$nf`.

``` r

objective_trace <- fit$progress$f
objective_trace[[1L]] <- problem$fn(problem$x0)

plot(
  seq_along(objective_trace) - 1L,
  objective_trace,
  type = "l",
  log = "y",
  xlab = "Outer iteration",
  ylab = "Objective"
)
```

![A line plot on a logarithmic vertical scale showing the Rosenbrock
objective falling by many orders of magnitude over the BFGS
iterations.](mize-integration_files/figure-html/objective-trace-1.png)

BFGS objective trace for Rosenbrock. The initial objective is followed
by the value at each outer iteration.

The plot uses outer iteration to show the trajectory. The result table
above reports objective, gradient, Hessian, and inverse-Hessian
evaluations because methods can do different amounts of work inside one
iteration.

## Repeat runs consistently

The direct call is all that is needed for one problem. Repeating it
across problems also calls for consistent stopping controls and result
summaries. `funconstrain` therefore installs a small sourceable example
harness:

``` r

harness_path <- system.file(
  "examples",
  "mize-harness.R",
  package = "funconstrain"
)
stopifnot(nzchar(harness_path))
source(harness_path)
```

Declare the problems and dimensions you want to run. Use `NA` to accept
a default dimension:

``` r

problem_set <- data.frame(
  name = c("rosen", "ex_rosen", "jenn_samp"),
  n = c(NA_integer_, 4L, NA_integer_),
  m = c(NA_integer_, NA_integer_, 20L)
)

suite <- run_mize_suite(
  problem_set,
  method = "BFGS",
  controls = list(
    max_iter = 50,
    abs_tol = NULL,
    rel_tol = NULL,
    ginf_tol = 1e-6,
    step_tol = NULL
  )
)

suite_table <- do.call(rbind, lapply(suite$runs, function(run) {
  data.frame(
    problem = run$problem$name,
    n = run$problem$configuration$effective$n,
    m = run$problem$configuration$effective$m,
    status = run$outcome$status,
    objective_evaluations = run$work$native$nf,
    gradient_evaluations = run$work$native$ng,
    final_objective = signif(run$quality$final$objective, 3),
    gradient_inf_norm = signif(
      run$quality$final$gradient_inf_norm,
      3
    )
  )
}))
rownames(suite_table) <- NULL
suite_table
#>     problem n  m    status objective_evaluations gradient_evaluations
#> 1     rosen 2  2 converged                    48                   48
#> 2  ex_rosen 4  4 converged                    54                   54
#> 3 jenn_samp 2 20 converged                    61                   61
#>   final_objective gradient_inf_norm
#> 1        5.59e-16          6.23e-07
#> 2        9.65e-24          1.02e-11
#> 3        1.45e+03          4.70e-07
```

All three runs use the same method and stopping rule. The same harness
handles a fixed problem, a chosen parameter dimension, and a chosen
residual count. Objective values are problem-specific, so compare each
run with its own start or applicable reference rather than comparing
final objectives across rows.

The table uses the stopping status and evaluation totals reported by
`mize`. The harness also retains `nh` and `nhi` in `run$work$native`;
both are zero for these BFGS runs. The final objective and gradient norm
are checked afterward with the original problem callbacks, outside the
optimizer’s counts, so every row receives the same final-quality
calculation.

## Reproduce a run

Each harness result retains the supplied controls, problem dimensions,
package versions, platform, and development `mize` commit. A compact
summary is often enough for a report:

``` r

rosen_run <- suite$runs$rosen

data.frame(
  field = c(
    "Problem",
    "Dimensions",
    "Method",
    "R version",
    "funconstrain",
    "mize",
    "mize commit",
    "R platform"
  ),
  value = c(
    rosen_run$problem$title,
    paste0(
      "n = ", rosen_run$problem$configuration$effective$n,
      ", m = ", rosen_run$problem$configuration$effective$m
    ),
    rosen_run$method,
    rosen_run$provenance$r_version,
    rosen_run$provenance$funconstrain_version,
    rosen_run$provenance$mize_version,
    substr(rosen_run$provenance$mize_commit, 1L, 12L),
    rosen_run$provenance$platform
  ),
  row.names = NULL
)
#>          field                        value
#> 1      Problem          Rosenbrock Function
#> 2   Dimensions                 n = 2, m = 2
#> 3       Method                         BFGS
#> 4    R version R version 4.6.1 (2026-06-24)
#> 5 funconstrain                        0.1.1
#> 6         mize                   0.2.5.9002
#> 7  mize commit                 5ec90ec3dbc5
#> 8   R platform          x86_64-pc-linux-gnu
```

The result also retains the complete supplied control list in
`rosen_run$controls`.

For optimizer behavior and method selection, continue with the [`mize`
documentation](https://jlmelville.github.io/mize/). For problem
dimensions and stored references, see [Getting started with
funconstrain](https://jlmelville.github.io/funconstrain/articles/getting-started.md)
and [Mathematical definitions of the test
problems](https://jlmelville.github.io/funconstrain/articles/problem-definitions.md).
