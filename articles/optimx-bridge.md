# Using the optimx Bridge

The `optimx` bridge is optional. It is useful when you want to run the
`funconstrain` test functions through the solvers supported by `optimx`,
while keeping the test-function definitions in this package.

For the full contributed write-up, including background, the RFO
specification format, Hessian examples, and bounded-parameter examples,
see [Running funconstrain tests in package
optimx](https://jlmelville.github.io/funconstrain/articles/RunFunconstrainTests.md).

## Entry points

[`fufn()`](https://jlmelville.github.io/funconstrain/reference/fufn.md)
returns the data needed to run one test problem through `optimx`,
including the objective, gradient, Hessian, starting point, bounds, and
selected method metadata.

[`fufnrun()`](https://jlmelville.github.io/funconstrain/reference/fufnrun.md)
reads a four-line RFO specification file and runs the selected test
problems and methods. The package includes `vignettes/RFO.txt` as a
compact example:

``` text
testsink240408A.txt
1, 9, 9, 1, 6:8, 35
c("L-BFGS-B", "lbfgs", "lbfgsb3c", "lbfgs")
FALSE
```

The lines specify the sink file name, test problem numbers, `optimx`
methods, and whether experimental bounds constraints should be used.

## Optional packages

Install `optimx` before using
[`fufn()`](https://jlmelville.github.io/funconstrain/reference/fufn.md)
or
[`fufnrun()`](https://jlmelville.github.io/funconstrain/reference/fufnrun.md):

``` r

install.packages("optimx")
```

Some methods in the example RFO file also require `lbfgs` or `lbfgsb3c`:

``` r

install.packages(c("lbfgs", "lbfgsb3c"))
```

Then run the example specification from a local checkout:

``` r

library(funconstrain)
fufnrun("vignettes/RFO.txt")
```
