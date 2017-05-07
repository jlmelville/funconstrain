# funconstrain

An R Package of Functions for Testing Unconstrained Numerical Optimization.

`funconstrain` is a pure R implementation of the 35 test functions in the paper
by [Moré, Garbow, and Hillstrom](https://doi.org/10.1145/355934.355936) useful
(to varying degrees) for testing unconstrained numerical optimization methods,
e.g. those implementing the likes of steepest descent, Newton, BFGS, L-BFGS, 
conjugate gradient and so on.

## Installing

```R
# if needed, install devtools:
# install.packages("devtools")
devtools::install_github("jlmelville/funconstrain")
library(funconstrain)
```

## Documentation

```R
package?funconstrain
```

## Examples

It's pretty simple. You call a function named after the test problem at hand,
and get back a list. That list contains: a function that implements the 
objective function; a function that implements the gradient; a function that
calculates the objective and the gradient in one go (if your fancy optimization
routine supports that); and a suggested starting point, which is also a function
if the test problem supports different dimensionalities, and is a plain numeric
vector otherwise.

```R
# The famous Rosenbrock function is a problem with two parameters
rbrock <- rosen()

# rbrock is a list containing function (fn), gradient (gr) and starting point
# (x0, a 2D numeric vector)
# Pass them to an optimization method:
res <- stats::optim(par = rbrock$x0, fn = rbrock$fn, gr = rbrock$gr, method = "L-BFGS-B")
# Or feel free to ignore the suggested starting point and use your own:
res <- stats::optim(par = c(1.2, 1.2), fn = rbrock$fn, gr = rbrock$gr, method = "L-BFGS-B")

# The Chebyquad function is defined for multiple parameters (any n > 0)
cheby <- chebyquad()

# To use different values of n, we provide it to the starting point x0, which
# is a function now that n can take multiple values for this test set.
# A five-parameter version:
res_n5 <- stats::optim(par = cheby$x0(n = 5), fn = cheby$fn, gr = cheby$gr, method = "L-BFGS-B")
# And a 10-parameter version:
res_n10 <- stats::optim(par = cheby$x0(n = 10), fn = cheby$fn, gr = cheby$gr, method = "L-BFGS-B")
```
The package and function documentation contain more examples.

## Why do this?

For testing numerical optimization routines, the go-to set of test problems is
[CUTEst](https://ccpforge.cse.rl.ac.uk/gf/project/cutest/wiki). However, if you
aren't compiling and linking to it directly, you'd have to write a parser for 
the SIF file format it uses. Neither of those possibilities appealed.

Instead, I re-implemented the functions as provided in the paper and also 
calculated analytical gradients for them. I did look to see if all 35 problems
were implemented in one place in R, but failed to find such a package.

## Are the functions correct?

There are unit tests for each test problem which ensure that:

* The analytical gradients match finite difference estimates at the suggested
starting point.
* If the location of a minima was given in the paper, that the analytical
gradient is close to zero at that location.
* If the location of a minima was given in the paper, that the objective
function has the correct value at that location.
* Running the L-BFGS or BFGS method as implemented in the `stats::optim`
function gets to the specified minimum.

## Is the implementation efficient?

Not really. My goal was correctness, and to make the code clear. Also perhaps, 
to be useful if anyone ever wants to translate these into other languages 
without having to know a lot of idiomatic R.

I have made use of vectorized arithmetic operation rather than explicit `for`
loops where possible as well as using functions like `sum`. Also, I am pretty
profligate in storing pre-computed vectors, trading off memory consumption for
clarity and potentially fast vectorized computations (I have not done any 
profiling). But I consciously eschewed the use of  `apply` `sweep` or other 
cleverness. 

I think I have elided the most gratuitious inefficiencies, such as unnecessary
recomputation of values inside loops.

## See also

* The aforementioned 
[CUTEst](https://ccpforge.cse.rl.ac.uk/gf/project/cutest/wiki). I believe all
or nearly all of the test problems in this package are implemented in CUTEst, 
but I make no representation that you will get the same results (if there are
any differences, assume it's a bug in `funconstrain`).
* I made ample use of the excellent 
[Derivative Calculator](http://www.derivative-calculator.net) to calculate the
analytical gradients.
* Shameless plug: I wrote this package to test 
[mize](https://github.com/jlmelville/mize), an R package for doing numerical
optimization.

## License

[MIT License](http://opensource.org/licenses/MIT).
