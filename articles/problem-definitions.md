# Mathematical definitions of the test problems

Use this page to find the residual equations, supported dimensions, and
standard starting point for any problem implemented by `funconstrain`.
The definitions follow Section 3 of [Moré, Garbow, and Hillstrom
(1981)](https://doi.org/10.1145/355934.355936) using one notation
throughout. The authors’ [Algorithm 566
source](https://netlib.org/toms/566.gz) is a second primary source for
resolving typographical and indexing ambiguities.

The formulas follow those sources; the dimensions and defaults describe
the current package. Where `funconstrain` supports a narrower range,
that range is shown. To run a problem, see [Getting started with
funconstrain](https://jlmelville.github.io/funconstrain/articles/getting-started.md)
or the [`funconstrain_problem()`
reference](https://jlmelville.github.io/funconstrain/reference/funconstrain_catalog.md).

For $`x=(x_1,\ldots,x_n)\in\mathbb{R}^n`$, every problem below defines
$`m`$ residuals. The scalar callback returned by the package is

``` math
F(x)=\sum_{i=1}^{m}r_i(x)^2,
\qquad \texttt{problem\$fn(x)}=F(x).
```

Indices in the equations are one-based. A sum over an empty index range
is zero. $`x^{(0)}`$ denotes the standard starting point.

## Problem index

| Number | Problem | Default n | Default m |
|---:|:---|---:|---:|
| 1 | [Rosenbrock Function](#problem-01-rosen) | 2 | 2 |
| 2 | [Freudenstein and Roth Function](#problem-02-freud-roth) | 2 | 2 |
| 3 | [Powell Badly Scaled Function](#problem-03-powell-bs) | 2 | 2 |
| 4 | [Brown Badly Scaled Function](#problem-04-brown-bs) | 2 | 3 |
| 5 | [Beale Function](#problem-05-beale) | 2 | 3 |
| 6 | [Jennrich and Sampson Function](#problem-06-jenn-samp) | 2 | 10 |
| 7 | [Helical Valley Function](#problem-07-helical) | 3 | 3 |
| 8 | [Bard Function](#problem-08-bard) | 3 | 15 |
| 9 | [Gaussian Function](#problem-09-gauss) | 3 | 15 |
| 10 | [Meyer Function](#problem-10-meyer) | 3 | 16 |
| 11 | [Gulf Research and Development Function](#problem-11-gulf) | 3 | 99 |
| 12 | [Box Three-Dimensional Function](#problem-12-box-3d) | 3 | 20 |
| 13 | [Powell Singular Function](#problem-13-powell-s) | 4 | 4 |
| 14 | [Wood Function](#problem-14-wood) | 4 | 6 |
| 15 | [Kowalik and Osborne Function](#problem-15-kow-osb) | 4 | 11 |
| 16 | [Brown and Dennis Function](#problem-16-brown-den) | 4 | 20 |
| 17 | [Osborne 1 Function](#problem-17-osborne-1) | 5 | 33 |
| 18 | [Biggs EXP6 Function](#problem-18-biggs-exp6) | 6 | 13 |
| 19 | [Osborne 2 Function](#problem-19-osborne-2) | 11 | 65 |
| 20 | [Watson Function](#problem-20-watson) | 6 | 31 |
| 21 | [Extended Rosenbrock Function](#problem-21-ex-rosen) | 8 | 8 |
| 22 | [Extended Powell Function](#problem-22-ex-powell) | 20 | 20 |
| 23 | [Penalty Function I](#problem-23-penalty-1) | 25 | 26 |
| 24 | [Penalty Function II](#problem-24-penalty-2) | 25 | 50 |
| 25 | [Variably Dimensioned Function](#problem-25-var-dim) | 30 | 32 |
| 26 | [Trigonometric Function](#problem-26-trigon) | 30 | 30 |
| 27 | [Brown Almost-Linear Function](#problem-27-brown-al) | 30 | 30 |
| 28 | [Discrete Boundary Value Function](#problem-28-disc-bv) | 35 | 35 |
| 29 | [Discrete Integral Equation Function](#problem-29-disc-ie) | 35 | 35 |
| 30 | [Broyden Tridiagonal Function](#problem-30-broyden-tri) | 40 | 40 |
| 31 | [Broyden Banded Function](#problem-31-broyden-band) | 40 | 40 |
| 32 | [Linear Function - Full Rank](#problem-32-linfun-fr) | 45 | 100 |
| 33 | [Linear Function - Rank 1](#problem-33-linfun-r1) | 45 | 100 |
| 34 | [Linear Function - Rank 1 with Zero Columns and Rows](#problem-34-linfun-r1z) | 45 | 100 |
| 35 | [Chebyquad Function](#problem-35-chebyquad) | 50 | 50 |

## 1. Rosenbrock Function

**Dimensions:** $`n=2`$, $`m=2`$  
**Factory:**
[`rosen()`](https://jlmelville.github.io/funconstrain/reference/rosen.md)

``` math
\begin{aligned}
r_1(x)&=10(x_2-x_1^2),\\
r_2(x)&=1-x_1.
\end{aligned}
```

**Start:** $`x^{(0)}=(-1.2,1)`$

## 2. Freudenstein and Roth Function

**Dimensions:** $`n=2`$, $`m=2`$  
**Factory:**
[`freud_roth()`](https://jlmelville.github.io/funconstrain/reference/freud_roth.md)

``` math
\begin{aligned}
r_1(x)&=-13+x_1+\{(5-x_2)x_2-2\}x_2,\\
r_2(x)&=-29+x_1+\{(1+x_2)x_2-14\}x_2.
\end{aligned}
```

**Start:** $`x^{(0)}=(1/2,-2)`$

## 3. Powell Badly Scaled Function

**Dimensions:** $`n=2`$, $`m=2`$  
**Factory:**
[`powell_bs()`](https://jlmelville.github.io/funconstrain/reference/powell_bs.md)

``` math
\begin{aligned}
r_1(x)&=10^4x_1x_2-1,\\
r_2(x)&=e^{-x_1}+e^{-x_2}-1.0001.
\end{aligned}
```

**Start:** $`x^{(0)}=(0,1)`$

## 4. Brown Badly Scaled Function

**Dimensions:** $`n=2`$, $`m=3`$  
**Factory:**
[`brown_bs()`](https://jlmelville.github.io/funconstrain/reference/brown_bs.md)

``` math
\begin{aligned}
r_1(x)&=x_1-10^6,\\
r_2(x)&=x_2-2\times10^{-6},\\
r_3(x)&=x_1x_2-2.
\end{aligned}
```

**Start:** $`x^{(0)}=(1,1)`$

## 5. Beale Function

**Dimensions:** $`n=2`$, $`m=3`$  
**Factory:**
[`beale()`](https://jlmelville.github.io/funconstrain/reference/beale.md)

For $`i=1,2,3`$, let $`(c_1,c_2,c_3)=(1.5,2.25,2.625)`$. Then

``` math
r_i(x)=c_i-x_1(1-x_2^i).
```

**Start:** $`x^{(0)}=(1,1)`$

## 6. Jennrich and Sampson Function

**Dimensions:** $`n=2`$, $`m\ge2`$; default $`m=10`$  
**Factory:**
[`jenn_samp()`](https://jlmelville.github.io/funconstrain/reference/jenn_samp.md)

For $`i=1,\ldots,m`$,

``` math
r_i(x)=2+2i-e^{ix_1}-e^{ix_2}.
```

**Start:** $`x^{(0)}=(0.3,0.4)`$

## 7. Helical Valley Function

**Dimensions:** $`n=3`$, $`m=3`$  
**Factory:**
[`helical()`](https://jlmelville.github.io/funconstrain/reference/helical.md)

Away from $`(x_1,x_2)=(0,0)`$, define

``` math
\theta(x_1,x_2)=
\begin{cases}
\dfrac{1}{2\pi}\tan^{-1}(x_2/x_1),&x_1>0,\\[4pt]
\dfrac{1}{2\pi}\tan^{-1}(x_2/x_1)+\dfrac12,&x_1<0,\\[4pt]
\dfrac14\operatorname{sgn}(x_2),&x_1=0.
\end{cases}
```

The residuals are

``` math
\begin{aligned}
r_1(x)&=10\{x_3-10\theta(x_1,x_2)\},\\
r_2(x)&=10\{(x_1^2+x_2^2)^{1/2}-1\},\\
r_3(x)&=x_3.
\end{aligned}
```

**Implementation note:** The package rejects the origin, where the
objective is undefined.

**Start:** $`x^{(0)}=(-1,0,0)`$

## 8. Bard Function

**Dimensions:** $`n=3`$, $`m=15`$  
**Factory:**
[`bard()`](https://jlmelville.github.io/funconstrain/reference/bard.md)

For $`i=1,\ldots,15`$, set $`v_i=16-i`$, $`w_i=\min(i,v_i)`$, and

``` math
\begin{aligned}
(y_1,\ldots,y_{15})={}&(0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,\\
&0.37,0.58,0.73,0.96,1.34,2.10,4.39).
\end{aligned}
```

Then

``` math
r_i(x)=y_i-\left(x_1+\frac{i}{v_ix_2+w_ix_3}\right).
```

**Start:** $`x^{(0)}=(1,1,1)`$

## 9. Gaussian Function

**Dimensions:** $`n=3`$, $`m=15`$  
**Factory:**
[`gauss()`](https://jlmelville.github.io/funconstrain/reference/gauss.md)

For $`i=1,\ldots,15`$, set $`t_i=(8-i)/2`$ and

``` math
\begin{aligned}
(y_1,\ldots,y_{15})={}&(0.0009,0.0044,0.0175,0.0540,0.1295,\\
&0.2420,0.3521,0.3989,0.3521,0.2420,\\
&0.1295,0.0540,0.0175,0.0044,0.0009).
\end{aligned}
```

Then

``` math
r_i(x)=x_1\exp\left\{-\frac{x_2}{2}(t_i-x_3)^2\right\}-y_i.
```

**Start:** $`x^{(0)}=(0.4,1,0)`$

## 10. Meyer Function

**Dimensions:** $`n=3`$, $`m=16`$  
**Factory:**
[`meyer()`](https://jlmelville.github.io/funconstrain/reference/meyer.md)

For $`i=1,\ldots,16`$, set $`t_i=45+5i`$ and

``` math
\begin{aligned}
(y_1,\ldots,y_{16})={}&(34780,28610,23650,19630,16370,13720,11540,9744,\\
&8261,7030,6005,5147,4427,3820,3307,2872).
\end{aligned}
```

Then

``` math
r_i(x)=x_1\exp\left(\frac{x_2}{t_i+x_3}\right)-y_i.
```

**Start:** $`x^{(0)}=(0.02,4000,250)`$

## 11. Gulf Research and Development Function

**Dimensions:** $`n=3`$, $`3\le m\le100`$; default $`m=99`$  
**Factory:**
[`gulf()`](https://jlmelville.github.io/funconstrain/reference/gulf.md)

For $`i=1,\ldots,m`$, define

``` math
t_i=\frac{i}{100},\qquad
y_i=25+\{-50\log(t_i)\}^{2/3}.
```

The residual is

``` math
r_i(x)=\exp\left\{-\frac{|x_2-y_i|^{x_3}}{x_1}\right\}-t_i.
```

**Source note:** The published equation has a typographical glyph where
the subtraction in $`|x_2-y_i|`$ should appear; the corrected form above
agrees with the authors’ Algorithm 566 source and the package.

**Start:** $`x^{(0)}=(5,2.5,0.15)`$

## 12. Box Three-Dimensional Function

**Dimensions:** $`n=3`$, $`m\ge3`$; default $`m=20`$  
**Factory:**
[`box_3d()`](https://jlmelville.github.io/funconstrain/reference/box_3d.md)

For $`i=1,\ldots,m`$, let $`t_i=i/10`$. Then

``` math
r_i(x)=e^{-t_ix_1}-e^{-t_ix_2}-x_3(e^{-t_i}-e^{-i}).
```

**Start:** $`x^{(0)}=(0,10,20)`$

## 13. Powell Singular Function

**Dimensions:** $`n=4`$, $`m=4`$  
**Factory:**
[`powell_s()`](https://jlmelville.github.io/funconstrain/reference/powell_s.md)

``` math
\begin{aligned}
r_1(x)&=x_1+10x_2,\\
r_2(x)&=\sqrt5(x_3-x_4),\\
r_3(x)&=(x_2-2x_3)^2,\\
r_4(x)&=\sqrt{10}(x_1-x_4)^2.
\end{aligned}
```

**Start:** $`x^{(0)}=(3,-1,0,1)`$

## 14. Wood Function

**Dimensions:** $`n=4`$, $`m=6`$  
**Factory:**
[`wood()`](https://jlmelville.github.io/funconstrain/reference/wood.md)

``` math
\begin{aligned}
r_1(x)&=10(x_2-x_1^2),&r_2(x)&=1-x_1,\\
r_3(x)&=\sqrt{90}(x_4-x_3^2),&r_4(x)&=1-x_3,\\
r_5(x)&=\sqrt{10}(x_2+x_4-2),&
r_6(x)&=(x_2-x_4)/\sqrt{10}.
\end{aligned}
```

**Start:** $`x^{(0)}=(-3,-1,-3,-1)`$

## 15. Kowalik and Osborne Function

**Dimensions:** $`n=4`$, $`m=11`$  
**Factory:**
[`kow_osb()`](https://jlmelville.github.io/funconstrain/reference/kow_osb.md)

The data are

``` math
\begin{aligned}
(y_1,\ldots,y_{11})={}&(0.1957,0.1947,0.1735,0.1600,0.0844,0.0627,\\
&0.0456,0.0342,0.0323,0.0235,0.0246),\\
(u_1,\ldots,u_{11})={}&(4,2,1,0.5,0.25,0.167,0.125,0.1,\\
&0.0833,0.0714,0.0625).
\end{aligned}
```

For $`i=1,\ldots,11`$,

``` math
r_i(x)=y_i-x_1\frac{u_i(u_i+x_2)}{u_i(u_i+x_3)+x_4}.
```

**Start:** $`x^{(0)}=(0.25,0.39,0.415,0.39)`$

## 16. Brown and Dennis Function

**Dimensions:** $`n=4`$, $`m\ge4`$; default $`m=20`$  
**Factory:**
[`brown_den()`](https://jlmelville.github.io/funconstrain/reference/brown_den.md)

For $`i=1,\ldots,m`$, let $`t_i=i/5`$ and define

``` math
\begin{aligned}
a_i(x)&=x_1+t_ix_2-e^{t_i},\\
b_i(x)&=x_3+x_4\sin(t_i)-\cos(t_i),\\
r_i(x)&=a_i(x)^2+b_i(x)^2.
\end{aligned}
```

**Objective note:** This problem squares $`a_i^2+b_i^2`$ once more in
the shared objective $`F(x)=\sum_i r_i(x)^2`$.

**Start:** $`x^{(0)}=(25,5,-5,1)`$

## 17. Osborne 1 Function

**Dimensions:** $`n=5`$, $`m=33`$  
**Factory:**
[`osborne_1()`](https://jlmelville.github.io/funconstrain/reference/osborne_1.md)

Let $`t_i=10(i-1)`$ and define the observations in consecutive blocks:

``` math
\begin{aligned}
(y_1,\ldots,y_8)&=(0.844,0.908,0.932,0.936,0.925,0.908,0.881,0.850),\\
(y_9,\ldots,y_{16})&=(0.818,0.784,0.751,0.718,0.685,0.658,0.628,0.603),\\
(y_{17},\ldots,y_{24})&=(0.580,0.558,0.538,0.522,0.506,0.490,0.478,0.467),\\
(y_{25},\ldots,y_{32})&=(0.457,0.448,0.438,0.431,0.424,0.420,0.414,0.411),\\
y_{33}&=0.406.
\end{aligned}
```

For $`i=1,\ldots,33`$,

``` math
r_i(x)=y_i-\{x_1+x_2e^{-t_ix_4}+x_3e^{-t_ix_5}\}.
```

**Start:** $`x^{(0)}=(0.5,1.5,-1,0.01,0.02)`$

## 18. Biggs EXP6 Function

**Dimensions:** $`n=6`$, $`m\ge6`$; default $`m=13`$  
**Factory:**
[`biggs_exp6()`](https://jlmelville.github.io/funconstrain/reference/biggs_exp6.md)

For $`i=1,\ldots,m`$, let $`t_i=i/10`$ and

``` math
y_i=e^{-t_i}-5e^{-10t_i}+3e^{-4t_i}.
```

Then

``` math
r_i(x)=x_3e^{-t_ix_1}-x_4e^{-t_ix_2}+x_6e^{-t_ix_5}-y_i.
```

**Start:** $`x^{(0)}=(1,2,1,1,1,1)`$

## 19. Osborne 2 Function

**Dimensions:** $`n=11`$, $`m=65`$  
**Factory:**
[`osborne_2()`](https://jlmelville.github.io/funconstrain/reference/osborne_2.md)

Let $`t_i=(i-1)/10`$. The observations, in consecutive blocks, are

``` math
\begin{aligned}
(y_1,\ldots,y_8)&=(1.366,1.191,1.112,1.013,0.991,0.885,0.831,0.847),\\
(y_9,\ldots,y_{16})&=(0.786,0.725,0.746,0.679,0.608,0.655,0.616,0.606),\\
(y_{17},\ldots,y_{24})&=(0.602,0.626,0.651,0.724,0.649,0.649,0.694,0.644),\\
(y_{25},\ldots,y_{32})&=(0.624,0.661,0.612,0.558,0.533,0.495,0.500,0.423),\\
(y_{33},\ldots,y_{40})&=(0.395,0.375,0.372,0.391,0.396,0.405,0.428,0.429),\\
(y_{41},\ldots,y_{48})&=(0.523,0.562,0.607,0.653,0.672,0.708,0.633,0.668),\\
(y_{49},\ldots,y_{56})&=(0.645,0.632,0.591,0.559,0.597,0.625,0.739,0.710),\\
(y_{57},\ldots,y_{64})&=(0.729,0.720,0.636,0.581,0.428,0.292,0.162,0.098),\\
y_{65}&=0.054.
\end{aligned}
```

For $`i=1,\ldots,65`$,

``` math
\begin{aligned}
r_i(x)=y_i-&\{x_1e^{-t_ix_5}+x_2e^{-x_6(t_i-x_9)^2}\\
&+x_3e^{-x_7(t_i-x_{10})^2}+x_4e^{-x_8(t_i-x_{11})^2}\}.
\end{aligned}
```

**Start:** $`x^{(0)}=(1.3,0.65,0.65,0.7,0.6,3,5,7,2,4.5,5.5)`$

## 20. Watson Function

**Dimensions:** $`2\le n\le31`$, $`m=31`$; default $`n=6`$  
**Factory:**
[`watson()`](https://jlmelville.github.io/funconstrain/reference/watson.md)

For $`i=1,\ldots,29`$, define

``` math
\begin{aligned}
t_i&=i/29,\\
a_i(x)&=\sum_{j=2}^{n}(j-1)x_jt_i^{j-2},\\
b_i(x)&=\sum_{j=1}^{n}x_jt_i^{j-1}.
\end{aligned}
```

The residuals are

``` math
\begin{aligned}
r_i(x)&=a_i(x)-b_i(x)^2-1,&&i=1,\ldots,29,\\
r_{30}(x)&=x_1,\\
r_{31}(x)&=x_2-x_1^2-1.
\end{aligned}
```

**Start:** $`x^{(0)}=(0,\ldots,0)`$

## 21. Extended Rosenbrock Function

**Dimensions:** even $`n\ge2`$, $`m=n`$; default $`n=8`$  
**Factory:**
[`ex_rosen()`](https://jlmelville.github.io/funconstrain/reference/ex_rosen.md)

For $`k=1,\ldots,n/2`$,

``` math
\begin{aligned}
r_{2k-1}(x)&=10(x_{2k}-x_{2k-1}^2),\\
r_{2k}(x)&=1-x_{2k-1}.
\end{aligned}
```

**Start:** Repeat $`(-1.2,1)`$ exactly $`n/2`$ times.

## 22. Extended Powell Function

**Dimensions:** $`n\ge4`$ and divisible by 4, $`m=n`$; default
$`n=20`$  
**Factory:**
[`ex_powell()`](https://jlmelville.github.io/funconstrain/reference/ex_powell.md)

For $`k=1,\ldots,n/4`$, let $`a=4k-3`$. Then

``` math
\begin{aligned}
r_a(x)&=x_a+10x_{a+1},\\
r_{a+1}(x)&=\sqrt5(x_{a+2}-x_{a+3}),\\
r_{a+2}(x)&=(x_{a+1}-2x_{a+2})^2,\\
r_{a+3}(x)&=\sqrt{10}(x_a-x_{a+3})^2.
\end{aligned}
```

**Naming note:** MGH calls this the “Extended Powell singular” problem;
the package title is shortened.

**Start:** Repeat $`(3,-1,0,1)`$ exactly $`n/4`$ times.

## 23. Penalty Function I

**Dimensions:** $`n\ge1`$, $`m=n+1`$; default $`n=25`$  
**Factory:**
[`penalty_1()`](https://jlmelville.github.io/funconstrain/reference/penalty_1.md)

With $`a=10^{-5}`$,

``` math
\begin{aligned}
r_j(x)&=\sqrt{a}(x_j-1),&&j=1,\ldots,n,\\
r_{n+1}(x)&=\sum_{j=1}^{n}x_j^2-\frac14.
\end{aligned}
```

**Start:** $`x_j^{(0)}=j`$

## 24. Penalty Function II

**Dimensions:** $`n\ge1`$, $`m=2n`$; default $`n=25`$  
**Factory:**
[`penalty_2()`](https://jlmelville.github.io/funconstrain/reference/penalty_2.md)

With $`a=10^{-5}`$, the residuals are

``` math
\begin{aligned}
r_1(x)&=x_1-0.2,\\
r_i(x)&=\sqrt a\{e^{x_i/10}+e^{x_{i-1}/10}-e^{i/10}-e^{(i-1)/10}\},
&&i=2,\ldots,n,\\
r_{n+i-1}(x)&=\sqrt a\{e^{x_i/10}-e^{-0.1}\},
&&i=2,\ldots,n,\\
r_{2n}(x)&=\sum_{j=1}^{n}(n-j+1)x_j^2-1.
\end{aligned}
```

**Edge case:** For $`n=1`$, the two middle ranges are empty.

**Start:** $`x_j^{(0)}=1/2`$

## 25. Variably Dimensioned Function

**Dimensions:** $`n\ge1`$, $`m=n+2`$; default $`n=30`$  
**Factory:**
[`var_dim()`](https://jlmelville.github.io/funconstrain/reference/var_dim.md)

Define $`s(x)=\sum_{j=1}^{n}j(x_j-1)`$. Then

``` math
\begin{aligned}
r_j(x)&=x_j-1,&&j=1,\ldots,n,\\
r_{n+1}(x)&=s(x),\\
r_{n+2}(x)&=s(x)^2.
\end{aligned}
```

**Start:** $`x_j^{(0)}=1-j/n`$

## 26. Trigonometric Function

**Dimensions:** $`n\ge1`$, $`m=n`$; default $`n=30`$  
**Factory:**
[`trigon()`](https://jlmelville.github.io/funconstrain/reference/trigon.md)

For $`i=1,\ldots,n`$,

``` math
r_i(x)=n-\sum_{j=1}^{n}\cos(x_j)+i\{1-\cos(x_i)\}-\sin(x_i).
```

**Start:** $`x_j^{(0)}=1/n`$

## 27. Brown Almost-Linear Function

**Dimensions:** $`n\ge1`$, $`m=n`$; default $`n=30`$  
**Factory:**
[`brown_al()`](https://jlmelville.github.io/funconstrain/reference/brown_al.md)

``` math
\begin{aligned}
r_i(x)&=x_i+\sum_{j=1}^{n}x_j-(n+1),&&i=1,\ldots,n-1,\\
r_n(x)&=\prod_{j=1}^{n}x_j-1.
\end{aligned}
```

**Edge case:** For $`n=1`$, the first range is empty.

**Start:** $`x_j^{(0)}=1/2`$

## 28. Discrete Boundary Value Function

**Dimensions:** $`n\ge1`$, $`m=n`$; default $`n=35`$  
**Factory:**
[`disc_bv()`](https://jlmelville.github.io/funconstrain/reference/disc_bv.md)

Let $`h=1/(n+1)`$, $`t_i=ih`$, and set the boundary values
$`x_0=x_{n+1}=0`$. For $`i=1,\ldots,n`$,

``` math
r_i(x)=2x_i-x_{i-1}-x_{i+1}+\frac{h^2}{2}(x_i+t_i+1)^3.
```

**Start:** $`x_i^{(0)}=t_i(t_i-1)`$

## 29. Discrete Integral Equation Function

**Dimensions:** $`n\ge1`$, $`m=n`$; default $`n=35`$  
**Factory:**
[`disc_ie()`](https://jlmelville.github.io/funconstrain/reference/disc_ie.md)

Let $`h=1/(n+1)`$, $`t_i=ih`$, and $`q_j(x)=(x_j+t_j+1)^3`$. For
$`i=1,\ldots,n`$,

``` math
r_i(x)=x_i+\frac{h}{2}\left\{
(1-t_i)\sum_{j=1}^{i}t_jq_j(x)
+t_i\sum_{j=i+1}^{n}(1-t_j)q_j(x)
\right\}.
```

**Start:** $`x_i^{(0)}=t_i(t_i-1)`$

## 30. Broyden Tridiagonal Function

**Dimensions:** $`n\ge1`$, $`m=n`$; default $`n=40`$  
**Factory:**
[`broyden_tri()`](https://jlmelville.github.io/funconstrain/reference/broyden_tri.md)

Set $`x_0=x_{n+1}=0`$. For $`i=1,\ldots,n`$,

``` math
r_i(x)=(3-2x_i)x_i-x_{i-1}-2x_{i+1}+1.
```

**Start:** $`x^{(0)}=(-1,\ldots,-1)`$

## 31. Broyden Banded Function

**Dimensions:** $`n\ge1`$, $`m=n`$; default $`n=40`$  
**Factory:**
[`broyden_band()`](https://jlmelville.github.io/funconstrain/reference/broyden_band.md)

For each $`i=1,\ldots,n`$, define

``` math
B_i=\{j:\max(1,i-5)\le j\le\min(n,i+1),\ j\ne i\}.
```

Then

``` math
r_i(x)=x_i(2+5x_i^2)+1-\sum_{j\in B_i}x_j(1+x_j).
```

**Start:** $`x^{(0)}=(-1,\ldots,-1)`$

## 32. Linear Function - Full Rank

**Dimensions:** $`n\ge1`$, $`m\ge n`$; defaults $`n=45`$, $`m=100`$  
**Factory:**
[`linfun_fr()`](https://jlmelville.github.io/funconstrain/reference/linfun_fr.md)

Define $`c(x)=1+(2/m)\sum_{j=1}^{n}x_j`$. Then

``` math
r_i(x)=
\begin{cases}
x_i-c(x),&i=1,\ldots,n,\\
-c(x),&i=n+1,\ldots,m.
\end{cases}
```

**Start:** $`x^{(0)}=(1,\ldots,1)`$

## 33. Linear Function - Rank 1

**Dimensions:** $`n\ge1`$, $`m\ge n`$; defaults $`n=45`$, $`m=100`$  
**Factory:**
[`linfun_r1()`](https://jlmelville.github.io/funconstrain/reference/linfun_r1.md)

For $`i=1,\ldots,m`$,

``` math
r_i(x)=i\sum_{j=1}^{n}jx_j-1.
```

**Start:** $`x^{(0)}=(1,\ldots,1)`$

## 34. Linear Function - Rank 1 with Zero Columns and Rows

**Dimensions:** $`n\ge1`$, $`m\ge n`$; defaults $`n=45`$, $`m=100`$  
**Factory:**
[`linfun_r1z()`](https://jlmelville.github.io/funconstrain/reference/linfun_r1z.md)

Define $`s(x)=\sum_{j=2}^{n-1}jx_j`$. The residuals are

``` math
\begin{aligned}
r_1(x)&=-1,\\
r_i(x)&=(i-1)s(x)-1,&&i=2,\ldots,m-1,\\
r_m(x)&=-1.
\end{aligned}
```

**Edge case:** For $`n\le2`$, the sum defining $`s(x)`$ is empty.

**Start:** $`x^{(0)}=(1,\ldots,1)`$

## 35. Chebyquad Function

**Dimensions:** $`n\ge1`$, $`m=n`$; default $`n=50`$  
**Factory:**
[`chebyquad()`](https://jlmelville.github.io/funconstrain/reference/chebyquad.md)

Let $`T_i^*(u)=T_i(2u-1)`$ be the Chebyshev polynomial shifted to
$`[0,1]`$, with

``` math
\begin{aligned}
T_0^*(u)&=1,\\
T_1^*(u)&=2u-1,\\
T_i^*(u)&=2(2u-1)T_{i-1}^*(u)-T_{i-2}^*(u).
\end{aligned}
```

For $`i=1,\ldots,n`$,

``` math
r_i(x)=\frac1n\sum_{j=1}^{n}T_i^*(x_j)-c_i,
\qquad
c_i=
\begin{cases}
0,&i\text{ odd},\\
-1/(i^2-1),&i\text{ even}.
\end{cases}
```

**Package note:** MGH defines Chebyquad for $`m\ge n`$, with residual
indices through $`m`$. The current package has no independent `m`
control and implements only $`m=n`$; the narrower package domain is the
supported contract here.

**Start:** $`x_j^{(0)}=j/(n+1)`$
