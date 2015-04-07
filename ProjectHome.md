## General Description ##
The studentized range distribution is commonly used for multiple comparison testing such as single step procedures like Tukey's HSD as well as step-down procedures such as the Ryan or Ryan-Einot-Gabriel-Welsch Q (REGWQ). Step-down procedures often require calculating q values at fractional values of alpha. Surprisingly few algorithms exist to calculate non-tabled values, and to my knowledge no other algorithms to estimate studentized range values have been implemented in Python.

This implementation builds off of the work by Gleason (1998, 1999). It uses tabled values of the _q_ distribution for finding approximations at non-tabled values. According to Gleason this method should be more accurate than the AS190 FORTRAN algorithm of Lund and Lund (1983) and works across a larger range (The AS190 only works from .9 <= p <= .99).

It is more efficient compared to the Copenhaver & Holland (1988) iterative algorithm (used by the _qtukey_ R function). To give you an idea of the performance difference, _qtukey_ took about 45 minutes to generate a test dataset of 1,216,000 points. Evaluating the same points using _qsturng_ took 181 seconds.

This implementation is also stable in ranges where the Copenhaver & Holland is not. For instance The C-H can run into problems when integrating with small degrees of freedom. The C-H algorithm is often implemented with 128 bit floats to help avoid such issues. The C-H algorithm does not provide any approximations for 1 degree of freedom. Since Harter's (1960) table provides these values when p >= .9 _qsturng_ can also provide approximations for 1 degree of freedom. The C-H also runs into problems when r values are high and p-values are low (at least R's version).

The project described was supported by NIH Grant Number P20 RR016454 from the INBRE Program of the National Center for Research Resources.
![http://qsturng-py.googlecode.com/files/qsturng%2Cdpi100.png](http://qsturng-py.googlecode.com/files/qsturng%2Cdpi100.png)
> Figure 1. _The x markers denote estimates obtained by Harter (1960) and/or the
> Copenhaver & Holland (1988) algorithm. The black traces show how this algorithm
> interpolates over p between the tabled values at v = 2 (top), v = 18 (middle),
> and v = infinity (bottom). The blue traces interpolate q over both p and r where r
> ranges from 21 to 29._

---

## Table of Contents ##
> 

---

## Function Descriptions ##
Implementation of Gleason's (1999) non-iterative quantile
studentized range approximation.
> ### qsturng(_p_, _r_, _v_) ###
> > Returns the q-value of the Studentized Range q-distribution as a
> > function of the probability (p), number of sample means  (r), and
> > the degrees of freedom (v).
> > > #### _p_- probability ####
> > > #### _r_- sample size for range ####
> > > #### _v_- degrees of freedom ####


> ### psturng(_q_, _r_, _v_) ###
> > returns the probability for the Studentized q-distribution where
> > the value q cooresponds to qsturng(1 - p, r, v).
> > > #### _q_- q statistic ####
> > > #### _r_- sample size for range ####
> > > #### _v_- degrees of freedom ####


---

## Some Usage Notes ##
  * The arguments _p_, _r_, and _v_ can be scalars or list-like objects (it has been vectorized using numpy.vectorize)
  * _qsturng_ accepts a range for p from .1 to .999. _psturng_ uses scalar minimization to obtain _p_. Returned _p_ values are consequently bound between .001 and .9. If .001 is returned _p_ should be interpreted as p <= .001
  * When .1 < p < .5 estimates are first-order in _p_/_q_ space (not that precise). This is to ensure stability of the _psturng_'s scalar minimization and to save on memory
  * Estimates between .5 and .9 should have increased accuracy over Gleason's Stata implementation due to the inclusion of table values at .675, .8, and .85.
  * Values of _r_ have fairly high precision between 2 and 100 and acceptable precision to 200. Values over 200 are not recommended
  * For p-values equal to .95 and over, degrees of freedom can range from 1 to infinity
  * For p-values less than .95, degrees of freedom must be greater than or equal to 2
  * When comparing directly to the C-H algorithm values may often differ at the fourth significant digit. _qsturng_ fits values to Harter's 1960 table which is slightly more conservative on about 20% of the values.

---

## Usage Examples ##
```
>>> from qsturng import qsturng, psturng
>>> qsturng(.9,5,6)
4.4364554589956198
>>> qs = qsturng([.8932, .9345,.9827], [4, 4, 4], [6, 6, 6])
>>> qs
array([ 3.98832389,  4.56835318,  6.26400894])
>>> 1. - psturng(qs, [4, 4, 4], [6, 6, 6])
array([ 0.89320111,  0.93449991,  0.98269855])
>>> 
```

---

## Requirements ##
  * Python (tested with version 2.7)
  * SciPy (tested with version 0.9.0)
  * Numpy (tested with version 1.5.1)

---

## Installation ##

> Download [qsturng.py](http://code.google.com/p/qsturng-py/downloads/detail?name=qsturng.py&can=2&q=) and place in site-packages folder.

---

## Mathematics behind the approximation ##
The studentized range distribution _q_, and the student distribution _t_ share a great deal of resemblance with one another. In fact, when _r_=2 the following relationship holds,
> ![https://chart.googleapis.com/chart?cht=tx&chl=q_p(2%2Cv)%3D%5Csqrt{2}t_{p'}(v)%5C%3B%5Cmathrm{where}%5C%3Bp'%3D%5Cfrac{1%2Bp}{2}%5C%3B%5Cmathrm{for%5C%3Bany}%5C%3B0%3Cp%5C%3B%5Cmathrm{and%5C%3Ball}%5C%3Bv%5Cgeq1&.png](https://chart.googleapis.com/chart?cht=tx&chl=q_p(2%2Cv)%3D%5Csqrt{2}t_{p'}(v)%5C%3B%5Cmathrm{where}%5C%3Bp'%3D%5Cfrac{1%2Bp}{2}%5C%3B%5Cmathrm{for%5C%3Bany}%5C%3B0%3Cp%5C%3B%5Cmathrm{and%5C%3Ball}%5C%3Bv%5Cgeq1&.png).
Because _q_ and _t_ have such a strong relation we can use more commonly available _t_ inverse functions to find _q_ if we can account for the variability of _r_. This can be accomplished by forming the following variable,
> ![https://chart.googleapis.com/chart?cht=tx&chl=y_p(r%2Cv)%3Dq_p(r%2Cv)%2F%5Csqrt%7B2%7Dt_%7Bp'%7D(v)&.png](https://chart.googleapis.com/chart?cht=tx&chl=y_p(r%2Cv)%3Dq_p(r%2Cv)%2F%5Csqrt%7B2%7Dt_%7Bp'%7D(v)&.png).

Gleason (1998, 1999) showed that this is possible by applying ordinary least squares (OLS) regression to r at tabled values of p and v such that ![https://chart.googleapis.com/chart?cht=tx&chl=y_p(r%2Cv)&.png](https://chart.googleapis.com/chart?cht=tx&chl=y_p(r%2Cv)&.png) is approximated by,
> ![https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7By%7D_%7Bp%2Cv%7D(r)%3D%201%20%2B%20%5Calpha_1%20%5Clog%7B(r-1)%7D%20%2B%20%5Calpha_2%20%5Clog%5E2%7B(r-1)%7D%20%2B%20%5Calpha_3%20%5Clog%5E3%7B(r-1)%7D%20%2B%20%5Calpha_4%20%5Clog%5E4%7B(r-1)%7D&.png](https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7By%7D_%7Bp%2Cv%7D(r)%3D%201%20%2B%20%5Calpha_1%20%5Clog%7B(r-1)%7D%20%2B%20%5Calpha_2%20%5Clog%5E2%7B(r-1)%7D%20%2B%20%5Calpha_3%20%5Clog%5E3%7B(r-1)%7D%20%2B%20%5Calpha_4%20%5Clog%5E4%7B(r-1)%7D&.png).

Gleason's original algorithm (implemented in Stata) applied least-squares regression to 8 levels of _p_: .5, .75, .9, .95, .975, .99, .99, & .995, and at each probability 25 ( when _p_<.9) or 26 (when p>=.9) levels of v between 1 or 2 respectively and infinity. Here three additional levels of p have been added at .675, .8, and .85 using _q_ values obtained by R's _qtukey_ function. This reduces RMSD by 60% for interpolated estimates when _p_ is between .5 and .9 (see the following section). At each combination of _p_ and _v_ OLS is applied to 24 levels of _r_ between 2 and 100.

Approximations of ![https://chart.googleapis.com/chart?cht=tx&chl=q_p(r,v)&.png](https://chart.googleapis.com/chart?cht=tx&chl=q_p(r,v)&.png) then take the form of,
> ![https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7Bq_p%7D(r%2Cv)%3D%5Csqrt%7B2%7D(1%2B%5Chat%7Bf%7D)t_%7Bp'%7D&.png](https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7Bq_p%7D(r%2Cv)%3D%5Csqrt%7B2%7D(1%2B%5Chat%7Bf%7D)t_%7Bp'%7D&.png)<br>
<blockquote>where<br>
<blockquote><img src='https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7Bf%7D(p%2Cr%2Cv)%3D%5Csum_%7Bj%3D1%7D%5E%7B4%7Da_j%20%5Clog%5Ej(r-1)%2Bb(p)%2Bc(v)&.png' /><br>
where<br>
<blockquote><img src='https://chart.googleapis.com/chart?cht=tx&chl=b(p)%3D%5Cbegin%7Bcases%7D%20-0.002%2F(1%2B12%7Bz%5E2%7D_p)%2C%20%26%20%5Cmbox%7Bif%20%7D%20r%3D3%20%5C%5C%20%0A%5C%3B%20%5C%3B%20%5C%3B%200%2C%20%26%20%5Cmbox%7Botherwise%7D%20%5Cend%7Bcases%7D&.png' /><br><br>
<img src='https://chart.googleapis.com/chart?cht=tx&chl=c(v)%3D%5Cbegin%7Bcases%7D%201%2F517-1%2F(312v))%2C%20%26%20%5Cmbox%7Bif%20%7D%20r%3D3%20%5Cmbox%7B%20and%20%7D%20v%20%5Cleq%204.364%5C%5C%201%2F(191v)%2C%20%26%20%5Cmbox%7Bif%20%7D%20r%3D3%20%5Cmbox%7B%20and%20%7D%20v%20%3E%204.364%5C%5C%20%0A0%2C%20%26%20%5Cmbox%7Botherwise%7D%20%5Cend%7Bcases%7D&.png' /><br>
where<br>
<blockquote><img src='https://chart.googleapis.com/chart?cht=tx&chl=z_p%3D%5CPhi%5E%7B-1%7D(p)&.png' /> (pth quantile inverse norm)</blockquote></blockquote></blockquote></blockquote>

As the reader can infer,<br>
<blockquote><img src='https://chart.googleapis.com/chart?cht=tx&chl=y_p(r%2Cv)%5Capprox%20%5Chat%7By%7D_%7Bp%2Cv%7D(r)%3D%5Chat%7Bf%7D(p%2Cr%2Cv)%2B1&.png' />.</blockquote>

The functions <img src='https://chart.googleapis.com/chart?cht=tx&chl=b(p)&.png' /> and <img src='https://chart.googleapis.com/chart?cht=tx&chl=c(p)&.png' /> provide some minor adjustments when <img src='https://chart.googleapis.com/chart?cht=tx&chl=r%3D3&.png' />.<br>
<br>
This formulation allows us to find <i>q</i> at tabled values of <i>p</i> and <i>v</i> for any <i>r</i>, but it still doesn't allow us to find values of <i>q</i> at non-tabled values of  <i>p</i> an <i>v</i>. This requires interpolating <img src='https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7By%7D_%7Bp%2Cv%7D(r)&.png' />.<br>
<br>
<img src='http://qsturng-py.googlecode.com/files/error.png' />
<blockquote>Figure 2. <i>Histogram depicts normalized error in calculating the tabled values across<br>
the 7344 combinations of p, r, and v.</i>
<hr />
<h2>Mathematics behind the interpolations</h2>
<h3>Interpolating <i>r</i></h3>
Part of the beauty of the OLS regression is it automatically takes care of interpolation over <i>r</i>. The table values for <i>r</i> only extend to 100, but interpolation is stable to 200 with only a slight loss of accuracy.</blockquote>

<blockquote><h3>Interpolating <i>v</i></h3>
Interpolation of <img src='https://chart.googleapis.com/chart?cht=tx&chl=v&.png' /> requires transforming the ordinate to <img src='https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7By%7D%5E2&.png' /> and transforming the abscissa to <img src='https://chart.googleapis.com/chart?cht=tx&chl=1%2Fv&.png' />. These transforms effectively linearize the space between <img src='https://chart.googleapis.com/chart?cht=tx&chl=%5Chat%7By%7D&.png' /> and <img src='https://chart.googleapis.com/chart?cht=tx&chl=v&.png' /> values. For good measure 3 point quadratic interpolation is used to approximate between tabled values of <img src='https://chart.googleapis.com/chart?cht=tx&chl=v&.png' />.</blockquote>

<blockquote><h3>Interpolating <i>p</i></h3>
Interpolation of <img src='https://chart.googleapis.com/chart?cht=tx&chl=p&.png' /> is slightly more complicated. For all interpolations the ordinate is transformed by <img src='https://chart.googleapis.com/chart?cht=tx&chl=%5Clog%7B(%5Chat%7By%7D-r%2Fv)%7D&.png' />. For values of <img src='https://chart.googleapis.com/chart?cht=tx&chl=p%3E.85&.png' /> the abcissa is transformed by:<br>
<blockquote><img src="https://chart.googleapis.com/chart?cht=tx&chl=x(p)%3D-1%2F(1%2B1.5z_%7Bp'%7D)%20%5C%3B&.png" />
<blockquote>where <img src="https://chart.googleapis.com/chart?cht=tx&chl=z_%7Bp'%7D%3D%5CPhi%5E%7B-1%7D((1%2Bp)%2F2)&.png" />
<blockquote>where <img src='https://chart.googleapis.com/chart?cht=tx&chl=\Phi^{-1}(p)&.png' /> is the inverse of the normal distribution of <img src='https://chart.googleapis.com/chart?cht=tx&chl=p&.png' /><br>
</blockquote></blockquote></blockquote>From the plot above it is clear that below .85 <i>p</i> is fairly well behaved and the <i>x(p)</i> transform is not required. Gleason's implementation actually uses the transform, but I have found the RMSD values against R's <i>qtukey</i> function are slightly lower when the abcissa transformation is not applied below .85. When .5 < <i>p</i> < .85 only the ordinate transform is applied. When .1 < <i>p</i> < .5 linear interpolation in the <i>p</i>/<i>q</i> space seems to yield the best results. Interpolation below .5 is not recommended. It is implemented primarily to give rough approximations for when <i>psturng</i> is trying to find values above .5.</blockquote>

<blockquote><h3>Interpolating <i>p</i> and <i>v</i></h3>
When <i>p</i> and <i>v</i> must both be interpolated a 9 point hyperbolic interpolation scheme is applied (quadratic interpolation treating <i>p</i> and <i>v</i> as orthogonal.)</blockquote>

<img src='http://qsturng-py.googlecode.com/files/non-tbl_error.png' />
<blockquote>Figure 3. <i>Histogram depicts normalized "error" in approximations compared to the Copenhaver & Holland (1988) algorithm. The bias is expected since qsturng is fit mostly with Harter's 1960 table values. Gleason notes that about 20% of the table values are conservatively rounded. The primary point of the figure is to illustrate the stability of the interpolations. The histogram reflects 10,000 randomly selected non-tabled values where p, r, and v were chosen at random. qtukey (in R) tends to fail to converge in certain ranges of p, r, and v. To obtain a usable metric the ranges of p, r, and v were restricted to regions where qtukey is more stable and NaN and unreasonably large values were filtered before normalized errors were computed. p was restricted to .5 to .95, r was restricted to integers between 2 and 100, v was restricted to reals between 2 and 1000. As a side note, qtukey estimates at p = .999 often differed from Harter's table by more than 20% when degrees of freedom were low.</i></blockquote>

<b>Here our intention is introduce the basic concept behind the approximation, and point out implementation distinctions between the implemented algorithm and the algorithm as described by Gleason. For more details please see Gleason, J. R. (1999), or the source code.</b>
<hr />
<h2><i>psturng</i> algorithm</h2>
The <i>psturng</i> algorithm uses a scalar minization (<i>scipy.optimize.fminbound</i>) to calculate <i>p</i> values given a quantile, <i>r</i> value, and <i>v</i> value. To make sure the algorithm is not falling into local minima 10,000 combinations of <i>p</i> between .1 to .999, <i>r</i> between 2 to 1000, and <i>v</i> between 2 to 1000 were chosen at random. For each value <i>qsturng</i> was used to calculate the corresponding quantile and then <i>psturng</i> was used to find <i>p</i>. For all 10,000 values <i>psturng</i> was able to approximate the original <i>p</i> within <i>fminbound</i>'s tolerance of 1e-5.<br>
<blockquote><img src='http://qsturng-py.googlecode.com/files/psturng_error.png' />
Figure 4. <i>Histogram depicts approximation errors (estimate - actual) of psturng. The uniform distribution reflects the tolerance parameter of the scalar minimization algorithm used by psturng.</i></blockquote>

<h2>References</h2>
<blockquote>Copenhaver, M. D., & Holland, B. S. (1988) Multiple comparisons of simple effects in the<br>
<blockquote>two-way analysis of variance with fixed effects. Journal of Statistical Computation<br>
and Simulation, (30), 1–15.</blockquote></blockquote>

<blockquote>Gleason, J. R. (1999). An accurate, non-iterative approximation<br>
<blockquote>for studentized range quantiles. Computational Statistics &<br>
Data Analysis, (31), 147-158.</blockquote></blockquote>

<blockquote>Gleason, J. R. (1998). A table of quantile points of the<br>
<blockquote>Studentized range distribution.<br>
<a href='http://www.stata.com/stb/stb46/dm64/sturng.pdf'>http://www.stata.com/stb/stb46/dm64/sturng.pdf</a>