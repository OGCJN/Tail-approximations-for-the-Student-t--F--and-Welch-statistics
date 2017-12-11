# Tail approximations for the Student t−, F−, and Welch statistics for non-normal and not necessarily i.i.d. random variables

## Abstract
We present a detailed study of the asymptotic behavior of the distribution of the tails of these, perhaps, most commonly used statistical tests under non-standard conditions, that is, releasing the underlying assumptions of normality, independence and identical distribution and considering a more general case where one only assumes that the vector of data has a continuous joint density. We determine asymptotic expressions for P(T > u) as u tends to infinity for this case. The approximations are particularly accurate for small sample sizes and may be used, for example, in the analysis of High-Throughput Screening experiments, where the number of replicates can be as low as two to five and often extremely high significance levels are used. We give numerous examples and complement our results by a thorough investigation of the convergence speed - both theoretically, by deriving exact bounds for absolute and relative errors of the approximations, and by means of a simulation study.

## Supplementary materials

MATLAB

[OST /TST /WELCH /F ]+ComputeKg.m - compute K<sub>g</sub> for the Student one- and two- sample t−, Welch, and F− statistics using adaptive Simpson or Lobatto quadratures. Here g is an arbitrary multivariate density.<sup>1</sup>

[TST /WELCH /F ]+ComputeKgIS+.m - the same as above but for the case where samples are independent.<sup>2</sup>

[OST /TST /WELCH /F ]+ComputeKgIID+.m - the same as above but assuming that the samples consist of i.i.d. random variables.<sup>2</sup>

RunSimulation+[IID/MVN ]+.m - perform simulation study for i.i.d. and dependent/non-homogeneous cases, see Section 7 and Appendix B.

Wolfram Mathematica

[OST /TST /WELCH /F ]+ComputeKg.nb - compute the exact expression for K<sub>g</sub> for an arbitrary multivariate density g and given sample size(s). We include a number of examples, such as evaluation of K<sub>g</sub> for the zero-mean Gaussian case with an arbitrary covariance matrix **Σ**; the “unequal variances” case for the Student two-sample t− and Welch statistics; and evaluation of K<sub>g</sub>. for the densities considered in the simulation study.

OSTComputeKgIID.nb - veriﬁes the constants in Table 1 for the i.i.d. case of the Student one-sample t−statistic.

TSTExactPDF.nb and WELCHExactPDF.nb - the exact distribution for the Student two-sample t− and Welch statistics for odd sample sizes, see (Ray and Pitman, 1961).

Other Materials

Supplementary-Materials.pdf - Remarks on Theorem 1.1 and its application to real data; extended version of the literature review; comparison of the result of Theorem 1.1 with the exact distribution of the Welch statistic; proof of Theorem 5.1.
 
___
<sup>1</sup>For the F−statistic we use Monte Carlo integration. 
<sup>2</sup>For the F−statistic and n<sub>1</sub> > 3 we use Monte Carlo integration.

## Reference
Zholud, D. (2014). [**Tail approximations for the Student t−, F−, and Welch statistics for non-normal and not necessarily i.i.d. random variables**](http://www.zholud.com/articles/Tail-approximations-for-the-Student-t-,-F-,-and-Welch-statistics-for-non-normal-and-not-necessarily-i.i.d.-random-variables.pdf), **Bernoulli**, Vol. 20, No. 4, pp. 2102-2130

## BiBTeX

``` BiBTeX
@article{Zholud2014,
  Author = {Zholud, D.},
  Year = {2014},
  Title = {Tail approximations for the Student t−, F−, and Welch statistics for non-normal and not necessarily i.i.d. random variables},
  Journal = {Bernoulli},
  Volume = {20},
  Number = {4},
  Pages = {2102--2130}
}
``` 
