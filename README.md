# nps-models
Some code for calculating uncertainty in [Net Promoter Score](https://en.wikipedia.org/wiki/Net_Promoter) data based on something I did for an interview exercise. Uses R and [Stan](https://mc-stan.org). Should work with any response scale, and can, with a little modification, be used to calculate other scores of this type that don't throw away information like NPS.

## Contents:

* NPS-multinomial.stan: Stan code modeling response proportions as a multinomial distribution. Makes some incorrect assumptions (specifically that differences like 1 vs 2 are not distinct from differences like 1 vs 10), but faster and less prone to sampling problems.
* NPS-multinomial-groups.stan: Same as above, but estimates proportions conditional on some group membership predictor
* NPS-cumulative.stan: Stan code modeling responses using an ordered/cumulative logit model. Better assumption-wise and less of a hassle to incorporate non-categorical predictor variables. Throws some (usually harmless) warnings about rejected metropolis proposals.
* NPS-models.R: R code for running models. Will require some editing to use.

## Notes:

All of these models incorporate a Dirichlet prior on the counts in each response bin. This is equivalent to adding some small number of pseudo-observations to each count, and is useful for smoothing when counts are low or 0 in some set of bins.
