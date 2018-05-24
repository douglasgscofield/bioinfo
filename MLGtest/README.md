MLGtest
=======

MLGtest implements a test for asexuality based on the number of observed
multilocus genotypes (MLGs).  For each of a given number of iterations (default
1000), _n_ multilocus genotypes are constructed based on given allele frequency
distributions.  The number of unique MLGs and the MLG frequency count are
recorded each iteration.  These are used to construct two distributions
specific to a sample of size _n_: a cumulative distribution of the number of
unique MLGs expected; and the mean frequency counts of unique MLGs.  Both of
these distributions can be returned to the researcher, and the MLG test will
test an observed number of MLGs against the expected number of MLGs.  MLGtest
has been used by [Werth and Sork (<i>American Journal of Botany</i>
2008)][Werth2008] and by [Allen and Lynch (<i>Evolution</i> 2012)][Allen2012].

The problem of detecting asexuality via multilocus genotype frequencies has
been addressed at considerable depth by [Parks and Werth (1993)][Parks1993],
[Ivey and Richards (2001)][Ivey2001] and [Stenberg, Lundmark and Saura
(2003)][Stenberg2003].  MLGtest turns out to be largely a duplication of effort
by Stenberg, Lundmark and Saura, but lacking the considerable further advances
provided in those authors' `MLGsim`.  [Arnaud-Haond et al.
(2007)][ArnaudHaond2007] address some additional issues and give advice on
adjusting allele frequencies.

The major advance in these other treatments not included here is that other
methods consider the probability of the specific genotypes observed, not just
the probability of a given distribution of genotype counts and number of
genotypes observed.  As a result, MLGtest turns out to be quite conservative
and would not be recommended if a sensitive test is desired.  In the
applications above, the taxa involved (_Ramalina_ lichen and _Daphnia_
crustaceans) readily reproduce asexually.


[Werth2008]: https://doi.org/10.3732/ajb.2007024
[Allen2012]: https://doi.org/10.1111/j.1558-5646.2011.01488.x
[Parks1993]: http://www.jstor.org/stable/2445369
[Ivey2001]: https://doi.org/10.1086/320775
[Stenberg2003]: https://doi.org/10.1046/j.1471-8286.2003.00408.x
[ArnaudHaond2007]: https://doi.org/10.1111/j.1365-294X.2007.03535.x
