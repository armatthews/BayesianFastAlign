# BayesianFastAlign
A Bayesian implementation of cdec's fast_align

Dependencies:
- Boost
- https://github.com/redpony/cpyp

Note: Do not try to compile with -Ofast.
cpyp has some regrettable (but intentional) NaN usage, and enabling fast-math disables NaN checks that it depends on.
