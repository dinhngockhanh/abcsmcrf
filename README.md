#   Approximate Bayesian Computation sequential Monte Carlo via random forests

##  Installation

The ABC-SMC-(D)RF library can be installed with

```R
devtools::install_github("dinhngockhanh/abc-smc-rf")
```

##  Vignettes

The `vignettes` folder contains examples of using ABC-SMC-(D)RF to infer parameters in different mathematical models. These examples are showcased in the paper "Approximate Bayesian Computation sequential Monte Carlo via random forests".

Files `Coalescent model with stats=S.r`, `Coalescent model with stats=S+SFS.r`, and `Coalescent model with stats=SFS.r` examine the inference of mutation rate in a coalescent model, with different statistics from the observed DNA sequences (Figure 1; Section 2.2) using ABC-RF .

File `Hierarchical model.r` studies 

Fig 2, 3.1

Fig 3, 5.1
`Lotka-Volterra model.r`

Fig 4, 5.2
`Birth-death model.r`

Fig 5, 5.3
`Michaelis-Menten model.r`

##  References
1.  Dinh KN, Xiang Z, Liu Z, Tavaré S. Approximate Bayesian Computation sequential Monte Carlo via random forests. arXiv preprint arXiv:2406.15865. 2024 Jun 22.
2.  Raynal L, Marin JM, Pudlo P, Ribatet M, Robert CP, Estoup A. ABC random forests for Bayesian parameter inference. Bioinformatics. 2019 May 15;35(10):1720-8.
3.  Cevid D, Michel L, Näf J, Bühlmann P, Meinshausen N. Distributional random forests: Heterogeneity adjustment and multivariate distributional regression. Journal of Machine Learning Research. 2022;23(333):1-79.