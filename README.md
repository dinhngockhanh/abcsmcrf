#   Approximate Bayesian Computation sequential Monte Carlo via random forests

##  Installation

The ABC-SMC-(D)RF library can be installed with

```R
devtools::install_github("dinhngockhanh/abc-smc-rf")
```

##  Vignettes

The `vignettes` folder contains examples of using ABC-SMC-(D)RF to infer parameters in different mathematical models. These examples are showcased in the paper "Approximate Bayesian Computation sequential Monte Carlo via random forests" [1].

Files `Coalescent model with stats=S.r`, `Coalescent model with stats=S+SFS.r`, and `Coalescent model with stats=SFS.r` examine the inference of mutation rate in a coalescent model, with different statistics from the observed DNA sequences (Figure 1; Section 2.2; [1]) using ABC-RF [2].

File `Hierarchical model.r` uses ABC-DRF [3] to infer parameters in a hierarchical normal mean model (Figure 2; Section 3.1; [1]).

File `Lotka-Volterra model.r` uses ABC-SMC-(D)RF to infer the birth rates of prey and predators in the deterministic Lotka-Volterra model (Figure 3; Section 5.1; [1]).

File `Birth-death model.r` uses ABC-SMC-(D)RF to infer the birth and death rates in a linear birth-death branching process (Figure 4; Section 5.2; [1]).

File `Michaelis-Menten model.r`uses ABC-SMC-(D)RF to infer reaction rates in the Michaelis-Menten model (Figure 5; Section 5.3; [1]).

##  References
1.  Dinh KN, Xiang Z, Liu Z, Tavaré S. Approximate Bayesian Computation sequential Monte Carlo via random forests. arXiv preprint arXiv:2406.15865. 2024 Jun 22.
2.  Raynal L, Marin JM, Pudlo P, Ribatet M, Robert CP, Estoup A. ABC random forests for Bayesian parameter inference. Bioinformatics. 2019 May 15;35(10):1720-8.
3.  Cevid D, Michel L, Näf J, Bühlmann P, Meinshausen N. Distributional random forests: Heterogeneity adjustment and multivariate distributional regression. Journal of Machine Learning Research. 2022;23(333):1-79.