library(tidyverse)
NumSNPs <- 100000
NumSam <- 50 # number of haploid chromosomes

# here are the unnormalized fractions we expect to find
fracts <- 1 / 1:(NumSam - 1)

# normalize so they sum to one, the multiply by NumSNPs
Sin <- NumSNPs * fracts / sum(fracts)

# let's plot that:
tibble(
    i = 1:(NumSam - 1),
    E_of_S = Sin
) %>%
    ggplot(aes(x = i, y = E_of_S)) +
    geom_col() +
    theme_bw()
