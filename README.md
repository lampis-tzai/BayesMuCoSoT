# BayesMuCoSoT
R library for Bayesian multivariate common source testing

Bayesian multivariate normal modelling (family Normal-inverse-Wishart) for common source testing. This library can help you to test if two multivariate distribution comes from the same source based on logarithmic Bayes factor. 

## Dependence

Bayesian Software JAGS (Just Another Gibbs Sampler). It is a program for the statistical analysis of Bayesian hierarchical models by Markov Chain Monte Carlo.
Can be downloaded from [here](https://sourceforge.net/projects/mcmc-jags/) 

## Documentation

[Theory]()

[Library Site](https://lampis-tzai.github.io/BayesMuCoSoT/)

[Examples Markdown]()

## Installation
```
#install.packages("devtools")
devtools::install_github('lampis-tzai/BayesMuCoSoT')
```
## Usage

```
library(BayesMuCoSoT)

all_data = iris[iris$Species=='setosa',]

#create questioned data
questioned_data = all_data[1:(nrow(all_data)/2),]

#create known data
known_data = all_data[(nrow(all_data)/2+1):nrow(all_data),]

#create background data
background_data = iris[iris$Species!='setosa',]

#test
background_data_id = 'Species'
y = names(questioned_data)[1:4]
BayesMuCoSoT_fit(y, x=NA, questioned_data, known_data, background_data,background_data_id)
```


## References

Bozza, Taroni, Marquis, Schmittbuhl, “Probabilistic Evaluation of Handwriting Evidence: Likelihood Ratio for Authorship.” Journal of the Royal Statistical Society: Series C (Applied Statistics) 57, no. 3 (June 1, 2008): 329–41. 
