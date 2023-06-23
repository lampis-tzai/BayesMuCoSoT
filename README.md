# BayesMuCoSoT
R library for Bayesian multivariate common source testing

Bayesian multivariate normal modelling (family Normal-inverse-Wishart) for common source testing. This library can help you to test if two multivariate distribution comes from the same source based on logarithmic Bayes factor. 


## Documentation

[Theory]()

[Library Site](https://lampis-tzai.github.io/BayesMuCoSoT/)

[Kaggle Dataset Example Markdown]()

## Installation
```
#install.packages("devtools")
devtools::install_github('lampis-tzai/BayesMuCoSoT)
```
## Usage

#### Common Source testing Iris dataset

```
library(BayesMuCoSoT)
all_data = iris[iris$Species=='setosa',]
questioned_data = all_data[1:(nrow(all_data)/2),]
known_data = all_data[(nrow(all_data)/2+1):nrow(all_data),]
background_data = iris[iris$Species!='setosa',]
background_data_id = 'Species'
y = names(questioned_data)[1:4]
BayesMuCoSoT_fit(y, x=NA, questioned_data, known_data, background_data,background_data_id)
```

#### Different Source testing Iris dataset

```
library(BayesMuCoSoT)
set.seed(10)
questioned_data = iris[iris$Species=='versicolor',]
split_ind <- sample(seq_len(nrow(questioned_data)), size = floor(0.5 * nrow(questioned_data)))
background_data1 = questioned_data[-split_ind,]
questioned_data = questioned_data[split_ind,]


known_data = iris[iris$Species=='virginica',]
split_ind <- sample(seq_len(nrow(known_data)), size = floor(0.5 * nrow(known_data)))
background_data2 = known_data[-split_ind,]
known_data = known_data[split_ind,]

background_data = rbind(iris[iris$Species=='setosa',],
                        background_data1,
                        background_data2)

background_data_id = 'Species'
y = names(questioned_data)[1:4]
BayesMuCoSoT_fit(y,x=NA,questioned_data,known_data,
                              background_data,background_data_id)
```

#### Different Source testing Simulated dataset

```
library(BayesMuCoSoT)
library(LaplacesDemon) #for sampling for multivariate normal
set.seed(10)
questioned_data = as.data.frame(rmvn(10,rep(0,2),array(c(2,0.5,0.5,2),dim = c(2,2))))
known_data = as.data.frame(rmvn(8,rep(2,2),array(c(2,0.5,0.5,2),dim = c(2,2))))
background_data = as.data.frame(rbind(rmvn(10,rep(0.5,2),array(c(1,0.2,0.2,2),dim = c(2,2))),
                       rmvn(10,rep(1.5,2),array(c(6,0.3,0.3,4),dim = c(2,2))),
                       rmvn(10,rep(3,2),array(c(3,0.8,0.8,3),dim = c(2,2)))))
background_data['id'] = rep(1:3,each=10)
background_data_id = 'id'
y = names(questioned_data)[1:2]
BayesMuCoSoT_fit(y,x=NA,questioned_data,known_data,background_data,background_data_id)
```


## References

Bozza, Taroni, Marquis, Schmittbuhl, “Probabilistic Evaluation of Handwriting Evidence: Likelihood Ratio for Authorship.” Journal of the Royal Statistical Society: Series C (Applied Statistics) 57, no. 3 (June 1, 2008): 329–41. 
