#' Conjugate Marginal Likelihood
#'
#' @param jags_data list of the names of the data objects used by the model (like jags data).
#' @return the logarithmic marginal likelihood
#' @import CholWishart
#' @export
marginal_likelihood_conjugate<- function(jags_data){

  y = jags_data$y
  x = jags_data$x_pred
  m0 = jags_data$beta_mu
  v0 = jags_data$nw
  k0 = diag(0.5,jags_data$L,jags_data$L)
  U0 = jags_data$U


  n = nrow(x)
  d = ncol(y)
  vn = v0 + n

  kn = t(x) %*% x + k0
  mn = solve(kn)%*%(t(x)%*%y + k0%*%m0)
  Un = U0 + t(y) %*% y + t(m0)%*%k0%*%m0 - t(mn)%*%kn%*%mn


  logml = - ((d*n)/2)*log(2*pi) + (d/2)*determinant(k0,logarithm = TRUE)$modulus[1] -
    (d/2)*determinant(kn,logarithm = TRUE)$modulus[1] +
    (v0/2)*determinant(U0/2,logarithm = TRUE)$modulus[1] -
    (vn/2)*determinant(Un/2,logarithm = TRUE)$modulus[1] +
    CholWishart::lmvgamma(vn/2, d) - CholWishart::lmvgamma(v0/2, d)

  return(logml)

}

#' Laplace-Metropolis  Marginal Likelihood Approximation
#'
#' @param samples samples of jags output (output$BUGSoutput$sims.matrix).
#' @param jags_data list of the names of the data objects used by the model (like jags data).
#' @return the logarithmic marginal likelihood
#' @import CholWishart LaplacesDemon
#' @export
marginal_likelihood_laplace_metr <- function(samples, jags_data){

  samples = as.data.frame(samples)
  samples = samples[names(samples)!='deviance']


  beta_mu_hat <- array(0, dim=c(jags_data$L,jags_data$P))
  for (r in 1:jags_data$L){
    beta_mu_hat[r,]=unlist(colMeans(samples[paste0("beta[", r
                                                   ,",",seq(jags_data$P), "]")]))
  }


  W_hat_post <- array(0, dim=c(jags_data$P,jags_data$P))
  for (r in seq(jags_data$P)){
    for (c in seq(jags_data$P)){
      W_hat_post[r,c] = mean(samples[paste0("W[",r,",",c,"]")][,1])
    }}


  comp_var = t(combn(seq(jags_data$P),2))
  new_col_names = colnames(samples)

  for (v in 1:nrow(comp_var)){
    new_col_names = new_col_names[new_col_names!=paste0("W[",comp_var[v,1],",",comp_var[v,2],"]")]
  }
  samples_new = samples[new_col_names]

  #samples_new = samples[!duplicated(as.list(samples))]

  par_cor = cor(samples_new)

  likelihood = sum(LaplacesDemon::dmvn(jags_data$y,  jags_data$x_pred %*% beta_mu_hat ,W_hat_post, log = TRUE))

  prior_betas=c()
  for (l in 1:jags_data$L){
    prior_beta_prob = LaplacesDemon::dmvn(beta_mu_hat[l,], jags_data$beta_mu[l,],
                           jags_data$beta_cov[,,l], log = TRUE)
    prior_betas = c(prior_betas,prior_beta_prob)
  }
  prior_beta = sum(prior_betas)

  prior_W = CholWishart::dInvWishart(W_hat_post,jags_data$nw, jags_data$U, log = TRUE)


  logml = (ncol(samples_new)/2)*log(2*pi) + (1/2)*determinant(par_cor,logarithm = TRUE)$modulus[1] +
    sum(log(sapply(samples_new, sd))) +
    likelihood + prior_beta + prior_W


  return (logml)
}

#' Generalized Harmonic Mean  Marginal Likelihood Approximation
#'
#' @param samples samples of jags output (output$BUGSoutput$sims.matrix).
#' @param jags_data list of the names of the data objects used by the model (like jags data).
#' @return the logarithmic marginal likelihood
#' @import abind CholWishart LaplacesDemon
#' @export
ml_generalized_harmonic_mean <- function(samples, jags_data) {

  samples = as.data.frame(samples)

  smp_size <- floor(0.5 * nrow(samples))

  iter_fit_ind <- sample(seq_len(nrow(samples)), size = smp_size)

  samples_fit <- samples[iter_fit_ind, ]
  samples_iter <- samples[-iter_fit_ind, ]

  N1 = nrow(samples_fit)
  N2 = nrow(samples_iter)

  beta_mu_hat <- array(0, dim=c(jags_data$L,jags_data$P))
  beta_cov_hat <- list()
  for (r in 1:jags_data$L){
    beta_mu_hat[r,]=unlist(colMeans(samples_fit[paste0("beta[", r
                                                       ,",",seq(jags_data$P), "]")]))
    beta_cov_hat[[r]]= cov(unname(samples_fit[paste0("beta[", r
                                                     ,",",seq(jags_data$P), "]")]))
  }

  W_hat_post <- array(0, dim=c(jags_data$P,jags_data$P))
  for (r in seq(jags_data$P)){
    for (c in seq(jags_data$P)){
      W_hat_post[r,c] = mean(samples_fit[paste0("W[",r,",",c,"]")][,1])
    }}

  nw_post =  jags_data$N + jags_data$nw - jags_data$L
  U_post =  W_hat_post * (nw_post - jags_data$P  - 1)


  beta_df <- samples_iter[paste0("beta[", 1 ,",",seq(jags_data$P), "]")]
  if (jags_data$L>1){
  for (r in 2:jags_data$L){
    beta_df =cbind(beta_df, samples_iter[paste0("beta[", r ,",",seq(jags_data$P), "]")])
  }
  }


  beta = abind::abind(apply(beta_df,1,
                     function(x) t(matrix(x,nrow=jags_data$P,ncol=jags_data$L)),
                     simplify = FALSE), along=3)

  W = abind::abind(apply(data.matrix(samples_iter[,1:(jags_data$P*jags_data$P)]),1,
                  function(x) matrix(x,nrow=jags_data$P,ncol=jags_data$P),
                  simplify = FALSE), along=3)


  prior_betas = array(0, dim = c(N2,jags_data$L))
  for (l in 1:jags_data$L){
    prior_beta_prob = LaplacesDemon::dmvn(t(beta[l,,]),
                           jags_data$beta_mu[l,],
                           jags_data$beta_cov[,,l], log = TRUE)
    prior_betas[,l] = prior_beta_prob
  }
  prior_beta = rowSums(prior_betas)

  prior_W = CholWishart::dInvWishart(W,jags_data$nw, jags_data$U, log = TRUE)

  posterior_betas = array(0, dim = c(N2,jags_data$L))
  for (l in 1:jags_data$L){
    posterior_beta_prob = LaplacesDemon::dmvn(t(beta[l,,]),
                               beta_mu_hat[l,],
                               beta_cov_hat[[l]],
                               log = TRUE)

    posterior_betas[,l] = posterior_beta_prob
  }
  posterior_beta = rowSums(posterior_betas)

  posterior_W = CholWishart::dInvWishart(W,nw_post, U_post, log = TRUE)

  prob_list = numeric(nrow(samples_iter))
  for (i in 1:nrow(samples_iter)){

    likelihood = sum(LaplacesDemon::dmvn(jags_data$y, jags_data$x_pred %*% beta[,,i], W[,,i], log = TRUE))

    prob =  posterior_beta[i] + posterior_W[i] - likelihood -  prior_beta[i] -prior_W[i]

    prob_list[i] = prob
  }

  prob_star = median(prob_list)

  logml = log(mean(exp(prob_list - prob_star)))+prob_star

  return (-logml)

}


#' Bridge Sampling  Marginal Likelihood Approximation
#'
#' @param samples samples of jags output (output$BUGSoutput$sims.matrix).
#' @param jags_data list of the names of the data objects used by the model (like jags data).
#' @return the logarithmic marginal likelihood
#' @import abind CholWishart LaplacesDemon
#' @export
ml_bridge_sampling <- function(samples, jags_data) {

  samples = as.data.frame(samples)

  smp_size <- floor(0.5 * nrow(samples))

  iter_fit_ind <- sample(seq_len(nrow(samples)), size = smp_size)

  samples_fit <- samples[iter_fit_ind, ]
  samples_iter <- samples[-iter_fit_ind, ]

  N1 = nrow(samples_fit)
  N2 = nrow(samples_iter)


  beta_mu_hat <- array(0, dim=c(jags_data$L,jags_data$P))
  beta_cov_hat <- list()
  for (r in 1:jags_data$L){
    beta_mu_hat[r,]=unlist(colMeans(samples_fit[paste0("beta[", r
                                                       ,",",seq(jags_data$P), "]")]))
    beta_cov_hat[[r]]= cov(unname(samples_fit[paste0("beta[", r
                                                     ,",",seq(jags_data$P), "]")]))
  }

  W_hat_post <- array(0, dim=c(jags_data$P,jags_data$P))
  for (r in seq(jags_data$P)){
    for (c in seq(jags_data$P)){
      W_hat_post[r,c] = mean(samples_fit[paste0("W[",r,",",c,"]")][,1])
    }}

  nw_post =  jags_data$N + jags_data$nw - jags_data$L
  U_post =  W_hat_post * (nw_post - jags_data$P  - 1)

  beta <- array(0, dim=c(jags_data$L,jags_data$P,N1))
  for (l in 1:jags_data$L){
    beta[l,,] = t(LaplacesDemon::rmvn(N1,beta_mu_hat[l,],beta_cov_hat[[l]]))
  }

  W = CholWishart::rInvWishart(N1,nw_post, U_post)


  prior_betas = array(0, dim = c(N1,jags_data$L))
  for (l in 1:jags_data$L){
    prior_beta_prob = LaplacesDemon::dmvn(t(beta[l,,]),
                           jags_data$beta_mu[l,],
                           jags_data$beta_cov[,,l], log = TRUE)
    prior_betas[,l] = prior_beta_prob
  }

  prior_beta = rowSums(prior_betas)

  prior_W = CholWishart::dInvWishart(W,jags_data$nw, jags_data$U, log=TRUE)

  posterior_betas = array(0, dim = c(N1,jags_data$L))
  for (l in 1:jags_data$L){
    posterior_beta_prob = LaplacesDemon::dmvn(t(beta[l,,]),
                               beta_mu_hat[l,],
                               beta_cov_hat[[l]],
                               log = TRUE)

    posterior_betas[,l] = posterior_beta_prob
  }
  posterior_beta = rowSums(posterior_betas)

  posterior_W = CholWishart::dInvWishart(W,nw_post, U_post, log=TRUE)

  l2 = numeric(N1)
  for (i in 1:N1){

    likelihood = sum(LaplacesDemon::dmvn(jags_data$y, jags_data$x_pred %*% beta[,,i],W[,,i], log = TRUE))

    l2_numer = likelihood + prior_beta[i] + prior_W[i]

    l2_denumer = posterior_beta[i] + posterior_W[i]

    l2[i] = (l2_numer - l2_denumer)
  }

  beta_df <- samples_iter[paste0("beta[", 1 ,",",seq(jags_data$P), "]")]
  if (jags_data$L>1){
  for (r in 2:jags_data$L){
    beta_df =cbind(beta_df, samples_iter[paste0("beta[", r ,",",seq(jags_data$P), "]")])
  }
  }

  beta = abind::abind(apply(beta_df,1,
                     function(x) t(matrix(x,nrow=jags_data$P,ncol=jags_data$L)),
                     simplify = FALSE), along=3)

  W = abind::abind(apply(data.matrix(samples_iter[,1:(jags_data$P*jags_data$P)]),1,
                  function(x) matrix(x,nrow=jags_data$P,ncol=jags_data$P),
                  simplify = FALSE), along=3)

  prior_betas = array(0, dim = c(N2,(jags_data$L)))
  for (l in 1:jags_data$L){
    prior_beta_prob = LaplacesDemon::dmvn(t(beta[l,,]),
                           jags_data$beta_mu[l,],
                           jags_data$beta_cov[,,l], log = TRUE)
    prior_betas[,l] = prior_beta_prob
  }
  prior_beta = rowSums(prior_betas)

  prior_W = CholWishart::dInvWishart(W,jags_data$nw, jags_data$U, log = TRUE)

  posterior_betas = array(0, dim = c(N2,jags_data$L))
  for (l in 1:jags_data$L){
    posterior_beta_prob = LaplacesDemon::dmvn(t(beta[l,,]),
                               beta_mu_hat[l,],
                               beta_cov_hat[[l]],
                               log = TRUE)

    posterior_betas[,l] = posterior_beta_prob
  }
  posterior_beta = rowSums(posterior_betas)

  posterior_W = CholWishart::dInvWishart(W,nw_post, U_post, log = TRUE)

  l1 = numeric(nrow(samples_iter))
  for (i in 1:nrow(samples_iter)){

    likelihood = sum(LaplacesDemon::dmvn(jags_data$y,  jags_data$x_pred %*% beta[,,i],W[,,i], log = TRUE))

    l1_numer = likelihood  + prior_beta[i] + prior_W[i]

    l1_denumer = posterior_beta[i] +posterior_W[i]

    l1[i] = (l1_numer - l1_denumer)
  }


  tol = 1e-10 # tolerance criterion
  criterion_val = tol + 1 # criterion value

  s1 = N1/(N1 + N2)
  s2 = N2/(N1 + N2)

  lstar = median(l1)

  m_y = 0.5 # initial guess of marginal likelihood

  i = 0 # iteration counter
  while (criterion_val > tol & i<1000){
    m_y_old = m_y
    numerator = exp(l2 - lstar)/((s1*exp(l2 - lstar)) + (s2*m_y))
    denominator = 1/((s1*exp(l1 - lstar)) + (s2*m_y))
    m_y = (N1/N2) * (sum(numerator)/sum(denominator))
    logml = log(m_y) + lstar
    i = i+1
    criterion_val = abs((m_y - m_y_old)/m_y)
    if (is.na(criterion_val)){criterion_val=1}

  }
  if (i>=1000){print('The Bridge sampling does not converge')}
  return(logml)
}
