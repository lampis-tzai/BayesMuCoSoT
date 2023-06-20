library(R2jags)

#' Bayesian Multivariate Normal Modelling for Common Source Testing
#'
#' @param y Vector of response columns names (p=variables).
#' @param x Predictors columns names (d=predictors). If there are any.
#' @param questioned_data Data frame with the unknown data.
#' @param known_data Data frame of the known data.
#' @param background_data Background data frame to estimate the prior parameters.
#' Background data recommended.
#' @param background_data_id Column name of the background data frame for
#' discriminate the data per id.
#' @param approach "conjugate" or "independent prior" approach.
#' Default "independent prior"
#' @param BayesFactor_Approximation "Laplace-Metropolis" or
#' "Generalized Harmonic Mean" or "Bridge Sampling". Default "Bridge Sampling".
#'  For the conjugate approach the Bayes factor comes in close form.
#' @param prior_elicitation "Non-Informative" or "Maximum Likelihood"
#' or list of type list(U= the inverse Wishart's scale matrix p by p,
#' beta_mu= matrix mean of the matrix Normal d by p,
#' beta_cov= d covariances matrix (3D matrix p by p by d).
#' Default "Maximum Likelihood" estimations from the Background data.
#' @param DoF Degrees of freedom for the prior inverse Wishart distribution.
#' Default "min": p+2. Otherwise, integer bigger than p+2.
#' @param verbose messages while running. Default TRUE.
#' @param ... the other arguments of jags {R2jags} function.
#' @return the logarithmic Bayes factor
#' @export
BayesMuCoSoT_fit <- function(y, x = NA, questioned_data, known_data,
                             background_data = NA, background_data_id = NA,
                             approach = "independent prior",
                             BayesFactor_Approximation = "Bridge Sampling",
                             prior_elicitation = "Maximum Likelihood",
                             DoF = "min", verbose = TRUE, ...){


  if (!length(background_data)>1){
    prior_elicitation = "Non-Informative"
  } else if (length(background_data)>1 & is.na(background_data_id)){
    print("To use background data for the prior elicitation we need a column id (at least 2 ids)")
    stop()
  }


  if (!length(x)>1){
    x = 'x'
    questioned_data['x'] = known_data['x'] = 1
    if (length(background_data)>1){
      background_data['x'] = 1
    }
  }


  all_together = rbind(questioned_data, known_data)
  p = length(y)
  l = length(x)

  if (DoF=="min"){
    nw_hat = p + 2
  } else if(is.integer(DoF) & DoF>=p + 2){
      nw_hat = DoF
  }else {print("DoF must be an integer and bigger than p+1")
        stop()}


  if (prior_elicitation == "Maximum Likelihood"){
    background_y = as.matrix(background_data[,y])
    background_x = as.matrix(background_data[,x])

    beta = solve(t(background_x) %*% background_x) %*% t(background_x) %*% background_y

    cov_beta = array(0, dim=c(p,p,l))
    S_w = 0
    for (w in unique(background_data[,background_data_id])){
      df_id = background_data[(background_data[,background_data_id]==w),]

      bid_x = as.matrix(df_id[,x])

      bid_y =  as.matrix(df_id[,y])

      beta_w = solve(t(bid_x) %*% bid_x) %*% t(bid_x) %*% bid_y

      for (l_id in 1:l){
        S.this <- (beta_w[l_id,] - beta[l_id,]) %*% t(beta_w[l_id,] - beta[l_id,])
        cov_beta[,,l_id] <- cov_beta[,,l_id] + S.this
      }

      res = bid_y - bid_x %*% beta_w
      Cov.this = (t(res) %*% res)/(nrow(df_id)-l)
      S_w = S_w + Cov.this
    }

    beta_cov_all = cov_beta/(length(unique(background_data[,background_data_id])) - 1)
    W_hat <- S_w/length(unique(background_data[,background_data_id]))
    U_hat <- W_hat * (nw_hat - p - 1)

  }else if (prior_elicitation == "Non-Informative"){
    U_hat = diag(5,p)
    beta = array(0,c(l,p))
    beta_cov_all = array(rep(diag(5,p),l), dim =c(p,p,l))
  } else{
    U_hat = prior_elicitation$U
    beta = prior_elicitation$beta_mu
    beta_cov_all = prior_elicitation$beta_cov
    }




  jags_data_H0 <- list(
    "P" = p,
    "L" = l,
    "N" = nrow(all_together),
    "x_pred" =  unname(as.matrix(all_together[,x])),
    "y" = unname(as.matrix(all_together[,y])),
    "U"= as.matrix(unname(U_hat)),
    "nw" = nw_hat,
    "beta_mu" = beta,
    "beta_cov" = beta_cov_all
  )


  jags_data_H1_1 <- list(
    "P" = p,
    "L" = l,
    "N" = nrow(questioned_data),
    "x_pred" =  unname(as.matrix(questioned_data[,x])),
    "y" = unname(as.matrix(questioned_data[,y])),
    "U"= as.matrix(unname(U_hat)),
    "nw" = nw_hat,
    "beta_mu"=beta,
    "beta_cov"=beta_cov_all
  )


  jags_data_H1_2 <- list(
    "P" = p,
    "L" = l,
    "N" = nrow(known_data),
    "x_pred" =  unname(as.matrix(known_data[,x])),
    "y" = unname(as.matrix(known_data[,y])),
    "U"= as.matrix(unname(U_hat)),
    "nw" = nw_hat,
    "beta_mu"=beta,
    "beta_cov"=beta_cov_all
  )

  if (approach == "independent prior"){
    params <- c("beta", "W")

    if (l>1){
    model.string <-"
     model
    {

    for ( n in 1:N ) {

    y[n,1:P] ~ dmnorm(x_pred[n,1:L]%*%beta[1:L,1:P],W.inv[1:P,1:P]);
    }

    ### Define the priors
    W.inv[1:P,1:P] ~ dwish(U[1:P,1:P], nw);
    W[1:P,1:P] = inverse(W.inv[1:P,1:P]);
    for (l in 1:L){
    beta[l,1:P]~ dmnorm(beta_mu[l,1:P], inverse(beta_cov[1:P,1:P,l]));
    }
    }
    "
    }else {
      model.string <-"
     model
    {

    for ( n in 1:N ) {

    y[n,1:P] ~ dmnorm(x_pred[n,1:L]*beta[1:L,1:P],W.inv[1:P,1:P]);
    }

    ### Define the priors
    W.inv[1:P,1:P] ~ dwish(U[1:P,1:P], nw);
    W[1:P,1:P] = inverse(W.inv[1:P,1:P]);
    for (l in 1:L){
    beta[l,1:P]~ dmnorm(beta_mu[l,1:P], inverse(beta_cov[1:P,1:P,l]));
    }
    }
    "
    }

    if (verbose){
      print("Prior elicitation finished")
      bar = "text"
    }else{
      bar = "none"
      }

    if (verbose){print("Sampling from posterior considering common modeling per source")}

    samps_H0 <- jags(jags_data_H0, parameters.to.save = params,
                     model.file =  textConnection(model.string),
                     progress.bar = bar,quiet = TRUE, ...)

    if (verbose){print("Sampling from posterior considering different modeling per source")}

    samps_H1_1 <- jags(jags_data_H1_1, parameters.to.save = params,
                       model.file =  textConnection(model.string),
                       progress.bar = bar,quiet = TRUE, ...)


    samps_H1_2 <- jags(jags_data_H1_2, parameters.to.save = params,
                       model.file =  textConnection(model.string),
                       progress.bar = bar,quiet = TRUE, ...)

    if (verbose){print("Calculating Bayes factor")}

    if (BayesFactor_Approximation == "Laplace-Metropolis"){
      laplace_lik_H0 = marginal_likelihood_laplace_metr(samps_H0$BUGSoutput$sims.matrix,
                                                        jags_data_H0)

      laplace_lik_H1_1 = marginal_likelihood_laplace_metr(samps_H1_1$BUGSoutput$sims.matrix,
                                                          jags_data_H1_1)

      laplace_lik_H1_2 = marginal_likelihood_laplace_metr(samps_H1_2$BUGSoutput$sims.matrix,
                                                          jags_data_H1_2)

      laplace_logbf = laplace_lik_H0-laplace_lik_H1_1-laplace_lik_H1_2

      return(laplace_logbf)

    } else if (BayesFactor_Approximation == "Generalized Harmonic Mean"){
      ghm_lik_H0 = ml_generalized_harmonic_mean(samps_H0$BUGSoutput$sims.matrix,
                                                jags_data_H0)
      ghm_lik_H1_1 = ml_generalized_harmonic_mean(samps_H1_1$BUGSoutput$sims.matrix,
                                                  jags_data_H1_1)
      ghm_lik_H1_2 = ml_generalized_harmonic_mean(samps_H1_2$BUGSoutput$sims.matrix,
                                                  jags_data_H1_2)

      ghm_logbf = ghm_lik_H0-ghm_lik_H1_1-ghm_lik_H1_2

      return(ghm_logbf)
    } else if (BayesFactor_Approximation == "Bridge Sampling"){

      bs_lik_H0 = ml_bridge_sampling(samps_H0$BUGSoutput$sims.matrix,
                                     jags_data_H0)
      bs_lik_H1_1 = ml_bridge_sampling(samps_H1_1$BUGSoutput$sims.matrix,
                                       jags_data_H1_1)
      bs_lik_H1_2 = ml_bridge_sampling(samps_H1_2$BUGSoutput$sims.matrix,
                                       jags_data_H1_2)

      bs_logbf = bs_lik_H0-bs_lik_H1_1-bs_lik_H1_2
      return(bs_logbf)
    } else {
      return("Not right Bayes factor approximation")
      }


  }else if (approach == "conjugate"){
    conjugate_lik_H0 = marginal_likelihood_conjugate(jags_data_H0)
    conjugate_lik_H1_1 = marginal_likelihood_conjugate(jags_data_H1_1)
    conjugate_lik_H1_2 = marginal_likelihood_conjugate(jags_data_H1_2)

    conjugate_logbf = conjugate_lik_H0-conjugate_lik_H1_1-conjugate_lik_H1_2
    return(conjugate_logbf)
  } else{
    return("Not right approach")
    }

}

