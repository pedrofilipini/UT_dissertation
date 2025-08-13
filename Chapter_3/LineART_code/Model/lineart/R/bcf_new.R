.ident <- function(...){
  # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
  args <- c(...)
  if( length( args ) > 2L ){
    #  recursively call ident()
    out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
  }else{
    out <- identical( args[1] , args[2] )
  }
  return( all( out ) )
}

cuts_classic <- function(x){
  x <- sort(unique(x))
  n <- length(x)
  x_new <- (x[-c(1)]+x[-c(n)])/2
  return(x_new)
}

.cp_quantile = function(x, num=10000, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))

  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }

  return(ret)
}

#' Fit Bayesian Causal Forests
#'
#' @references Hahn, Murray, and Carvalho(2017). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.
#'  https://arxiv.org/abs/1706.09523. (Call citation("bcf") from the
#' command line for citation information in Bibtex format.)
#'
#' @details Fits a generalized version of the Bayesian Causal Forest model (Hahn et. al. 2018): For a response
#' variable y, treatment z, and covariates x,
#' \deqn{y_i = \mu(x_i, \hat z_i) + \tau(x_i, \pi_i)\omega(z_i) + \epsilon_i}
#' where \eqn{\z_i} is an (optional) estimate of \eqn{E(Z_i | X_i=x_i)} and
#' \eqn{\epsilon_i \sim N(0,\sigma^2)}
#'
#' Some notes:
#' \itemize{
#'    \item x_control and x_moderate must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame
#'    \item sd_control and sd_moderate are the prior SD(mu(x)) and SD(tau(x)) at a given value of x (respectively). If
#'    use_muscale = FALSE, then this is the parameter \eqn{\sigma_\mu} from the original BART paper, where the leaf parameters
#'    have prior distribution \eqn{N(0, \sigma_\mu/m)}, where m is the number of trees.
#'    If use_muscale=TRUE then sd_control is the prior median of a half Cauchy prior for SD(mu(x)). If use_tauscale = TRUE,
#'    then sd_moderate is the prior median of a half Normal prior for SD(tau(x)).
#'    \item By default the prior on \eqn{\sigma^2} is calibrated as in Chipman, George and McCulloch (2008).
#'
#'
#' }
#' @param y Response variable
#' @param Omega Desgin matrix computed from treatment variable
#' @param x_control Design matrix for the "prognostic" function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param pihat Length n estimates of
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param ntree_moderate Number of trees in tau(x)
#' @param sd_moderate SD(tau(x)) marginally at any covariate value (or its prior median if use_tauscale=TRUE)
#' @param base_moderate Base for tree prior on tau(x) trees (see details)
#' @param power_moderate Power for the tree prior on tau(x) trees (see details)
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to
#' include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#' @return A list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' @examples
#'\donttest{
#'
#' #TODO: Update this example
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=2000, nsim=2000)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'
#'}
#'\dontshow{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#'}
bcf_new <- function(y, z, x_control, x_moderate=x_control, pihat,
                    z_est, x_control_est, x_moderate_est=x_control_est, pihat_est,
                    nburn, nsim, nthin = 1, update_interval = 100,
                    ntree_control = 250,
                    sd_control = 2*sd(y),
                    base_control = 0.95,
                    power_control = 2,
                    ntree_moderate = 50,
                    sd_moderate = sd(y),
                    base_moderate = 0.25,
                    power_moderate = 3,
                    nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                    include_pi = "control", linear = "moderate",
                    use_muscale=TRUE, use_tauscale=TRUE,
                    dart = FALSE, gcon = 1, gmod = 1,
                    save_trees = TRUE
) {

  pihat = as.matrix(pihat)
  pihat_est = as.matrix(pihat_est)

  if( !.ident(#length(y),
    #nrow(Omega),
    nrow(x_control_est),
    nrow(x_moderate_est),
    nrow(pihat_est)
  )
  ) {

    stop("Data size mismatch. The following should all be equal:",
         #"length(y): ", length(y), "\n",
         #"nrow(Omega): ", nrow(Omega), "\n",
         "nrow(x_control_est): ", nrow(x_control_est), "\n",
         "nrow(x_moderate_est): ", nrow(x_moderate_est), "\n",
         "nrow(pihat_est): ", nrow(pihat_est),"\n"
    )
  }

  if( !.ident(length(y),
              #nrow(Omega),
              nrow(x_control),
              nrow(x_moderate),
              nrow(pihat)
  )
  ) {

    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         #"nrow(Omega): ", nrow(Omega), "\n",
         "nrow(x_control): ", nrow(x_control), "\n",
         "nrow(x_moderate): ", nrow(x_moderate), "\n",
         "nrow(pihat): ", nrow(pihat),"\n"
    )
  }

  if(any(is.na(y))) stop("Missing values in y")
  #if(any(is.na(Omega))) stop("Missing values in Omega")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if(any(is.na(pihat))) stop("Missing values in pihat")

  if(any(is.na(x_control_est))) stop("Missing values in x_control_est")
  if(any(is.na(x_moderate_est))) stop("Missing values in x_moderate_est")
  if(any(is.na(pihat_est))) stop("Missing values in pihat_est")

  if(any(!is.finite(y))) stop("Non-numeric values in y")
  #if(any(!is.finite(Omega))) stop("Non-numeric values in Omega")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")

  if(any(!is.finite(x_control_est))) stop("Non-numeric values in x_control_est")
  if(any(!is.finite(x_moderate_est))) stop("Non-numeric values in x_moderate_est")
  if(any(!is.finite(pihat_est))) stop("Non-numeric values in pihat_est")

  #if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")

  if(length(unique(y))<5) warning("y appears to be discrete")

  if(nburn<0) stop("nburn must be positive")
  if(nsim<0) stop("nsim must be positive")
  if(nthin<0) stop("nthin must be positive")
  if(nthin>nsim+1) stop("nthin must be < nsim")
  if(nburn<100) warning("A low (<100) value for nburn was supplied")

  ### TODO range check on parameters

  ###
  x_c = matrix(x_control, ncol=ncol(x_control))
  x_m = matrix(x_moderate, ncol=ncol(x_moderate))
  if(include_pi=="both" | include_pi=="control") {
    x_c = cbind(x_control, pihat)
  }
  if(include_pi=="both" | include_pi=="moderate") {
    x_m = cbind(x_moderate, pihat)
  }

  x_c_est = matrix(x_control_est, ncol=ncol(x_control_est))
  x_m_est = matrix(x_moderate_est, ncol=ncol(x_moderate_est))
  if(include_pi=="both" | include_pi=="control") {
    x_c_est = cbind(x_control_est, pihat_est)
  }
  if(include_pi=="both" | include_pi=="moderate") {
    x_m_est = cbind(x_moderate_est, pihat_est)
  }

  #####
  n = length(y)
  n_est = length(z_est)
  if(linear=="both" | linear=="moderate") {
    Omega_mod=cbind(rep(1,n),x_m)*z
    Omega_mod_est=cbind(rep(1,n_est),x_m_est)*z_est
    linear_mod = T
  }else{
    Omega_mod=as.matrix(z)
    Omega_mod_est=as.matrix(z_est)
    linear_mod = F
  }
  if(linear=="both" | linear=="control") {
    Omega_con=cbind(rep(1,n),x_c)
    Omega_con_est=cbind(rep(1,n_est),x_c_est)
    linear_con = T
  }else{
    Omega_con=matrix(rep(1,n), ncol=1)
    Omega_con_est=matrix(rep(1,n_est), ncol=1)
    linear_con = F
  }

  #####

  cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
  cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))

  #cutpoint_list_c = lapply(1:ncol(x_c), function(i) cuts_classic(x_c[,i]))
  #cutpoint_list_m = lapply(1:ncol(x_m), function(i) cuts_classic(x_m[,i]))

  yscale = scale(y)
  sdy = sd(y)
  muy = mean(y)

  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~as.matrix(x_c))
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  dir = tempdir()

  con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy)
  mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1) # if HN make sd_moderate the prior median

  if(linear_con){
    Sigma0_con = matrix(0, nrow=2, ncol = 2)
    diag(Sigma0_con) = diag(Sigma0_con) + con_sd*con_sd/(ntree_control) #Just for the edge cases of dummy variables
  }else{
    Sigma0_con = matrix(con_sd*con_sd/(ntree_control), nrow=1)
  }

  if(linear_mod){
    Sigma0_mod = matrix(0, nrow=2, ncol = 2)
    diag(Sigma0_mod) = diag(Sigma0_mod) + mod_sd*mod_sd/(ntree_moderate) #Just for the edge cases of dummy variables
  }else{
    Sigma0_mod = matrix(mod_sd*mod_sd/(ntree_moderate), nrow=1)
  }

  print(Sigma0_con)
  print(Sigma0_mod)
  fitbcf = bcfoverparRcppClean(yscale, t(Omega_con), t(Omega_mod),
                               t(Omega_con_est), t(Omega_mod_est),
                               t(x_c), t(x_m), t(x_c_est), t(x_m_est),
                               cutpoint_list_c, cutpoint_list_m,
                               random_des = matrix(1),
                               random_var = matrix(1),
                               random_var_ix = matrix(1),
                               random_var_df = 3,
                               nburn, nsim, nthin,
                               ntree_moderate, ntree_control, lambda, nu,
                               Sigma0_con, Sigma0_mod,
                               base_control, power_control,
                               base_moderate, power_moderate,
                               tempdir(), prior_sample = FALSE,
                               use_muscale, use_tauscale,
                               -1, -1,
                               status_interval = update_interval,
                               linear_mod, linear_con, save_trees, dart, 
                               g_con = gcon, g_mod = gmod,
                               0,0)

  #B = drop(fit$post_B)
  #B0 = fit$b0
  #EYs = fit$post_yhat

  #return(List::create(_["m_post"] = m_post, _["b_post"] = b_post, _["b_est_post"] = b_est_post,
  #                     _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post,
  #                     _["gamma"] = gamma_post, _["random_var_post"] = random_var_post

  m_post = muy + sdy*fitbcf$m_post
  tau_post = sdy*fitbcf$b_post
  #yhat_post = muy + sdy*fitbcf$m_post
  #yhat_post[,z==1] = yhat_post[,z==1] + fitbcf$b_post
  #yhat_post = (muy + sdy*fitbcf$m_post)
  #yhat_post[,z[perm]==1] = yhat_post[,z[perm]==1] + sdy*fitbcf$b_post
  #yhat_post = yhat_post[,order(perm)]


  list(ymean = mean(y), 
       linear = linear,
       sigma = sdy*fitbcf$sigma,
       yhat = muy + sdy*fitbcf$yhat_post,
       mu_post = muy + sdy*fitbcf$m_post,
       mu_est_post = muy + sdy*fitbcf$m_est_post,
       b_post = sdy*fitbcf$b_post,
       b_est_post = sdy*fitbcf$b_est_post,
       #       mu  = m_post,
       varprobs_con = fitbcf$probs_con_post,
       varprobs_mod = fitbcf$probs_mod_post,
       #varalpha_con = fitbcf$alpha_con_post,
       #varalpha_mod = fitbcf$alpha_mod_post,
       tau = sdy*fitbcf$b_post,
       gamma_con = fitbcf$gamma_con_post,
       eta_con = sdy*fitbcf$eta_con_post,
       gamma_mod = fitbcf$gamma_mod_post,
       eta_mod = sdy*fitbcf$eta_mod_post,
       mu_scale = fitbcf$msd*sdy,
       tau_scale = fitbcf$bsd*sdy,
       coefs_con = fitbcf$coefs_con*sdy,
       coefs_con_est = fitbcf$coefs_con_est*sdy,
       coefs_mod = fitbcf$coefs_mod*sdy,
       coefs_mod_est = fitbcf$coefs_mod_est*sdy,
       treedraws_con = fitbcf$treedraws_con,
       treedraws_mod = fitbcf$treedraws_mod
  )

}
