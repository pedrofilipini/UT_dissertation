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
bcf_core <- function(y, tvar_con, tvar_mod = tvar_con, z_train = NULL, z_out = z_train,
                            Omega_con=matrix(rep(1,length(y)), ncol=1), Omega_mod=Omega_con,
                            x_control, x_moderate = x_control,
                            tvar_con_out, tvar_mod_out = tvar_con_out,
                            Omega_con_out = Omega_con, Omega_mod_out = Omega_mod,
                            control_fits = TRUE, moderate_fits = TRUE, # Only returns coefficients for out of sample x's
                            x_control_out = x_control, x_moderate_out = x_moderate,
                            pihat = rep(0.5, length(y)), pihat_out = rep(0.5, nrow(x_control_out)),
                            randeff_design = matrix(1),
                            randeff_variance_component_design = matrix(1),
                            randeff_scales = 1,
                            randeff_df = 3,
                            nburn, nsim, nthin = 1, update_interval = 100,
                            ntree_control = 250,
                            sd_control = 2*sd(y),
                            base_control = 0.95,
                            power_control = 2,
                            ntree_moderate = 50,
                            sd_moderate = sd(y),
                            base_moderate = 0.25,
                            power_moderate = 3,
                            j_con = 11, j_mod = 8,
                            enc_mod_prior = 1, #This should be lower than the lower bound
                            enc_con_prior = 1,
                            enc_con = 3,
                            enc_mod = 3, #Expected number of crossings
                            dart = FALSE, var_sizes_con = 1, var_sizes_mod = 1,
                            update_ls = TRUE,
                            nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                            include_pi = "control",
                            use_muscale=TRUE, use_tauscale=TRUE,
                            vanilla = all((Omega_con-1)<1e-8),
                            debug = FALSE
) {

  # TODO: checks for random effects matrices

  pihat = as.matrix(pihat)
  if( !.ident(length(y),
              #nrow(Omega_con),
              #nrow(Omega_mod),
              nrow(x_control),
              nrow(x_moderate),
              nrow(pihat)
  )
  ) {

    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         #"nrow(Omega_con): ", nrow(Omega_con), "\n",
         #"nrow(Omega_mod): ", nrow(Omega_mod), "\n",
         "nrow(x_control): ", nrow(x_control), "\n",
         "nrow(x_moderate): ", nrow(x_moderate), "\n",
         "nrow(pihat): ", nrow(pihat),"\n"
    )
  }

  if( !.ident(#nrow(Omega_con_out),
              #nrow(Omega_mod_out),
              nrow(x_control_out),
              nrow(x_moderate_out),
              nrow(pihat_out)
  )
  ) {

    stop("Data size mismatch. The following should all be equal:",
         #"nrow(Omega_con_out): ", nrow(Omega_con_out), "\n",
         #"nrow(Omega_mod_out): ", nrow(Omega_mod_out), "\n",
         "nrow(x_control_out): ", nrow(x_control_out), "\n",
         "nrow(x_moderate_out): ", nrow(x_moderate_out), "\n",
         "nrow(pihat_out): ", nrow(pihat_out),"\n"
    )
  }

  if(any(is.na(y))) stop("Missing values in y")
  #if(any(is.na(Omega_con))) stop("Missing values in Omega_con")
  #if(any(is.na(Omega_mod))) stop("Missing values in Omega_mod")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if(any(is.na(pihat))) stop("Missing values in pihat")

  if(any(!is.finite(y))) stop("Non-numeric values in y")
  #if(any(!is.finite(Omega_con))) stop("Non-numeric values in Omega_con")
  #if(any(!is.finite(Omega_mod))) stop("Non-numeric values in Omega_mod")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")

  #if(any(is.na(Omega_con_out))) stop("Missing values in Omega_con_out")
  #if(any(is.na(Omega_mod_out))) stop("Missing values in Omega_mod_out")
  if(any(is.na(x_control_out))) stop("Missing values in x_control_out")
  if(any(is.na(x_moderate_out))) stop("Missing values in x_moderate_out")
  if(any(is.na(pihat_out))) stop("Missing values in pihat_out")

  #if(any(!is.finite(Omega_con_out))) stop("Non-numeric values in Omega_con_out")
  #if(any(!is.finite(Omega_mod_out))) stop("Non-numeric values in Omega_mod_out")
  if(any(!is.finite(x_control_out))) stop("Non-numeric values in x_control_out")
  if(any(!is.finite(x_moderate_out))) stop("Non-numeric values in x_moderate_out")
  if(any(!is.finite(pihat_out))) stop("Non-numeric values in pihat_out")

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
  x_c_o = matrix(x_control_out, ncol=ncol(x_control_out))
  x_m_o = matrix(x_moderate_out, ncol=ncol(x_moderate_out))

  if(include_pi!="neither"){
    if(length(unique(pihat)) == 1){
      warning("All values of pihat are equal. Pihat will not be included among covariates")
      include_pi="neither"
    }
  }

  if(include_pi=="both" | include_pi=="control") {
    x_c = cbind(x_control, pihat)
    x_c_o = cbind(x_control_out, pihat_out)
  }
  if(include_pi=="both" | include_pi=="moderate") {
    x_m = cbind(x_moderate, pihat)
    x_m_o = cbind(x_moderate_out, pihat_out)
  }


  cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
  cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))

  yscale = scale(y)
  sdy = sd(y)
  muy = mean(y)
  n = length(y)

  tvar_con_out = (tvar_con_out-(max(tvar_con_out)-min(tvar_con_out))/2-min(tvar_con_out))/((max(tvar_con_out)-min(tvar_con_out))/2)
  tvar_con = (tvar_con-(max(tvar_con)-min(tvar_con))/2-min(tvar_con))/((max(tvar_con)-min(tvar_con))/2)
  tvar_mod_out = (tvar_mod_out-(max(tvar_mod_out)-min(tvar_mod_out))/2-min(tvar_mod_out))/((max(tvar_mod_out)-min(tvar_mod_out))/2)
  tvar_mod = (tvar_mod-(max(tvar_mod)-min(tvar_mod))/2-min(tvar_mod))/((max(tvar_mod)-min(tvar_mod))/2)

  ###################

  #Starting here we gather everything for the mod
  softmin_mod <- (max(tvar_mod)-min(tvar_mod))/(pi*enc_mod) #put an if for softmin
  l_mod <- softmin_mod
  mc_mod <- m_select(a=1, mmin=5, mmax=100,cmin=1, cmax=10, step = 0.01, softmin=l_mod,xmin=min(tvar_mod),xmax=max(tvar_mod),e=0.01)
  print(paste0("Selected j_mod: ", mc_mod[1]))
  lmax_mod <- lmax_select(a=1, m=mc_mod[1], c=5, step = 0.01, lmin=0.01, xmin=min(tvar_mod),xmax=max(tvar_mod))

  l_mod_aux <- seq(0.01, lmax_mod[1], length.out = 100)
  c_mod_aux <- numeric(0)
  for(ll_mod in 1:length(l_mod_aux)){
    c_mod_aux[ll_mod] <- c_select(a=1, m=mc_mod[1], cmin=1, cmax=10, step = 0.01, l=l_mod_aux[ll_mod], xmin=min(tvar_mod), xmax=max(tvar_mod))[1]
    #print(paste0("C_mod: ", c_mod_aux[ll_mod], "; l: ", l_mod_aux[ll_mod]))
  }

  #plot(l_mod_aux, c_mod_aux)
  L_mod <- mc_mod[2]*(max(tvar_mod)-min(tvar_mod))/2
  Phi_mod = t(Omega(tvar_mod-(max(tvar_mod)+min(tvar_mod))/2, mc_mod[1], L_mod))
  Phi_mod_out = t(Omega(tvar_mod_out-(max(tvar_mod_out)+min(tvar_mod_out))/2, mc_mod[1], L_mod))

  #Now the same for the con
  if(!vanilla){
    softmin_con <- (max(tvar_con)-min(tvar_con))/(pi*enc_con)
    l_con <- softmin_con
    mc_con <- m_select(a=1, mmin=5, mmax=100,cmin=1, cmax=10, step = 0.01, softmin=l_con,xmin=min(tvar_con),xmax=max(tvar_con),e=0.01)
    print(paste0("Selected j_con: ", mc_con[1]))
    lmax_con <- lmax_select(a=1, m=mc_con[1], c=5, step = 0.01, lmin=0.01, xmin=min(tvar_con),xmax=max(tvar_con))

    l_con_aux <- seq(0.01, lmax_con[1], length.out = 100)
    c_con_aux <- numeric(0)
    for(ll_con in 1:length(l_con_aux)){
      c_con_aux[ll_con] <- c_select(a=1, m=mc_con[1], cmin=1, cmax=10, step = 0.01, l=l_con_aux[ll_con], xmin=min(tvar_con), xmax=max(tvar_con))[1]
      #print(paste0("C_con: ", c_con_aux[ll_con], "; l: ", l_con_aux[ll_con]))
    }

    #plot(l_con_aux, c_con_aux)
    L_con <- mc_con[2]*(max(tvar_con)-min(tvar_con))/2
    Phi_con = (t(Omega(tvar_con-(max(tvar_con)+min(tvar_con))/2, mc_con[1], L_con)))
    Phi_con_out = (t(Omega(tvar_con_out-(max(tvar_con_out)+min(tvar_con_out))/2, mc_con[1], L_con)))

    ab <- c(1,1) #placeholder

    Omega_con <- t(t(Omega(tvar_con-(max(tvar_con)+min(tvar_con))/2, mc_con[1], L_con) %*% sqrt(SDiag(1, l_con, mc_con[1], L_con))))
    Omega_con_out <- t(t(Omega(tvar_con_out-(max(tvar_con)+min(tvar_con))/2, mc_con[1], L_con) %*% sqrt(SDiag(1, l_con, mc_con[1], L_con))))
  }

  if(is.null(z_train)){
    Omega_mod <- Omega(tvar_mod-(max(tvar_mod)+min(tvar_mod))/2, mc_mod[1], L_mod) %*% sqrt(SDiag(1, l_mod, mc_mod[1], L_mod))
    Omega_mod_out <- Omega(tvar_mod_out-(max(tvar_mod_out)+min(tvar_mod_out))/2, mc_mod[1], L_mod) %*% sqrt(SDiag(1, l_mod, mc_mod[1], L_mod))
  }else{
    Omega_mod <- (Omega(tvar_mod-(max(tvar_mod)+min(tvar_mod))/2, mc_mod[1], L_mod) %*% sqrt(SDiag(1, l_mod, mc_mod[1], L_mod)))*z_train
    Omega_mod_out <- (Omega(tvar_mod_out-(max(tvar_mod_out)+min(tvar_mod_out))/2, mc_mod[1], L_mod) %*% sqrt(SDiag(1, l_mod, mc_mod[1], L_mod)))*z_out
  }
  ####################
  #length-scale priors
  ls_prior_mod <- numeric(0)
  ls_prior_mod[1] <- 5
  ls_prior_mod[2] <- ((max(tvar_mod)-min(tvar_mod))/(pi*enc_mod_prior))/(sqrt(2)*pracma::erfinv(0.9))

  ls_prior_con <- numeric(0)
  ls_prior_con[1] <- 5
  ls_prior_con[2] <- ((max(tvar_con)-min(tvar_con))/(pi*enc_con_prior))/(sqrt(2)*pracma::erfinv(0.9))

  ####################

  con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy)
  mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1) # if HN make sd_moderate the prior median

  if(vanilla){
    Omega_con = matrix(rep(1,n), ncol=1)
    Omega_con_out = matrix(rep(1,nrow(x_control_out)), ncol=1)
    Sigma0_con = matrix(con_sd*con_sd/(ntree_control), nrow=1)
    mc_con <- 1
    l_con <- 1
    l_con_aux <- 1
    c_con_aux <- 1
    Phi_con <- 1
    Phi_con_out <- 1
  }else{
    Sigma0_con = matrix(0, nrow=ncol(Omega_con), ncol = ncol(Omega_con))
    diag(Sigma0_con) = diag(Sigma0_con) + 1#con_sd*con_sd/(ntree_control*ncol(Omega_con))
  }


  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~Omega_con+as.matrix(x_c))
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  dir = tempdir()


  if(length(unique(Omega_mod))==1) stop("Treatment has one unique value")
  Sigma0_mod = matrix(0, nrow=ncol(Omega_mod), ncol = ncol(Omega_mod))
  diag(Sigma0_mod) = diag(Sigma0_mod) + 1#mod_sd*mod_sd/(ntree_moderate*ncol(Omega_mod))

  print("Finished here. /n")

  const_args =
    list(y_ = yscale, z_ = z_train, z_est_ = z_out,
                      tvar_con_ = tvar_con, tvar_mod_ = tvar_mod,
                      tvar_con_est_ = tvar_con_out, tvar_mod_est_ = tvar_mod_out,
                      Omega_con = t(Omega_con), Omega_mod = t(Omega_mod),
                      Omega_con_est = t(Omega_con_out), Omega_mod_est = t(Omega_mod_out),
                      Phi_con = Phi_con, Phi_mod = Phi_mod,
                      Phi_con_est = Phi_con_out, Phi_mod_est = Phi_mod_out,
                      l_con = l_con, l_mod = l_mod, L = L_mod, j_con_ = mc_con[1], j_mod_ = mc_mod[1],
                      lvec_con_ = l_con_aux, cvec_con_ = c_con_aux, lvec_mod_ = l_mod_aux, cvec_mod_ = c_mod_aux,
                      t(x_c), t(x_m), t(x_c_o), t(x_m_o),
                      cutpoint_list_c, cutpoint_list_m,
                      random_des = randeff_design,
                      random_var = matrix(rep(0.00000001, ncol(randeff_variance_component_design)), ncol=1),
                      random_var_ix = randeff_variance_component_design,
                      random_var_df = randeff_df, randeff_scales = randeff_scales/sdy,
                      lprior_con_ = ls_prior_con, lprior_mod_ = ls_prior_mod,
                      burn = nburn, nd = nsim, thin = nthin,
                      ntree_mod = ntree_moderate, ntree_con = ntree_control,
                      lambda = lambda, nu = nu,
                      Sigma0_con = Sigma0_con, Sigma0_mod =Sigma0_mod,
                      #Sigma0_con=2*matrix(1/250), Sigma0_mod=1*matrix(0.02),
                      con_alpha = base_control, con_beta = power_control,
                      mod_alpha = base_moderate, mod_beta = power_moderate,
                      treef_name = tempdir(),
                      #est_mod_fits = moderate_fits, est_con_fits = control_fits,
                      prior_sample = FALSE,
                      use_con_scale = use_muscale, use_mod_scale = use_tauscale,
                      con_scale_df = 1, mod_scale_df = -1,
                      status_interval = update_interval,
                      vanilla = vanilla,
                      dart = dart, update_ls= update_ls, var_sizes_con = var_sizes_con, var_sizes_mod = var_sizes_mod
    )

  fitbcf = do.call(bcfCore, const_args)

  out = list(sigma = sdy*fitbcf$sigma,
       yhat = muy + sdy*fitbcf$yhat_post,
       mu_post = muy + sdy*fitbcf$m_post,
       mu_est_post = muy + sdy*fitbcf$m_est_post,
       tau_omega_post = sdy*fitbcf$b_post,
       tau_omega_est_post = sdy*fitbcf$b_est_post,
       #       mu  = m_post,
       #tau = tau_post,
       mu_scale = fitbcf$msd*sdy,
       tau_scale = fitbcf$bsd*sdy,
       coefs_mod = fitbcf$coefs_mod*sdy,
       coefs_con = fitbcf$coefs_con*sdy,
       l_con_post = fitbcf$l_con_post,
       l_mod_post = fitbcf$l_mod_post,
       L_con_post = fitbcf$L_con_post,
       L_mod_post = fitbcf$L_mod_post,
       eta_con = fitbcf$eta_con_post,
       eta_mod = fitbcf$eta_mod_post,
       #_["coefs_mod_est"] = scoefs_mod_est,
       #_["coefs_con_est"] = scoefs_con_est,
       coefs_mod_est = fitbcf$coefs_mod_est*sdy,
       coefs_con_est = fitbcf$coefs_con_est*sdy,
       random_effects = fitbcf$gamma*sdy,
       random_effects_sd = fitbcf$random_sd_post*sdy,
       splitting_prob_con = fitbcf$var_prob_con,
       splitting_prob_mod = fitbcf$var_prob_mod
  )

  if(debug) {
    out$const_args = const_args
    out$raw_fit = fitbcf
  }


  return(out)
}
