##Version for working with basis ts-bart

##################################
# Useful Functions#
##################################
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

##################################
##################################


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
#' @param x_control Design matrix for the function mu(x)
#' @param pihat Length n estimates of
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_pi Takes values "control" or "none".
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @return A list with elements
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' @examples

lbart_new <- function(y, #response
                    Omega, #basis functions
                    x_control, #covariates
                    Omega_est, #basis functions
                    x_control_est, #covariates
                    #pihat, #propensity score
                    nburn = 1000, #burn-in
                    nsim = 1000, #posterior samples
                    nthin = 1, #thinning
                    update_interval = 100, #print every
                    ntree_control = 200, #number of trees
                    sd_control = 2*sd(y), #Initial guess for variance
                    base_control = 0.95, #alpha hyperparameter
                    power_control = 2, #beta hyperparameter
                    nu = 3, lambda = NULL, sigq = .9, sighat = NULL, #for the sigma prior
                    include_pi = "control", #include propensity score?
                    use_muscale=FALSE #use half-Cauchy?
) {

  #pihat = as.matrix(pihat) #set matrix of propensity score to facilitate the calculations
  if( !.ident(length(y), #Testing if those are equal (more information on the beginning of the function)
              nrow(Omega),
              nrow(x_control) #,
              # nrow(pihat)
  )
  ) {
    #If are equal, then stop the whole thing
    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         "nrow(Omega): ", nrow(Omega), "\n",
         "nrow(x_control): ", nrow(x_control), "\n" #,
         #"nrow(pihat): ", nrow(pihat),"\n"
    )
  }

  #Look for NAs
  if(any(is.na(y))) stop("Missing values in y")
  if(any(is.na(Omega))) stop("Missing values in Omega")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  #if(any(is.na(pihat))) stop("Missing values in pihat")

  #Look for infinity values on the input data
  if(any(!is.finite(y))) stop("Non-numeric values in y")
  if(any(!is.finite(Omega))) stop("Non-numeric values in Omega")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  #if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")

  #Verify if Y is discrete (only a warning, might be a signal of an error for the user)
  if(length(unique(y))<5) warning("y appears to be discrete")


  #Verify if the number of samples make sense
  if(nburn<0) stop("nburn must be positive")
  if(nsim<0) stop("nsim must be positive")
  if(nthin<0) stop("nthin must be positive")
  if(nthin>nsim+1) stop("nthin must be < nsim")
  if(nburn<100) warning("A low (<100) value for nburn was supplied")

  #####################################
  ### TODO range check on parameters###
  #####################################

  #Define covariate matrix and include pi, if needed
  x_c = matrix(x_control, ncol=ncol(x_control))
  x_c_est = matrix(x_control_est, ncol=ncol(x_control_est))
  #if(include_pi=="control") {
  #  x_c = cbind(x_control, pihat)
  #}


  #Defining the cutpoints
  cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))

  #Defining some new variable with statistics that are interesting
  yscale = scale(y)
  sdy = sd(y)
  muy = mean(y)

  #Calculate the specs for the sigma prior
  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~Omega+as.matrix(x_c))
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  #define temporary directory for (maybe) saving stuff
  dir = tempdir()

  #scale of mu
  con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy)

  #Variance over the number of trees
  Sigma0_con = matrix(0, nrow=2, ncol = 2)
  diag(Sigma0_con) = diag(Sigma0_con) + con_sd*con_sd/(ntree_control*2)


  #Calling the Rcpp function
  fitlbart = lbartRcpp(yscale,
                       t(Omega),
                       t(Omega_est), #For prediction
                       t(x_c),
                       t(x_c_est), #For prediction
                       cutpoint_list_c,
                       nburn, nsim, nthin,
                       ntree_control, lambda, nu,
                       Sigma0_con,
                       base_control, power_control,
                       tempdir(), prior_sample = FALSE,
                       use_muscale,
                       1,
                       status_interval = update_interval)



  ##################################################
  ##########TODO: WORK ON THE OUTPUT!!##############
  ##################################################
  mu_est_post = muy + sdy*fitlbart$m_est_post

  m_post = muy + sdy*fitlbart$m_post

  list(sigma = sdy*fitlbart$sigma,
       yhat = muy + sdy*fitlbart$yhat_post,
       mu_post = muy + sdy*fitlbart$m_post,
       mu_est_post = muy + sdy*fitlbart$m_est_post,
       #mu_scale = fitlbart$msd*sdy,
       coefs = fitlbart$coefs_con_est*sdy,
       treedraws = fitlbart$treedraws
  )

}
