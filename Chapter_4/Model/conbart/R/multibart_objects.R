make_bart_designs = function(X_list, basis_matrix_list) {
  num_designs = length(X_list)
  if(num_designs != length(basis_matrix_list)) {
    stop("X_list and basis_matrix_list must be the same size")
  }
  
  lapply(1:num_designs, function(i) make_bart_design(X_list[[i]], basis_matrix_list[[i]], i-1))
  
}

make_bart_design = function(
  X,
  basis_matrix,
  index = -1
  ) {
  
  cutpoint_list = lapply(1:ncol(X), function(i) multibart:::.cp_quantile(X[,i]))
  
  list(X=t(X),
       Omega = t(basis_matrix),
       info = cutpoint_list,
       index = index
  )
}

make_bart_spec = function(
  design,
  ntree,
  Sigma0,
  mu0 = NA,
  scale_df = -1,
  vanilla = FALSE,
  alpha = 0.95, 
  beta = 2,
  dart = FALSE
  ) {
  
  if(is.na(mu0)) mu0 = rep(0, ncol(Sigma0))
  if(ncol(design$Omega)==1 & all( (design$Omega - 1)<1e-8)) vanilla = TRUE
  
  list(
    design_index = design$index,
    ntree = ntree,
    scale_df = scale_df,
    mu0 = mu0,
    Sigma0 = Sigma0,
    vanilla = vanilla,
    alpha = alpha,
    beta = beta,
    dart = dart
  )
}