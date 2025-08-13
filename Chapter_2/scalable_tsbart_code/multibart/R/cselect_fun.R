#' Testing the c_select function
#'
c_select_fun <- function(x, y, cmax, value) {

  const_args = list(x_ = x, y_ = y, c_max= cmax, value = value)

  cselect_res = do.call(cSelect, const_args)

  return(cselect_res)
}
