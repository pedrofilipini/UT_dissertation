#include "arma_config.h"
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

#include "cselect.h"

using namespace Rcpp;

// [[Rcpp::export]]
List cSelect(NumericVector x_, NumericVector y_, double c_max, double value)
{
  /*****************************************************************************
  * Read, format y
  *****************************************************************************/
  std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY;

  for(NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) {
    y.push_back(*it);
  }
  size_t n = y.size();

  std::vector<double> x; //storage for y
  double minx = INFINITY, maxx = -INFINITY;

  for(NumericVector::iterator it=x_.begin(); it!=x_.end(); ++it) {
    x.push_back(*it);
  }

  double result;
  c_select_interp csel(x, y, n);
  result = csel.val(value);

  return(List::create(_["x"] = x,
                       _["y"] = y,
                       _["x_length"] = n,
                       _["result"] = result
  ));
}
