#include "slice.h"

double slice(double x0, //current x value
             logdensity* g, //log density of distribution
             double w, //Defines step size for interval creation around x0 for new x
             double m, //Useless, not used anywhere
             double lower, //Lower bound of x that can be sampled
             double upper) { //Upper bound
  double x1; //new x value

  //Rcpp::Rcout << "Before evaluation" << endl;
  double gx0 = g->val(x0); //Evaluate log(f(x0)) - Using log avoids underflow
  //Rcpp::Rcout << "l: " << x0 << endl;
  //Rcpp::Rcout << "gx0: " << gx0 << endl;
  double logy = gx0 - R::rexp(1.); //because if y ~ u (a,b), c ~ exp (k), then y = a +(b-a) exp (-k c), so just apply log
  double u = R::runif(0., w); //Sampling the step size
  double L = x0 - u; //now for a lower bound
  double R = x0 + (w-u); //and an upper bound

  while(true) {
    R_CheckUserInterrupt();
    if(L<=lower) { break; } //If L gets below the values that x can assume, stop
    if(g->val(L)<=logy) { break; } // If the density of the point L gets lower than f(x0), then I do not sample from there anymore, so stop
    L -= w; //If the above requirements were not fulfilled, push L even further to the left
  }
  //Rcpp::Rcout << "Left " << endl;
  while(true) {
    R_CheckUserInterrupt();
    if(R>=upper) { break; } //If R gets above the values that x can assume, stop
    if(g->val(R)<=logy) { break; } // If the density of the point R gets lower than f(x0), then I do not sample from there anymore, so stop
    R += w; //If the above requirements were not fulfilled, push R even further to the right
  }
  //Rcpp::Rcout << "Right " << endl;

  if(L<lower) {L=lower;} //If L is below the minimum value allowed for x, use the lower bound
  if(R>upper) {R=upper;} //If R is over the maximum value allowed for x, use the upper bound

  while(true) {
    R_CheckUserInterrupt();
    x1 = R::runif(L, R); // Sample new x value from uniform (L,R)
    double gx1 = g->val(x1); //Calculate the log(f(x new))
    if(gx1>=logy) { break; } //If f(x new) is greater than y, stop and accept x new
    if(x1>x0) { //If f(x new) is less than y, then update the sampling interval
      R = x1; //Update R if x new is to the right of x old
    } else {
      L = x1; //Update L if x new is to the left of x old
    }
  }
  //Rcpp::Rcout << "lnew: " << x1 << endl;


  return(x1);
}
