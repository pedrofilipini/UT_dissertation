#ifndef slice_h
#define slice_h

#include "arma_config.h"
#include <RcppArmadillo.h>

#include <vector>
#include "funs.h"
#include "info.h"
#include <tgmath.h>
#include "cselect.h"

class logdensity {
  public:
  virtual double val(double y) = 0;
};

class ld_norm: public logdensity {
  public:
  double mu;
  double sigma;
  double val(double y) { return(R::dnorm(y, mu, sigma, 1)); }
  ld_norm(double mu_, double sigma_) { mu=mu_; sigma=sigma_; }
};

class ld_basis_lengthscale: public logdensity {
public:
  double* data;
  arma::mat coefs;
  pinfo pi;
  dinfo di;
  std::vector<double>* lprior;
  c_select_interp csel;
  bool bcf;

  ld_basis_lengthscale(double* data_, arma::mat &coefs_, pinfo& pi_, dinfo& di_, std::vector<double>* lprior_, c_select_interp& csel_, bool bcf_) {
    data = data_;
    coefs = coefs_;
    pi = pi_;
    di = di_;
    lprior = lprior_;
    csel = csel_;
    bcf = bcf_;
  }
  double val(double y) { //y is the length-scale parameters, a.k.a. l
    //Rcpp::Rcout << "inside function" << endl;
    //Rcpp::Rcout << "Pi: " << PI << endl;
    //Rcpp::Rcout << "Basis dimension: " << basis_dim << endl;
    //Rcpp::Rcout << "l sampled: " << y << endl;
    //Rcpp::Rcout << "Q: " << Q[0] << endl;
    //Rcpp::Rcout << "Sigma: " << sigma << endl;
    //Rcpp::Rcout << "Gamma A: " << gamma_a << endl;
    //Rcpp::Rcout << "Gamma B: " << gamma_b << endl;
    //double L = di.tlen/2 + std::max(di.tlen*0.025, 2*y); // L by using Murray's Law
    //double L = di.L_a + di.L_b*y;
    //double L = std::max((di.tlen/2)*(di.L_a + di.L_b*y),1.0);
    double L = (di.tlen/2)*csel.val(y);
    std::vector<double> Q(di.basis_dim*di.n,0); //allocate outside

    if(bcf){
      for(size_t i=0; i<di.n; i++){
        for(size_t j=0; j<di.basis_dim; j++){
          di.phi[(i*(di.basis_dim) + j)] = sqrt(1/L)*sin(PI*(j+1)*(di.t[i]+L)/(2*L))*di.z[i];
        }
      }
    }else{
      for(size_t i=0; i<di.n; i++){
        for(size_t j=0; j<di.basis_dim; j++){
          di.phi[(i*(di.basis_dim) + j)] = sqrt(1/L)*sin(PI*(j+1)*(di.t[i]+L)/(2*L));
        }
      }
    }

    for(size_t j=0;j<(di.n*di.basis_dim);++j){
      Q[j] = di.phi[j]*pi.eta*coefs[j]; //coefs as argument and join on the lklh
    }


    std::vector<double> lambda_vec(di.basis_dim,0); //creating the lambda vector with the value of y provided
    for(size_t k=0; k<di.basis_dim; ++k){
      lambda_vec[k] = sqrt(sqrt(2*PI)*y*exp(-y*y*((PI*(k+1)/(2*L))*(PI*(k+1)/(2*L)))/2)); //calculate the vector using the new value of l
      //Rcpp::Rcout << "Lambda vector: " << lambda_vec[k] << endl;
    }

    double ret=0;
    if((*lprior)[0]==1.){//Half-Cauchy
      double param_a, param_b; //hyperparameters of p(l)
      param_b = (*lprior)[1];

      ret = log(2.)-log(PI)-log(param_b)-log(1.+pow(y/param_b,2.));
    }else if((*lprior)[0]==2){//Beta-Prime
      double param_a, param_b; //hyperparameters of p(l)
      double aux;
      param_a = (*lprior)[1];
      param_b = (*lprior)[2];
      aux = R::lbeta(param_a,param_b);

      ret = (param_a-1)*log(y)+(-param_a-param_b)*log(1+y)-aux;
    }else if((*lprior)[0]==3){//Inverse-Gamma
      double param_a, param_b; //hyperparameters of p(l)
      param_a = (*lprior)[1];
      param_b = (*lprior)[2];

      ret = param_a*log(param_b)- lgamma(param_a) + (-param_a-1)*log(y)-param_b/y;
    }else if((*lprior)[0]==4){//Weibull
      double param_a, param_b; //hyperparameters of p(l)
      param_a = (*lprior)[1];
      param_b = (*lprior)[2];

      ret = R::dweibull(y, param_a, param_b, 1);
    }else if((*lprior)[0]==5){//Half-Normal
      double param_a, param_b; //hyperparameters of p(l)
      param_b = (*lprior)[1];

      ret = (1/2)*log(2)-log(param_b)-(1/2)*log(PI)-pow(y,2)/(2*pow(param_b,2));
    }else if((*lprior)[0]==6){//Log-Normal
      double param_a, param_b; //hyperparameters of p(l)
      param_a = (*lprior)[1];
      param_b = (*lprior)[2];

      ret = -log(y)-log(param_b)-(1/2)*log(2*PI)-pow((log(y)-param_a),2)/(2*pow(param_b,2));
    }else if((*lprior)[0]==7){//Weibull-Geometric Distribution
      double param_a, param_b, param_c; //hyperparameters of p(l)
      param_a = (*lprior)[1];
      param_b = (*lprior)[2];
      param_c = (*lprior)[3];

      ret = log(param_a)-param_a*log(param_b)+log(1-param_c)+(param_a-1)*log(y)-pow((y/param_b),param_a)-2*log(1-param_c*exp(-pow((y/param_b),param_a)));
    }else{//Gamma
      double param_a, param_b; //hyperparameters of p(l)
      param_a = (*lprior)[1];
      param_b = (*lprior)[2];

      ret = R::dgamma(y, param_a, param_b, 1);
    }

    //Rcpp::Rcout << "Prior d: " << ret << endl;
    double totsum = 0.0;
    double resid;

    for(size_t i=0; i<di.n; ++i) { //size of the mean vector in th multivariate normal
      double indsum = 0.0; //will be the sum of each term in the mean vector

      double term;
      for(size_t j=0; j<di.basis_dim; ++j) { //size of the basis functions
        term = Q[i*di.basis_dim + j] * lambda_vec[j]; //computing lambda*phi*beta
        indsum += term; //summing the j terms that make the basis function approx.
      }

      resid = data[i] - indsum; //taking the difference between the observed and the term of the mean vector
      totsum += resid*resid; //square the term and sum with the other terms (log prior + other terms)
    }
    ret -= totsum/(2*pi.sigma*pi.sigma);

    return(ret);
  }
};

class ld_dart_alpha: public logdensity {
  public:
  std::vector<double> counts;
  std::vector<double> var_sizes;
  std::vector<double> probs;
  double dart_a, dart_b, dart_rho;
  ld_dart_alpha(std::vector<double> counts_, std::vector<double> var_sizes_, std::vector<double> probs_) {
    ld_dart_alpha(counts_, var_sizes_, probs_, dart_a = 0.5, dart_b = 1, dart_rho = counts_.size());
  }
  ld_dart_alpha(std::vector<double> counts_, std::vector<double> var_sizes_, std::vector<double> probs_,
                double dart_a_, double dart_b_, double dart_rho_) {
    counts = counts_;
    var_sizes = var_sizes_;
    probs  = probs_;
    dart_a = dart_a_;
    dart_b = dart_b_;
    dart_rho = dart_rho_;
  }
  double val(double y) {
    double ret = R::dbeta(y/(y+dart_rho), dart_a, dart_b, 1) + log(dart_rho) - 2*log(y+dart_rho); //log prior
    double csum = 0.0;
    for(size_t i=0; i<counts.size(); ++i) {
      double totct = counts[i] + y/(counts.size()*var_sizes[i]);
      csum += totct;
      ret += (totct-1)*log(probs[i]) - lgamma(totct);
    }
    ret += lgamma(csum);
    return(ret);
  }
};


double slice(double x0, logdensity* g, double w=1., double m=INFINITY,
             double lower=-INFINITY, double upper=INFINITY);
#endif slice_h
