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

using namespace Rcpp;


// [[Rcpp::export]]
List gpCore(arma::vec y_, arma::mat Omega_con,
             arma::mat Phi_con, double l_con, double L,
             int burn, int nd, int thin, double nu, double lambda,
             arma::mat Sigma0_con)
{

  RNGScope scope;
  RNG gen; //this one random number generator is used in all draws

  /*****************************************************************************/
  /* Read, format y */
  /*****************************************************************************/
  std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.

  for(NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it; // sum of y
    allys.sy2 += (*it)*(*it); // sum of y^2
  }
  size_t n = y.size();
  allys.n = n;

  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation
  //Rcout << "One" << std::endl;

  //--------------------------------------------------
  //prior parameters


  //Rcout << "Two" << std::endl;
  pinfo pi_con;
  pi_con.mu0 = zeros(Omega_con.n_rows);
  pi_con.Sigma0 = Sigma0_con;
  pi_con.Prec0 = pi_con.Sigma0.i();
  pi_con.logdetSigma0 = log(det(pi_con.Sigma0));
  pi_con.eta = 1;
  pi_con.gamma = 1;

  pi_con.sigma = 1;//shat;
  pi_con.scale_df = 1;

  pi_con.length_scale = l_con;
  pi_con.L = L;

  double sigma = 1;//shat;
  //Rcout << "Three" << std::endl;

  //--------------------------------------------------
  //dinfo for control function m(x)
  //Rcout << "ybar " << ybar << endl;
  double* allfit_con = new double[n]; //sum of fit of all trees
  //TESTING
  //Rcout << "Three One" << std::endl;
  for(size_t i=0;i<n;i++) allfit_con[i] = ybar;// ybar;
  //Rcout << "Three Two" << std::endl;
  //END TESTING
  dinfo di_con;
  di_con.n=n;
  double* y_data = new double[n];
  for(size_t k=0; k<n; k++){
    y_data[k]=y[k];
  }
  di_con.y=y_data; //
  //Rcout << "Three Three" << std::endl;
  di_con.basis_dim = Omega_con.n_rows;
  di_con.omega = &Omega_con[0];
  di_con.phi = &Phi_con[0];

  //Rcout << "Four" << std::endl;

  //--------------------------------------------------
  //storage for the fits
  double* allfit = new double[n]; //yhat
  for(size_t i=0;i<n;i++) {
    allfit[i] = ybar;
  }
  double* ftemp  = new double[n]; //fit of current tree

  //Rcout << "Eight" << std::endl;
  NumericVector l_con_post(nd);
  NumericVector sigma_post(nd);
  NumericMatrix m_post(nd,di_con.basis_dim);
  NumericVector eta_con_post(nd);

  double* Q_con  = new double[(di_con.basis_dim*di_con.n)];

  //Rcout << "Basis size: " << di_con.basis_dim  << endl;

  /*****************************************************************************/
  /* MCMC */
  /*****************************************************************************/

  Rcout << "\nBeginning MCMC:\n";
  time_t tp;
  int time1 = time(&tp);

  size_t save_ctr = 0;

  double* beta_post = new double[(di_con.basis_dim)];


  for(size_t i=0;i<(nd*thin+burn);i++) {
    Rcout << "Iteration: " << i << endl;

    //vector for suff stats
    std::vector<sinfo> sv;

    //get sufficient statistics
    allsuff_basis_gp(di_con, sv, pi_con.eta);

    //sample mu
    drmu_basis_gp(beta_post, di_con,pi_con, sv, gen);

    //Rcout << "Sampled betas" << endl;

    //for(size_t j=0; j<di_con.basis_dim; ++j) {
    //  Rcout << "beta[" << j << "]" << beta_post[j] <<endl;
    //}

    //Constructing Q with current beta and eta
    for(size_t k=0; k<n; k++){
      for(size_t j=0;j<(di_con.basis_dim);++j){
        Q_con[k*di_con.basis_dim + j] = di_con.phi[k*di_con.basis_dim + j]*pi_con.eta*beta_post[j];
      }
    }

    //for(size_t j=0; j<di_con.basis_dim; ++j) {
    //  Rcout << "Phi[" << j << "]" << di_con.phi[j] <<endl;
    //}

    //for(size_t j=0; j<2*di_con.basis_dim; ++j) {
    //  Rcout << "Q[" << j << "]" << Q_con[j] <<endl;
    //}

    //Rcout << "Constructed Q" << endl;

    //creating allfit_con
    double *lambda_vec = new double[di_con.basis_dim]; //creating the lambda vector with the value of y provided
    for(size_t k=0; k<di_con.basis_dim; ++k){
      lambda_vec[k] = sqrt(sqrt(2*PI)*pi_con.length_scale*exp(-0.5*pi_con.length_scale*pi_con.length_scale*(PI*(k+1)/(2*L))*(PI*(k+1)/(2*L)))); //calculate the vector using the new value of l
      //Rcpp::Rcout << "Lambda vector: " << lambda_vec[k] << endl;
    }


    for(size_t k=0; k<n; ++k) { //size of the mean vector in th multivariate normal
      double indsum = 0.0; //will be the sum of each term in the mean vector

      double term;
      for(size_t j=0; j<di_con.basis_dim; ++j) { //size of the basis functions
        term = Q_con[k*di_con.basis_dim + j] * lambda_vec[j]; //computing lambda*phi*beta
        indsum += term; //summing the j terms that make the basis function approx.
      }
      allfit_con[k] = indsum;
      //Rcout << "Allfit[" << k << "]: " << allfit_con[k] << endl;
      allfit[k] = allfit_con[k];
    }
    //getting sigma back
    pi_con.sigma = sigma;
    //Rcout << "Before l sample" << endl;
    //update_length_scale(y_data, allfit_con, n, Q_con, pi_con, di_con);

    //Rcout << "After l sample" << endl;

    for(size_t k=0; k<n; ++k) { //size of the mean vector in th multivariate normal
      //Rcout << "Allfit[" << k << "]: " << allfit_con[k] << endl;
      allfit[k] = allfit_con[k];
    }




    //Update the eta
    double eta_old;

    for(size_t k=0;k<n;k++) {
      ftemp[k] = y[k];
    }
    pi_con.sigma = sigma/fabs(pi_con.eta);
    eta_old = pi_con.eta; // Save previous eta before drawing new one, for adjusting scaling.
    pi_con.scale_df=0;
    update_scale(ftemp, allfit_con, n, sigma, pi_con, gen);

    //Rcout << "Eta" << endl;

    // Update fits to have new pi.eta scaling.
    for(size_t k=0; k<n; ++k) {
      allfit[k] -= allfit_con[k];
      allfit_con[k] = allfit_con[k] * pi_con.eta / eta_old;
      allfit[k] += allfit_con[k];
    }
    //
    //Rcout << pi_mod.eta << endl;

    //draw sigma
    double rss = 0.0;
    double restemp = 0.0;
    for(size_t k=0;k<n;k++) {
      restemp = y[k]-allfit[k];
      rss += restemp*restemp;
    }
    sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));

    pi_con.sigma = sigma/fabs(pi_con.eta);
    //pi_con.sigma = pi_con.sigma*eta_old/pi_con.eta;


    //place where I am saving the betas
    for(size_t k=0;k<di_con.basis_dim;k++) {
      m_post(save_ctr, k) = pi_con.eta*beta_post[k];
    }
    //Rcout << "Saved beta" << endl;
    sigma_post(save_ctr) = sigma;
    l_con_post(save_ctr) = pi_con.length_scale;
    eta_con_post(save_ctr) = pi_con.eta;

    save_ctr += 1;

    //Rcout << "Finished iteration " << i << endl;
  }

  //Rcout << "Before time " << endl;
  int time2 = time(&tp);
  Rcout << "Time for loop: " << time2 - time1 << endl;


  return(List::create(_["m_post"] = m_post,
                      _["eta_con_post"] = eta_con_post,
                      _["sigma"] = sigma_post,
                      _["l_con_post"] = l_con_post
  ));
}
