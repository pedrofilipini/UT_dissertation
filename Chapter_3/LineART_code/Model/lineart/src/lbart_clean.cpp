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

// y = m(x) + e, e~N(0, sigma^2_y)

//x_con is the design matrix for m. It should have n = rows

// [[Rcpp::export]]
List lbartRcpp(arma::vec y_, //response
               arma::mat Omega_con, //basis matrix
               arma::mat Omega_con_est, //for prediction
               NumericVector x_con_, //covariates
               NumericVector x_con_est_, //for prediction
               List x_con_info_list, //cutpoints
               int burn, int nd, int thin, //Draw nd*thin + burn samples, saving nd draws after burn-in
               int ntree_con, //number of tress
               double lambda, double nu, //prior pars for sigma^2_y
               arma::mat Sigma0_con, //variance estimate
               double con_alpha, double con_beta, //hyperparameters for alpha and beta
               CharacterVector treef_name_, //directory where I can save the trees
               bool prior_sample = false, // I DO NOT KNOW WHAT THIS MEANS
               bool use_con_scale = true, // Half-Cauchy prior for mu(x)
               double con_scale_df = 1, // I do not know what this thing do
               int status_interval=100) //printing intervals
{

  //for saving the trees
  std::string treef_name = as<std::string>(treef_name_);
  std::ofstream treef(treef_name.c_str());

  //RNG stuff
  RNGScope scope;
  RNG gen; //this one random number generator is used in all draws

  /*****************************************************************************
    /* Read, format y
  *****************************************************************************/
    std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY; //setting an initial range
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.

  //finding the min and the max of y
  for(NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) { //run through y_
    y.push_back(*it); //copy the whole vector on y (*it access y_)
    if(*it<miny) miny=*it; //compare the actual value with the actual min. Update the min.
    if(*it>maxy) maxy=*it; //compare the actual value with the actual max. Update the max.
    allys.sy += *it; // sum of y
    allys.sy2 += (*it)*(*it); // sum of y^2
  }
  size_t n = y.size(); //length of y
  allys.n = n; //include n among the suff stats

  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation


  //Now let us do the same thing for the covariate matrix
  /*****************************************************************************
    /* Read, format X_con
  *****************************************************************************/
    //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  std::vector<double> x_con; //storage for x_con_
  for(NumericVector::iterator it=x_con_.begin(); it!= x_con_.end(); ++it) { //run through x_con_
    x_con.push_back(*it); //copy the whole vector on x_con (*it access x_con_)
  }
  size_t p_con = x_con.size()/n; //since you know n from before, dividing the length of x_con by n will give you p (# covariates)

    Rcout << "Using " << p_con << " control variables." << std::endl; //printing information

    //x cutpoints
    xinfo xi_con; //creating variable for cutpoints. vector of vectors, will be split rules (defined on tree.h)

    xi_con.resize(p_con); //letting xi_con having size p
    for(int i=0; i<p_con; ++i) { //for each variable available
      NumericVector tmp = x_con_info_list[i]; //x_con_info_list is the original. We copy the file on tmp for variable [i]
      std::vector<double> tmp2; //creating another auxiliary vector
      for(size_t j=0; j<tmp.size(); ++j) { //for the size of this auxiliary vector
        tmp2.push_back(tmp[j]); //copy the content on the new vector (is this really necessary?)
      }
      xi_con[i] = tmp2; //include this vector in the vector of vectors xi_con[i] access the cutpoints of variable [i]
    }





    /*****************************************************************************
      /* Setup the model
    *****************************************************************************/
      //--------------------------------------------------
      //trees

    //Initialize beta in a naive way (using y)
    arma::vec betahat = zeros(2);//solve(Omega_con*Omega_con.t()+ //Omega^T Omega
                                         //0.05*eye(Omega_con.n_rows,Omega_con.n_rows), //Adding noise (assure PSD)
                                         //Omega_con*y_); //Omega^T y
    //Rcout<<betahat<<endl;

    std::vector<tree> t_con(ntree_con); //vector of trees, size ntree_con

    for(size_t i=0;i<ntree_con;i++) t_con[i].setm(zeros(2)); //setting parameters on every node (as zero)
    for(size_t i=0;i<ntree_con;i++) t_con[i].setv(floor((p_con)*(i)/(ntree_con)));

    //--------------------------------------------------
      //prior parameters

    pinfo pi_con; //object that holds the prior information
    pi_con.pbd = 1.0; //prob of birth/death move
    pi_con.pb = .5; //prob of birth given  birth/death

    pi_con.alpha = con_alpha; //include alpha received from function
    pi_con.beta  = con_beta; //include beta received from function

    //pi_con.mu0 = zeros(Omega_con.n_rows); //prior for mu0 (betas) is 0
    pi_con.mu0 = zeros(2); //prior for mu0 (betas) is 0
    pi_con.Sigma0 = Sigma0_con; //prior for variance of mu0 (betas)
    pi_con.Prec0 = pi_con.Sigma0.i(); //inverse of variance is precision
    pi_con.logdetSigma0 = log(det(pi_con.Sigma0)); //log of determinant of sigma prior
    pi_con.eta = 1;
    pi_con.gamma = 1;
    pi_con.scale_df = con_scale_df;
    pi_con.ntree = ntree_con;


    //Rcout << "Precision " << pi_con.Prec0 << endl;

    //Vector of mis for the precision
    //pi_con.mi = zeros(p_con); //initializing
    //size_t v_aux;
    //Filling the mis
    //for(size_t i=0;i<ntree_con;i++) {
    //  v_aux = t_con[i].getv();
    //  pi_con.mi[v_aux] += 1;
    //}

    //Rcout << pi_con.mi << endl;


    pi_con.sigma = shat; //sample sd

    double sigma = shat; //sample sd is the initial guess for sigma

    //--------------------------------------------------
      //dinfo for control function m(x)
    double* allfit_con = new double[n]; //sum of fit of all trees

    //TESTING
    for(size_t i=0;i<n;i++) allfit_con[i] = 0.0;// ybar;
    //END TESTING

    double* r_con = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp, will be the residual
    dinfo di_con; //data information
    di_con.n=n; //sample size
    di_con.p=p_con; //number of covariates
    di_con.x = &x_con[0]; //address of the covariates
    di_con.y=r_con; //the y for each draw will be the residual
    //di_con.basis_dim = Omega_con.n_rows; //dimension of basis
    di_con.basis_dim = 2; //dimension of basis
    di_con.omega = &Omega_con[0]; //address of the basis
    di_con.ntree = ntree_con; //number of trees



    dinfo di_con_est; //data information for prediction
    std::vector<double> x_con_est;     //stored like x
    size_t n_con_est;

    for(NumericVector::iterator it=x_con_est_.begin(); it!=x_con_est_.end(); ++it) {
      x_con_est.push_back(*it); //copy the whole vector (*it access the values)
    }
    n_con_est = x_con_est.size()/p_con; //dividing by p, so it is the size of n

    if(x_con_est.size() != n_con_est*p_con) stop("error, wrong number of elements in effect estimate data set\n");

    di_con_est.n=n_con_est; //size of n
    di_con_est.p=p_con; //size of p
    di_con_est.x = &x_con_est[0]; //address of matrix for prediction
    di_con_est.y=0; //there are no y's!
    di_con_est.ntree = ntree_con; //number of trees
    //di_con_est.basis_dim = Omega_con.n_rows; //size of basis
    di_con_est.basis_dim = 2; //size of basis
    di_con_est.omega = &Omega_con_est[0]; //placeholder, wopn't be touched
    //  }


//--------------------------------------------------
  //storage for ouput


//--------------------------------------------------
  //storage for the fits
double* allfit = new double[n]; //yhat total
for(size_t i=0;i<n;i++) { //run through the observations
  allfit[i] = allfit_con[i]; //only the mu(x) tree
}
double* ftemp  = new double[n]; //fit of current tree

// init node pointers
std::vector<std::vector<tree::tree_cp> > node_pointers_con(ntree_con); //node root pointers for each tree

for(size_t j=0; j<ntree_con; ++j) { //for each tree
  node_pointers_con[j].resize(n);
  fit_linear(t_con[j],xi_con,di_con,ftemp,node_pointers_con[j],true);
}

NumericVector sigma_post(nd);
NumericVector msd_post(nd);
NumericMatrix m_post(nd,n);
NumericMatrix yhat_post(nd,n);
NumericMatrix m_est_post(nd,n_con_est);
NumericMatrix yhat_est_post(nd,n_con_est);

arma::cube scoefs_con_est(di_con_est.basis_dim, di_con_est.n, nd);
arma::mat coefs_con_est(di_con_est.basis_dim, di_con_est.n);


/*
  //save stuff to tree file
treef << xi << endl; //cutpoints
treef << m << endl;  //number of trees
treef << p << endl;  //dimension of x's
     treef << (int)(nd/thin) << endl;
     */
std::stringstream treess;  //string stream to write trees to
treess.precision(10);
treess << nd << " " << ntree_con << " " << p_con << endl;

    //*****************************************************************************
    /* MCMC
     * note: the allfit objects are all carrying the appropriate scales
     */
    //*****************************************************************************
    Rcout << "\nBeginning MCMC:\n";
    time_t tp;
    int time1 = time(&tp);

    size_t save_ctr = 0;
    for(size_t i=0;i<(nd*thin+burn);i++) { //for the number of iterations

       if(prior_sample) {
          for(int k=0; k<n; k++) y[k] = gen.normal(allfit[k], sigma);
       }


       Rcpp::checkUserInterrupt();
       if(i%status_interval==0) {
          Rcout << "iteration: " << i << " sigma: "<< sigma << endl;
       }

       //draw trees for m(x)
       for(size_t j=0;j<ntree_con;j++) {

          fit_linear(t_con[j],xi_con,di_con,ftemp,node_pointers_con[j],false);
          for(size_t k=0;k<n;k++) {
             if(ftemp[k] != ftemp[k]) {
                Rcout << "control tree " << j <<" obs "<< k<<" "<< endl;
                Rcout << t_con[j] << endl;
                stop("nan in ftemp");
             }
             allfit[k] = allfit[k]-pi_con.eta*ftemp[k];
             allfit_con[k] = allfit_con[k]-pi_con.eta*ftemp[k];
             r_con[k] = (y[k]-allfit[k])/pi_con.eta;
             if(r_con[k] != r_con[k]) {
                Rcout << (y[k]-allfit[k]) << endl;
                Rcout << pi_con.eta << endl;
                Rcout << r_con[k] << endl;
                stop("NaN in resid");
             }
          }
          double aa = bd_linear(t_con[j],xi_con,di_con,pi_con,gen,node_pointers_con[j]);
          drmu_linear(t_con[j],xi_con,di_con,pi_con,gen);
          fit_linear(t_con[j],xi_con,di_con,ftemp,node_pointers_con[j],false);

          for(size_t k=0;k<n;k++) {
             allfit[k] += pi_con.eta*ftemp[k];
             allfit_con[k] += pi_con.eta*ftemp[k];
          }
       }


       for(size_t k=0;k<n;k++) {
          ftemp[k] = y[k] - allfit[k];
       }



       pi_con.sigma = sigma/fabs(pi_con.eta); //take a look


      //draw sigma
      double rss = 0.0;
      double restemp = 0.0;
      for(size_t k=0;k<n;k++) {
         restemp = y[k]-allfit[k];
         rss += restemp*restemp;
      }
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
      pi_con.sigma = sigma/fabs(pi_con.eta);

      //Rcout << sigma << endl;

      if( ((i>=burn) & (i % thin==0)) )  {
        for(size_t j=0;j<ntree_con;j++) treess << t_con[j]; //saving trees in a file

         sigma_post(save_ctr) = sigma;

         for(size_t k=0;k<n;k++) {
            m_post(save_ctr, k) = allfit_con[k];
            yhat_post(save_ctr, k) = allfit[k];
         }

         for(size_t k=0;k<di_con_est.n;k++) {

           m_est_post(save_ctr, k) = fit_i_linear(k, t_con, xi_con, di_con_est);
         }


         coefs_con_est.zeros();
         for(size_t j=0; j<ntree_con; ++j) {
            coefs_con_est += pi_con.eta*coef_basis(t_con[j], xi_con, di_con_est);
         }
         scoefs_con_est.slice(save_ctr) = coefs_con_est;

         save_ctr += 1;
      }
  }

  int time2 = time(&tp);
  Rcout << "time for loop: " << time2 - time1 << endl;

  t_con.clear();
  delete[] allfit;
  delete[] allfit_con;
  delete[] r_con;
  delete[] ftemp;

  treef.close();

  //Save all cutpoints in a list
  Rcpp::List xiret(xi_con.size());
  for(size_t i=0;i<xi_con.size();i++) {
    Rcpp::NumericVector vtemp(xi_con[i].size());
    std::copy(xi_con[i].begin(),xi_con[i].end(),vtemp.begin());
    xiret[i] = Rcpp::NumericVector(vtemp);
  }
  //Save all trees in a list
  Rcpp::List treesL;
  treesL["cutpoints"] = xiret;
  treesL["trees"]=Rcpp::CharacterVector(treess.str());

  return(List::create(_["yhat_post"] = yhat_post,
                      _["coefs_con_est"] = scoefs_con_est,
                      _["m_post"] = m_post,
                      _["m_est_post"] = m_est_post,
                      _["sigma"] = sigma_post,
                      _["treedraws"] = treesL
  ));
}
