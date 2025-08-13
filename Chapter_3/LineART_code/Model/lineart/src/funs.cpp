#include "arma_config.h"
#include <RcppArmadillo.h>

#include <cmath>
#include "funs.h"
#include "rng.h"
#include <map>
#ifdef MPIBART
#include "mpi.h"
#endif

#include "slice.h"

using Rcpp::Rcout;
using namespace arma;
using namespace Rcpp;
#define PI 3.1415926535897931

//find variables n can split on, put their indexes in goodvars
void getbadvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& badvars)
{
  int L,U;
  for(size_t v=0;v!=xi.size();v++) {//try each variable
    L=0;
    U = xi[v].size()-1; //checking number of elements for each variable
    n->rg(v,&L,&U); //find region [L,U] for var v.
    if(U<L) badvars.push_back(v); //no valid region between [L,U] for var v.
  }
}

void update_dart(std::vector<tree>& trees, pinfo& pi, dinfo& di, xinfo& xi, RNG gen) {

  std::vector<double> var_use_counts(di.p, 0); //counting how many times a variable has been used

  for(size_t j=0;j<trees.size();j++) { //for each tree
    std::vector<tree*> tj_nodes;
    trees[j].getnobots(tj_nodes); //get all nodes except leaves


    for(size_t b=0; b<tj_nodes.size(); ++b) {
      var_use_counts[tj_nodes[b]->getv()] += 1.0; //counting on how many nodes a variable has been used
      std::vector<size_t> badvars; //variables this tree couldn't split on
      getbadvars(tj_nodes[b], xi, badvars);

      if(badvars.size()>0) {
        double total_prob = 0;
        vector<double> sample_probs(badvars.size()); //vector of probs for badvars

        for(size_t r = 0; r<badvars.size(); ++r) {
          total_prob += pi.var_probs[badvars[r]]; //taking the probs that have been sampled for badvars and summing
          sample_probs[r] = pi.var_probs[badvars[r]]; //taking the probs that have been sampled for badvars
        }

        for(size_t r = 0; r<badvars.size(); ++r) sample_probs[r]/=total_prob; //make between 0 and 1
        double numaug = R::rnbinom(1, 1-total_prob); //binomial with p=probability of goodvars

        int tt = 0;
        while( (numaug>0.5) & (tt<50000) ) { //just to make sure it does not get stuck in here
          var_use_counts[badvars[rdisc(&sample_probs[0], gen)]] += 1;
          numaug -= 1; tt += 1;
        }

      }

    }
  }

  //Rcpp::Rcout << "Out of first for" << endl;

  //Rcpp::Rcout << "di.p: " << di.p <<endl;

  double tot = 0;

  for(size_t v = 0; v < di.p; ++v) {
    //Rcpp::Rcout << "pi.var_sizes[v]: " << pi.var_sizes[v] <<endl;
    //Rcpp::Rcout << "pi.alpha: " << pi.alpha <<endl;
    //Rcpp::Rcout << "di.p: " << di.p <<endl;

    //Rcpp::Rcout << "var_use_counts[v]: " << var_use_counts[v] <<endl;
    pi.var_probs[v] = gen.gamma(var_use_counts[v] + pi.alpha/(di.p * pi.var_sizes[v]), 1); //sampling dirichlet by using Gammas (see wikipedia)
    //Rcpp::Rcout << "pi.var_probs: " << pi.var_probs[v] <<endl;
    tot += pi.var_probs[v];
  }

  //Rcpp::Rcout << "Out of second for" << endl;

  for(size_t v = 0; v < di.p; ++v) pi.var_probs[v] /= tot; //just constructing in form of probability for the Dirichlet sample

  //New class defined at slice.h
  ld_dart_alpha logdens(var_use_counts, pi.var_sizes, pi.var_probs); //receive counts, varsizes and probabilities

  double a0 = pi.dart_alpha;

  //Rcpp::Rcout << "Before slice" << endl;
  pi.dart_alpha = slice(a0, &logdens, 1.0, INFINITY, 0.0, INFINITY); //sample new alpha through a slice sampler
  //Rcpp::Rcout << "after slice" << endl;
}

void update_scale(double* r, double* fits, size_t n, double sigma, pinfo& pi, RNG& gen) {

  // assumes that the scale multiplier is contained in fits!

  double a = 0.0; double b = 0.0;
  //double sigma = pi.sigma*fabs(pi.eta);
  for(size_t k=0;k<n;k++) {
    //double r = y[k] - allfit_con[k];
    //if(randeff) r -= allfit_random[k];

    a += fits[k] * r[k]/pi.eta;
    b += fits[k] * fits[k]/(pi.eta*pi.eta);
  }

  double var = 1 / (1 / (pi.gamma * pi.gamma) + b/(sigma*sigma));
  double mean = var * (0.0/(pi.gamma*pi.gamma) + a/(sigma*sigma));

  // Rcpp::Rcout << "postmean " << mean << std::endl;
  // Rcpp::Rcout << "var " << var << std::endl;
  
  pi.eta = mean + gen.normal(0.,1.) * sqrt(var);

  if(pi.scale_df>0) {
    pi.gamma = sqrt(0.5*(pi.scale_df+pi.eta*pi.eta)/gen.gamma((1.0 + pi.scale_df)/2.0, 1.0));
  }

}

//-------------------------------------------------------------
// Squared Exponential Covariance Function for two time vectors (x,y).
//-------------------------------------------------------------
mat cov_se(vec t1, vec t2, double ls, double var)
{
   //-------------------------------------------------------------
   // INPUTS:	   x,y   = two vectors from the same space.
   //				ls    = b 		= length (of period)
   //				var   = tau1.sq = variance of function
   //-------------------------------------------------------------
   // OUTPUT:	The squared exponential covariance matrix.
   //-------------------------------------------------------------
   double arg;
   int n1 = t1.size();
   int n2 = t2.size();
   mat C(n1,n2);
   for(int i = 0; i < n1; i++) {
      for(int j=0; j < n2; j++) {
         arg = (t1[i] - t2[j])/ls;
         C(i,j) = var*exp(-0.5*arg*arg);
      }
   }

   // Add jitter to diagonal to ensure is pos semidef.
   C.diag() = C.diag() + .000001;

   return C;
}

//-------------------------------------------------------------
// Utility function for calculating posterior MVN params for N(Phi^(-1)*m, Phi^(-1))
//-------------------------------------------------------------
List mvn_post_util(double sigma, vec mu0, mat Prec0, vec n_vec, vec sy_vec){

   // Initialize matrices and vectors.
   double s2 = sigma * sigma;
   mat Lam = diagmat(n_vec) / s2;
   vec m = (1/s2) * sy_vec + Prec0 * mu0;  // K = Prec0

   return List::create(
      _["Phi"] = Prec0 + Lam,
      _["m"]   = m
   ) ;
}

List mvn_post_util_het(vec mu0, mat Prec0, vec n0_vec, vec n_vec, vec sy_vec){
   //n0_vec for het is vector of nt for each time (not weighted).
   //n_vec for het is vector of sum(w_i) precisions at each time.
   //sy_vec for het is vector of sum(y_i*w_i) at each time.

   // Initialize matrices and vectors.
   mat Lam = diagmat(n_vec) + Prec0;
   vec m = sy_vec + Prec0 * mu0;

   return List::create(
      _["Phi"] = Lam,
      _["m"]   = m
   ) ;
}

//-------------------------------------------------------------
// Generates realizations from multivariate normal.
//-------------------------------------------------------------
mat rmvnormArma(int n, vec mu, mat sigma) {
   //-------------------------------------------------------------
   // INPUTS:	   n = sample size
   //				   mu = vector of means
   //				   sigma = covariance matrix
   //-------------------------------------------------------------
   // OUTPUT:	n realizations of the specified MVN.
   //-------------------------------------------------------------
   int ncols = sigma.n_cols;
   mat Y = randn(n, ncols);
   mat result = (repmat(mu, 1, n).t() + Y * chol(sigma)).t();
   return result;
}

//-------------------------------------------------------------------------------
// log of the integrated likelihood for tsbart, for a given tree/leaf.
//-------------------------------------------------------------------------------
double lil_ts(vec nt, vec sy_vec, double sy2, double sigma, vec mu0, mat Prec0){
   // nt = vector of number of obs in each time point for the given tree/leaf. nt = [nl_{t=1}, ..., nl_{t=T}]
   // sy = vector of sums of y's at each time point.  sy = [ sum(y in t=1), ..., sum(y in t=T) ]
   // sigma = error sd, sqrt of sigma^2
   // mu0 = vector of prior means.
   // Prec0 = prior precision matrix for means (from sq exp kernel)

   // For computational efficiency, we leave out the -.5*(mu_0^T K mu_0) term, since
   // we let mu_0=0.  Add this term if mu0 != 0.

   double sig2 = sigma*sigma;    //sigma^2
   double nl = sum(nt);          // Total number of yl in leaf.

   // Precache a few terms to make calculation faster.
   vec b = sy_vec/sig2 + Prec0*mu0;
   mat C = Prec0;
   C.diag() = C.diag() + nt/sig2;

   // Calculate log-likelihood.  Note: mu0.t()*K*mu0 excluded as mu0=0.
   double ll = -.5*nl*log(2*PI*sig2) + .5*log(det(Prec0)) - .5*log(det(C)) -
                .5*as_scalar(sy2/sig2 - b.t()*C.i()*b);

   return(ll);
}

// For het variances.
double lilhet_ts(double n0, double n, vec n_vec, vec sy_vec, double sy2, vec mu0, mat Prec0){
   // n0 = number of obs in given tree/leaf.
   // n = sum of log-precisions phi in given tree/leaf.
   // n_vec = vector of sum of phi's for het, at each time point, for given tree/leaf.
   // sy = vector of sums of y's * phi's at each time point.  sy = [ sum(phi*y in t=1), ..., sum(phi*y in t=T) ]
   // sy2 = scalar, sum of y*y*phi for all obs in given tree/leaf.
   // mu0 = vector of prior means.
   // Prec0 = prior precision matrix for means (from sq exp kernel)

   // For computational efficiency, we leave out the -.5*(mu_0^T K mu_0) term, since
   // we let mu_0=0.  Add this term if mu0 != 0.

   // Precache a few terms to make calculation faster.
   mat C = Prec0; // K = Prec0
   C.diag() += n_vec; // Add sums of precisions to diagonal.
   vec b = sy_vec + Prec0 * mu0;

   // Calculate log-likelihood.  Note: mu0.t()*K*mu0 excluded as cancels in ratios.
   double ll = - .5*n0*log(2*PI)
               + .5*n // This is the .5 * log(det(Lambda)) term, where Lambda=diag(w).
               + .5*log(det(Prec0))
               - .5*log(det(C))
               - as_scalar(.5*(sy2 - b.t()*C.i()*b));

   return(ll);
}

//Basis

void allsuff_basis(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, std::vector<tree::tree_cp>& node_pointers)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
  // or long int, etc.  We don't have to know that ahead of time by using this notation.
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    //sv[i].n_vec.resize(di.tlen);
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
    y = di.y[i];           // Resolves to r.

    //tbn = x.bn(xx,xi); // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    // get the variable we split on similar to allsuff_basis_birth

    ++(sv[ni].n0);

    if(di.basis_dim == 1) {
      omega = di.omega[i];
      sv[ni].sy += omega*y;
      sv[ni].n += omega*omega;
    } else {

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      sv[ni].sy_vec += y*omega_i;

      for(size_t j=0; j<di.basis_dim; ++j) {
        sv[ni].WtW.at(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
        for(size_t g=0; g<j; ++g) {
          double a = omega_i_tmp[j]*omega_i_tmp[g];
          sv[ni].WtW.at(g,j) += a;
          //sv[ni].WtW(j,g) += a; //this is faster than below?
          //sv[ni].WtW[i,j] = sv[ni].WtW[i,j] + a;
          //sv[ni].WtW[j,i] = sv[ni].WtW[j,i] + a; //get this outside obs loop later
        }
      }
    }

  }

  if(di.basis_dim>1) {
    for(size_t q=0; q<sv.size(); ++q) {
      for(size_t j=0; j<di.basis_dim; ++j) {
        for(size_t g=0; g<j; ++g) {
          sv[q].WtW.at(j,g) = sv[q].WtW.at(g,j);
        }
      }
    }
  }

  if(di.basis_dim<2) {
    for(size_t j=0; j<sv.size(); ++j) {
      sv[j].WtW.at(0,0) = sv[j].n;
      sv[j].sy_vec.at(0) = sv[j].sy;
    }
  }


}

void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv,
                         std::vector<sinfo>& sv, sinfo& sl, sinfo& sr, std::vector<tree::tree_cp>& node_pointers)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    //sv[i].n_vec.resize(di.tlen);
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].n = 0.0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  double *omega_i_tmp;
  arma::vec omega_i;
  bool in_candidate_nog, left, right;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    y = di.y[i];

    //tbn = x.bn(xx,xi); // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    // Get the split variable by calling int split_var = (*tbn->getp()).getv()

    left = false; right = false;

    in_candidate_nog = (tbn == nx);

    //For the proposed split we know the variable is v so when we compute
    //summary stats for sl ad sr, use variable v

    left  = in_candidate_nog & (xx[v] < xi[v][c]);
    right = in_candidate_nog & !(xx[v] < xi[v][c]);

    ++(sv[ni].n0);
    if(left) sl.n0 += 1;
    if(right) sr.n0 += 1;

    if(di.basis_dim == 1) {
      omega = di.omega[i];
      sv[ni].sy += omega*y;
      sv[ni].n += omega*omega;
      if(left) {
        sl.sy += omega*y;
        sl.n += omega*omega;
      }
      if(right) {
        sr.sy += omega*y;
        sr.n += omega*omega;
      }

    } else {
      omega_i_tmp = di.omega + i*di.basis_dim;
      //get design vector
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      sv[ni].sy_vec += y*omega_i;
      if(left) sl.sy_vec += y*omega_i;
      if(right) sr.sy_vec += y*omega_i;

      for(size_t j=0; j<di.basis_dim; ++j) {
        double a = omega_i_tmp[j]*omega_i_tmp[j];
        sv[ni].WtW.at(j,j) += a;
        if(left) sl.WtW.at(j,j) += a;
        if(right) sr.WtW.at(j,j) += a;
        for(size_t g=0; g<j; ++g) {
          double a = omega_i_tmp[j]*omega_i_tmp[g];
          sv[ni].WtW.at(g,j) += a;
          if(left) sl.WtW.at(g,j) += a;
          if(right) sr.WtW.at(g,j) += a;
        }
      }
    }
  }

  if(di.basis_dim>1) {

    //Rcout << "L" << endl << sl.WtW << endl << sl.sy_vec << endl;
    //Rcout << "R" << endl << sr.WtW << endl << sr.sy_vec << endl;

    for(size_t q=0; q<sv.size(); ++q) {
      for(size_t j=0; j<di.basis_dim; ++j) {
        for(size_t g=0; g<j; ++g) {
          sv[q].WtW(j,g) = sv[q].WtW(g,j);
          sl.WtW(j,g) = sl.WtW(g,j);
          sr.WtW(j,g) = sr.WtW(g,j);
        }
      }
    }

    //Rcoutt << "L" << endl << sl.WtW << endl;
    //Rcoutt << "R" << endl << sr.WtW << endl;
  }

  if(di.basis_dim<2) {
    sl.WtW.at(0,0) = sl.n;
    sl.sy_vec.at(0) = sl.sy;
    sr.WtW.at(0,0) = sr.n;
    sr.sy_vec.at(0) = sr.sy;
    for(size_t j=0; j<sv.size(); ++j) {
      sv[j].WtW.at(0,0) = sv[j].n;
      sv[j].sy_vec.at(0) = sv[j].sy;
    }
  }




}


void getsuff_basis(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double *xx;//current x
  double y;  //current y
  double t;  //current t

  
  sl.n0=0;sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sl.n0=0;sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  sl.WtW = zeros(di.basis_dim,di.basis_dim); sl.sy_vec = zeros(di.basis_dim);
  sr.WtW = zeros(di.basis_dim,di.basis_dim); sr.sy_vec = zeros(di.basis_dim);

  bool orig = false;

  double omega;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;


    if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

      y = di.y[i];   // extract current yi.  resolves to r.

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      if(xx[v] < xi[v][c]) { // Update left.

        ++(sl.n0);

        if(di.basis_dim==1) {
          omega = di.omega[i];
          sl.sy += omega*y;
          sl.n += omega*omega;

        } else {

          sl.sy_vec += y*omega_i;
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sl.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sl.WtW(i,j) = sl.WtW(i,j) + a;
              //sl.WtW(j,i) += a;
            }
          }

        }


      } else { //Update right.

        ++(sr.n0);

        if(di.basis_dim==1) {
          omega = di.omega[i];
          sr.sy += omega*y;
          sr.n += omega*omega;

        } else {

          sr.sy_vec += y*omega_i;
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sr.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sr.WtW(i,j) = sr.WtW(i,j) + a;
              //sr.WtW(j,i) += a;
            }
          }

        }


      }
    }
  }


  if(di.basis_dim<2) {
    sl.WtW(0,0) = sl.n;
    sl.sy_vec[0] = sl.sy;
    sr.WtW(0,0) = sr.n;
    sr.sy_vec[0] = sr.sy;
  }

  // this is faster without bounds checking, but incrementing w/ +=above is faster *with* bounds checking?
  for(size_t j=0; j<di.basis_dim; ++j) {
    for(size_t i=0; i<j; ++i) {
      sl.WtW(j,i) = sl.WtW(i,j);
      sr.WtW(j,i) = sr.WtW(i,j);
    }
  }

}

/*
void getsuff_basis(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double *xx;//current x
  double y;  //current y
  double t;  //current t

  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  sl.WtW = zeros(di.basis_dim,di.basis_dim); sl.sy_vec = zeros(di.basis_dim);
  sr.WtW = zeros(di.basis_dim,di.basis_dim); sr.sy_vec = zeros(di.basis_dim);

  bool orig = false;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;


    if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

      y = di.y[i];   // extract current yi.  resolves to r.

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);


      if(xx[v] < xi[v][c]) { // Update left.

        ++(sl.n);
        sl.sy_vec += y*omega_i;

        if(orig) {
          sl.WtW += omega_i*omega_i.t();
        } else {
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sl.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sl.WtW[i,j] = sl.WtW[i,j] + a;
              //sl.WtW(j,i) += a;
            }
          }
        }



      } else { //Update right.
        ++(sr.n);
        sr.sy_vec += y*omega_i;

        if(orig) {
          sr.WtW += omega_i*omega_i.t();
        } else {
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sr.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sr.WtW[i,j] = sr.WtW[i,j] + a;
              //sr.WtW(j,i) += a;
            }
          }
        }

      }
    }
  }

  // this is faster without bounds checking, but incrementing w/ +=above is faster *with* bounds checking
  for(size_t j=0; j<di.basis_dim; ++j) {
    //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
    for(size_t i=0; i<j; ++i) {
      sl.WtW[j,i] = sl.WtW[i,j];
      sr.WtW[j,i] = sr.WtW[i,j];
    }
  }

}
*/
void getsuff_basis(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double *xx;//current x
  double y;  //current y
  double t;  //current t

  bool orig = false;

  sl.n0=0;sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sl.n0=0;sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  double omega;

  sl.WtW = zeros(di.basis_dim,di.basis_dim); sl.sy_vec = zeros(di.basis_dim);
  sr.WtW = zeros(di.basis_dim,di.basis_dim); sr.sy_vec = zeros(di.basis_dim);

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    tree::tree_cp bn = x.bn(xx,xi);

    y = di.y[i];   // extract current yi.

    if(bn==nl) {

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      ++(sl.n0);

      if(di.basis_dim==1) {
        omega = di.omega[i];
        sl.sy += omega*y;
        sl.n += omega*omega;

      } else {

        sl.sy_vec += y*omega_i;
        for(size_t j=0; j<di.basis_dim; ++j) {
          //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
          sl.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
          for(size_t i=0; i<j; ++i) {
            double a = omega_i_tmp[j]*omega_i_tmp[i];
            sl.WtW(i,j) = sl.WtW(i,j) + a;
            //sl.WtW(j,i) += a;
          }
        }

      }

    }

    if(bn==nr) {

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      ++(sr.n0);

      if(di.basis_dim==1) {
        omega = di.omega[i];
        sr.sy += omega*y;
        sr.n += omega*omega;

      } else {

        sr.sy_vec += y*omega_i;
        for(size_t j=0; j<di.basis_dim; ++j) {
          //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
          sr.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
          for(size_t i=0; i<j; ++i) {
            double a = omega_i_tmp[j]*omega_i_tmp[i];
            sr.WtW(i,j) = sr.WtW(i,j) + a;
            //sr.WtW(j,i) += a;
          }
        }

      }

    }
  }


  if(di.basis_dim<2) {
    sl.WtW(0,0) = sl.n0;
    sl.sy_vec[0] = sl.sy;
    sr.WtW(0,0) = sr.n0;
    sr.sy_vec[0] = sr.sy;
  }

  for(size_t j=0; j<di.basis_dim; ++j) {
    for(size_t i=0; i<j; ++i) {
      sl.WtW(j,i) = sl.WtW(i,j);
      sr.WtW(j,i) = sr.WtW(i,j);
    }
  }

}

void drmu_basis(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{

  bool debug = false;
  tree::npv bnv;
  std::vector<sinfo> sv; //will be resized in allsuff
  tree::npv bnv0; std::vector<sinfo> sv0;
  //debug is broken, would need to eat node_pointers. suff stat code works anyhow
//  if(debug) allsuff_basis(t,xi,di,bnv0,sv0,node_pointers);

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  t.getbots(bnv);   // Save bottom nodes for x to bnv variable.
  tree::npv::size_type nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  double prior_prec = pi.Prec0(0,0);
  double prior_mean = pi.mu0(0);

  for(tree::npv::size_type i=0; i<nb; i++) {
    sv[i] = bnv[i]->s;
    if(debug) {
      if(abs(sv[i].sy - sv0[i].sy)>1e-8) {

        Rcout << "depth "<< bnv[i]->depth() << endl;
        Rcout << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
        Rcout << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
        stop("shit");
      }
     else {
        //Rcout << "depth "<< bnv[i]->depth() << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec <<" new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << endl;
        //Rcoutt << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
       //Rcoutt << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
      }
     sv[i] = sv0[i];
    }

  }

  if(di.basis_dim<2) {
    vec beta_draw(1);
    double tt = prior_prec*prior_mean;
    double s2 = (pi.sigma*pi.sigma);
    
    for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      double Phi = 1.0/(sv[i].n/s2 + prior_prec);
      double m = tt + sv[i].sy_vec(0)/s2;
      beta_draw(0) = m*Phi + gen.normal(0,1)*sqrt(Phi); //rmvnorm_post(m, Phi);
      // Rcpp::Rcout<< "prior_prec " << prior_prec << endl;
      // Rcpp::Rcout<< "prior_mean " << prior_mean << endl;
      // Rcpp::Rcout<< "mu_a " << m/Phi << endl;
      // Rcpp::Rcout<< "Var_a " << (1/Phi) << endl;

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(beta_draw);
      
      
      // Rcout << "fcmean " << m*Phi << endl;
      // Rcout << "sqrt(fcvar) " << sqrt(Phi) << endl;

      // Check for NA result.
      if(beta_draw(0) != beta_draw(0)) {
        Rcpp::stop("drmu failed");
      }
    }
  } else {
    mat Phi;
    vec m;
    vec beta_draw;

    vec tt = pi.Prec0*pi.mu0;
    double s2 = (pi.sigma*pi.sigma);
    for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      //Rcoutt << "phi ";
      //Rcoutt << "WtW" << endl << sv[i].WtW << endl;
      Phi = sv[i].WtW/s2 + pi.Prec0;
      //Rcoutt << "m ";
      m = tt + sv[i].sy_vec/s2;
      //Rcoutt << "draw ";
      beta_draw = rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(beta_draw);

      // Check for NA result.
      if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
        Rcpp::stop("drmu failed");
      }
    }
  }

}


double lil_basis(sinfo& s, pinfo& pi){

//  Rcout << "log likelihood calc" << endl;

  double ll=0; double tt;

  double s2 = pi.sigma*pi.sigma;

  if(false) {//(s.sy_vec.n_elem == 1) {
    /*

    This is slow as hell??
     Instead: check for 1d when computing summary stats, avoid matrix ops** everywhere**

    double precpred = pi.Prec0.at(0,0) + s.n0/s2;

    //double dotp = dot(s.sy_vec/s2, solve(precpred, s.sy_vec/s2));
    double dotp = precpred*s.sy/(s2*s2);
    double tt = 0;

    tt = -0.5*log(precpred) - 0.5*pi.logdetSigma0;

    ll = -0.5*n*log(s2) + tt + 0.5*dotp;
    */


  } else {
    // double precpred = pi.Prec0(0,0) + s.n/s2;
    // 
    // double dotp = (s.sy/s2)*(s.sy/s2)/precpred;
    // double tt = 0;
    // 
    // tt = -0.5*log(precpred) - 0.5*pi.logdetSigma0;
    // 
    // ll = -0.5*s.n0*log(s2) + tt + 0.5*dotp;
    
    double tau = 1/pi.Prec0(0,0);
    double n = s.n;
    
    
    ll = 0.5*log(1/(1+tau*n/s2)) + 0.5*(tau/((1+tau*n/s2)))*(s.sy/s2)*(s.sy/s2);
    
    // Rcpp::Rcout<< "s2 " << s2 << endl;
    // Rcpp::Rcout<< "tau " << tau << endl;
    // Rcpp::Rcout<< "n " << n << endl;
    // Rcpp::Rcout<< "sy " << s.sy << endl;
    // Rcpp::Rcout<< "ll " << ll << endl;
    // 
    
    // mat precpred = pi.Prec0 + s.WtW/s2;
    // 
    // double dotp = dot(s.sy_vec/s2, solve(precpred, s.sy_vec/s2));
    // double tt = 0;
    // 
    // tt = -0.5*log(det(precpred)) - 0.5*pi.logdetSigma0;
    // 
    // ll = -0.5*s.n0*log(s2) + tt + 0.5*dotp;
    // Rcpp::Rcout<< "solve " << solve(precpred, s.sy_vec/s2) << endl;
    // Rcpp::Rcout<< "s.sy_vec " << s.sy_vec << endl;
    // Rcpp::Rcout<< "s.WtW " << s.WtW << endl;
    // Rcpp::Rcout<< "pi.Prec0 " << pi.Prec0 << endl;
    // Rcpp::Rcout<< "det(precpred) " << det(precpred) << endl;
    // Rcpp::Rcout<< "logdetSigma0 " << pi.logdetSigma0 << endl;
    // Rcpp::Rcout<< "s2 " << s2 << endl;
    // Rcpp::Rcout<< "dotp " << dotp << endl;
    // Rcpp::Rcout<< "ll " << ll << endl;
    
  }

  // Rcpp::Rcout<< "ll " << ll << endl;

  return(ll);
}

double lil_basis(sinfo& sl, sinfo& sr, pinfo& pi){

  sinfo st = sl;
  st.n += sr.n;
  st.n0 += sr.n0;
  st.sy += sr.sy;
  st.sy_vec += sr.sy_vec;
  st.WtW += sr.WtW;

  double ll = lil_basis(st, pi);

  //Rcpp::Rcout<< ll << endl;

  return(ll);
}

void fit_basis(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate, bool vanilla)
{
  double *xx;
  double *omega_i_tmp;
  tree::tree_cp bn;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    if(populate) {
      bn = t.bn(xx,xi);
      node_pointers[i] = bn;
    } else {
      bn = node_pointers[i];
    }

    omega_i_tmp = di.omega + i*di.basis_dim;
    //arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);
    //fv[i] = arma::dot(bn->getm(), omega_i);

    //ok
    if(vanilla) {
      fv[i] = bn->mu(0);
    } else {
      fv[i] = 0.0;
      for(size_t j=0; j<di.basis_dim; ++j) {
        fv[i] += omega_i_tmp[j]*bn->mu(j); //todo: do subset before taking inner product
        //using int var =  (*(node_pointer[i]->getp()).getv()
      }
    }


  }
}

double fit_i_basis(size_t& i, std::vector<tree>& t, xinfo& xi, dinfo& di, bool vanilla)
{
  double *xx;
  double fv = 0.0;
  double *omega_i_tmp;
  tree::tree_cp bn;
  
  for(size_t j=0;j<t.size();j++) {
    xx = di.x + i*di.p;
    bn = t[j].bn(xx,xi); //instead of dropping, using the node_pointer object
    
    omega_i_tmp = di.omega + i*di.basis_dim;
    
    //ok
    if(vanilla) {
      fv += bn->mu(0);
    } else {
      for(size_t k=0; k<di.basis_dim; ++k) {
        fv += omega_i_tmp[k]*bn->mu(k);
      }
    }
    
  }
  
  return fv;
}



//-------------------------------------------------------------------------------
// NEW: tsbart: get sufficients stats for all bottom nodes
//-------------------------------------------------------------------------------

void allsuff_ts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
{
   // Bottom nodes are written to bnv.
   // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
   // Initialize data structures
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y
   double t;          //current t

   bnv.clear();      // Clear the bnv variable if any value is already saved there.
   x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

   typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
                                       // or long int, etc.  We don't have to know that ahead of time by using this notation.
   bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
   sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

   // Resize vectors within sufficient stats to have di.tlen length.
   for(size_t i = 0; i < nb; ++i){
      sv[i].n_vec.zeros(di.tlen);
      sv[i].sy_vec.zeros(di.tlen);
   }

   // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
   //map looks like
   // bottom node 1 ------ 1
   // bottom node 2 ------ 2

   // Sum the y values (sy) and the y^2 values (sy2) for each node and store in sv.
   // Loop through each observation.  Push each obs x down the tree and find its bottom node,
   // then index into the suff stat for the bottom node corresponding to that obs.

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
      y=di.y[i];           // Resolves to r.
      t = di.t[i];         // NEW: Resolves to current t

      tbn = x.bn(xx,xi); // Find bottom node for this observation.
      ni = bnmap[tbn];   // Map bottom node to integer index

      // Update the sufficient stats for the bottom node to which that obs belongs.
      ++(sv[ni].n);
      sv[ni].sy += y;
      sv[ni].sy2 += y*y;

      // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
      uvec id = find(di.tref == t); // Idx of current obs t value.
      sv[ni].n_vec(id) += 1;
      sv[ni].sy_vec(id) += y;

   } // End obs loop.
}

// For het variances.
void allsuffhet_ts(tree& x, xinfo& xi, dinfo& di, double* phi, tree::npv& bnv, std::vector<sinfo>& sv)
{
   // phi are precisions for each observation.

   // Bottom nodes are written to bnv.
   // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
   // Initialize data structures
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y
   double t;          //current t

   bnv.clear();      // Clear the bnv variable if any value is already saved there.
   x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

   typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
   // or long int, etc.  We don't have to know that ahead of time by using this notation.
   bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
   sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

   // Resize vectors within sufficient stats to have di.tlen length.
   for(size_t i = 0; i < nb; ++i){
      sv[i].n0_vec.resize(di.tlen);
      sv[i].n_vec.resize(di.tlen);
      sv[i].sy_vec.resize(di.tlen);

      // Fill with zeros, in case.
      sv[i].n0_vec.fill(0);
      sv[i].n_vec.fill(0);
      sv[i].sy_vec.fill(0);
   }

   // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
   //map looks like
   // bottom node 1 ------ 1
   // bottom node 2 ------ 2

   // Sum the y values (sy) and the y^2 values (sy2) for each node and store in sv.
   // Loop through each observation.  Push each obs x down the tree and find its bottom node,
   // then index into the suff stat for the bottom node corresponding to that obs.

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
      y=di.y[i];           // Resolves to r.
      t = di.t[i];         // NEW: Resolves to current t

      tbn = x.bn(xx,xi); // Find bottom node for this observation.
      ni = bnmap[tbn];   // Map bottom node to integer index

      // Update the sufficient stats for the bottom node to which that obs belongs.
      sv[ni].n0 += 1;
      sv[ni].n += log(phi[i]); //sv[ni].n += phi[i];
      sv[ni].sy += phi[i]*y;
      sv[ni].sy2 += phi[i]*y*y;

      // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
      uvec id = find(di.tref == t); // Idx of current obs t value.
      sv[ni].n0_vec(id) +=1;
      sv[ni].n_vec(id) += phi[i];
      sv[ni].sy_vec(id) += phi[i]*y;

   } // End obs loop.
}

//-------------------------------------------------------------------------------
// NEW: tsbart: get sufficient stats for children of node nx in tree x
// (for birth proposal)
//-------------------------------------------------------------------------------
// Birth proposal, homog variances.
void getsuff_ts(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

         y = di.y[i];   // extract current yi.  resolves to r.
         t = di.t[i];   // extract current ti

         if(xx[v] < xi[v][c]) { // Update left.
            sl.n++;
            sl.sy += y;
            sl.sy2 += y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sl.n_vec(id) += 1;
            sl.sy_vec(id) += y;

         } else { //Update right.
            sr.n++;
            sr.sy += y;
            sr.sy2 += y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sr.n_vec(id) += 1;
            sr.sy_vec(id) += y;

         }
      }
   }
}

// Birth proposal, het variances.
void getsuffhet_ts(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n0=0;sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n0=0;sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n0_vec = zeros(tlen); sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n0_vec = zeros(tlen); sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

         y = di.y[i];   // extract current yi.
         t = di.t[i];   // extract current ti

         if(xx[v] < xi[v][c]) { // Update left.
            sl.n0 += 1;
            sl.n += log(phi[i]); //+= phi[i];
            sl.sy += phi[i]*y;
            sl.sy2 += phi[i]*y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sl.n0_vec(id) += 1;
            sl.n_vec(id) += phi[i];
            sl.sy_vec(id) += phi[i]*y;

         } else { //Update right.
            sr.n0 += 1;
            sr.n += log(phi[i]); //phi[i];
            sr.sy += phi[i]*y;
            sr.sy2 += phi[i]*y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sr.n0_vec(id) += 1;
            sr.n_vec(id) += phi[i];
            sr.sy_vec(id) += phi[i]*y;
         }
      }
   }
}

//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
// (for death proposal)
//--------------------------------------------------
// Death proposal, homog variance.
void getsuff_ts(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);

      y = di.y[i];   // extract current yi.
      t = di.t[i];   // extract current ti

      uvec id = find(di.tref == t); // Idx of current obs t value.

      if(bn==nl) {
         sl.n++;
         sl.sy += y;
         sl.sy2 += y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         sl.n_vec(id) += 1;
         sl.sy_vec(id) += y;
      }

      if(bn==nr) {
         sr.n++;
         sr.sy += y;
         sr.sy2 += y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         sr.n_vec(id) += 1;
         sr.sy_vec(id) += y;
      }
   }
}

// For death proposal, het variances.
void getsuffhet_ts(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n0=0;sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n0=0;sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n0_vec = zeros(tlen); sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n0_vec = zeros(tlen); sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);

      y = di.y[i];   // extract current yi.
      t = di.t[i];   // extract current ti

      if(bn==nl) {
         y = di.y[i];
         sl.n0 += 1;
         sl.n += log(phi[i]); //phi[i];
         sl.sy += phi[i]*y;
         sl.sy2 += phi[i]*y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         uvec id = find(di.tref == t); // Idx of current obs t value.
         sl.n0_vec(id) += 1;
         sl.n_vec(id) += phi[i];
         sl.sy_vec(id) += phi[i]*y;
      }

      if(bn==nr) {
         y = di.y[i];
         sr.n0 += 1;
         sr.n += log(phi[i]); //phi[i];
         sr.sy += phi[i]*y;
         sr.sy2 += phi[i]*y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         uvec id = find(di.tref == t); // Idx of current obs t value.
         sr.n0_vec(id) += 1;
         sr.n_vec(id) += phi[i];
         sr.sy_vec(id) += phi[i]*y;
      }
   }
}

//--------------------------------------------------
// draw all the bottom node mu's

// For homog variances.
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{
   tree::npv bnv;
   std::vector<sinfo> sv(di.tlen);
   allsuff_ts(t,xi,di,bnv,sv);

   List post_pars;
   mat Phi;
   vec m;
   vec mu_draw;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      // Draw new mu value from MVN.
      post_pars = mvn_post_util(pi.sigma, pi.mu0, pi.Prec0, sv[i].n_vec, sv[i].sy_vec);

      Phi = Rcpp::as<arma::mat>(post_pars["Phi"]);
      m = Rcpp::as<arma::vec>(post_pars["m"]);
      mu_draw = rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(mu_draw);

      // Check for NA result.
      if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
          Rcpp::stop("drmu failed");
      }
   }
}

// For heterogeneous variances.
void drmuhet(tree& t, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen)
{
   tree::npv bnv;
   std::vector<sinfo> sv(di.tlen);
   allsuffhet_ts(t,xi,di,phi,bnv,sv);

   List post_pars;
   mat Phi;
   vec m;
   vec mu_draw;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      // Draw new mu value from MVN.
      post_pars = mvn_post_util_het(pi.mu0, pi.Prec0, sv[i].n0_vec, sv[i].n_vec, sv[i].sy_vec);

      Phi = Rcpp::as<arma::mat>(post_pars["Phi"]);
      m = Rcpp::as<arma::vec>(post_pars["m"]);
      mu_draw = rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(mu_draw);

      // Check for NA result.
      if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
         for(size_t i=0; i<di.n; ++i) Rcout << *(di.x + i*di.p) <<" "; //*(x + p*i+j)
         Rcpp::stop("drmu failed");
      }
   }
}

//--------------------------------------------------
// normal density N(x, mean, variance)
double pn(double x, double m, double v)
{
	double dif = x-m;
	return exp(-.5*dif*dif/v)/sqrt(2*PI*v);
}

//--------------------------------------------------
// draw from discrete distributin given by p, return index
int rdisc(double *p, RNG& gen)
{

	double sum;
	double u = gen.uniform();

    int i=0;
    sum=p[0];
    while(sum<u) {
		i += 1;
		sum += p[i];
    }
    return i;
}

//--------------------------------------------------
//evalute tree tr on grid given by xi and write to os
void grm(tree& tr, xinfo& xi, std::ostream& os)
{
	size_t p = xi.size();
	if(p!=2) {
		cout << "error in grm, p !=2\n";
		return;
	}
	size_t n1 = xi[0].size();
	size_t n2 = xi[1].size();
	tree::tree_cp bp; //pointer to bottom node
	double *x = new double[2];
	for(size_t i=0;i!=n1;i++) {
		for(size_t j=0;j!=n2;j++) {
			x[0] = xi[0][i];
			x[1] = xi[1][j];
			bp = tr.bn(x,xi);
			os << x[0] << " " << x[1] << " " << bp->getm() << " " << bp->nid() << endl;
		}
	}
	delete[] x;
}

//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
	int L,U;
	bool v_found = false; //have you found a variable you can split on
	size_t v=0;
	while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
		L=0; U = xi[v].size()-1;
		n->rg(v,&L,&U);
		if(U>=L) v_found=true;
		v++;
	}
	return v_found;
}

//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots)
{
	double pb;  //prob of birth to be returned
	tree::npv bnv; //all the bottom nodes
	t.getbots(bnv);
	for(size_t i=0;i!=bnv.size();i++)
		if(cansplit(bnv[i],xi)) goodbots.push_back(bnv[i]);
	if(goodbots.size()==0) { //are there any bottom nodes you can split on?
		pb=0.0;
	} else {
		if(t.treesize()==1) pb=1.0; //is there just one node?
		else pb=pi.pb;
	}
	return pb;
}

//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
	int L,U;
	for(size_t v=0;v!=xi.size();v++) {//try each variable
		L=0; U = xi[v].size()-1;
		n->rg(v,&L,&U);
		if(U>=L) goodvars.push_back(v);
	}
}

//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi)
{
	if(cansplit(n,xi)) {
		return pi.alpha/pow(1.0+n->depth(),pi.beta);
	} else {
		return 0.0;
	}
}

//--------------------------------------------------

//get counts for all bottom nodes
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
	size_t ni;         //the  index into vector of the current bottom node
	double *xx;        //current x
	double y;          //current y

	bnv.clear();
	x.getbots(bnv);

	typedef tree::npv::size_type bvsz;
//	bvsz nb = bnv.size();

  std::vector<int> cts(bnv.size(), 0);

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;

	for(size_t i=0;i<di.n;i++) {
		xx = di.x + i*di.p;
		y=di.y[i];

		tbn = x.bn(xx,xi);
		ni = bnmap[tbn];

    cts[ni] += 1;
	}
  return(cts);
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   tree::npv& bnv, //vector of pointers to bottom nodes
                   int sign)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
	double *xx;        //current x
	double y;          //current y

	typedef tree::npv::size_type bvsz;
//	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index

	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   std::map<tree::tree_cp,size_t>& bnmap,
                   int sign)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
	double y;          //current y
  /*
	typedef tree::npv::size_type bvsz;
	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index
	*/
	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}


void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   std::map<tree::tree_cp,size_t>& bnmap,
                   int sign,
                   tree::tree_cp &tbn
                   )
{
  //tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  /*
	typedef tree::npv::size_type bvsz;
	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index
	*/
	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}

bool min_leaf(int minct, std::vector<tree>& t, xinfo& xi, dinfo& di) {
  bool good = true;
  tree::npv bnv;
  std::vector<int> cts;
  int m = 0;
  for (size_t tt=0; tt<t.size(); ++tt) {
    cts = counts(t[tt], xi, di, bnv);
    m = std::min(m, *std::min_element(cts.begin(), cts.end()));
    if(m<minct) {
      good = false;
      break;
    }
  }
  return good;
}

//--------------------------------------------------
//fit for multiple data points, not by reference.
void fit(tree& t, xinfo& xi, dinfo& di, vec& fv)
{
   double *xx;
   tree::tree_cp bn;
   fv.resize(di.n);

   arma::uvec id; //idx of current obs t value.

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);

      // Find index of mu (time point) corresponding to each obs.  Use this mu.
      fv[i] = bn->getm(di.t[i], di.tref);
   }
}

//--------------------------------------------------
//fit for multiple data points, by reference.
void fit(tree& t, xinfo& xi, dinfo& di, double* fv)
{
   double *xx;
   tree::tree_cp bn;

   arma::uvec id; //idx of current obs t value.


   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);

      // Find index of mu (time point) corresponding to each obs.  Use this mu.
      fv[i] = bn->getm(di.t[i], di.tref);
   }
}

mat coef_basis(tree& t, xinfo& xi, dinfo& di)
{
  double *xx;

  tree::tree_cp bn;

  mat out(di.basis_dim,di.n);

  for(size_t i=0; i<di.n; i++) {
    xx = di.x + i*di.p;
    bn = t.bn(xx,xi);
    out.col(i) = bn->getm();
    /*
    for(size_t j=0; j<di.basis_dim; ++j) {
      out.at(i,j) = bn->mu.at(j);
    }
     */
  }
  return(out);
}

//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv)
{
	double *xx;
	tree::tree_cp bn;
	pv.resize(di.n);
	for(size_t i=0;i<di.n;i++) {
		xx = di.x + i*di.p;
		bn = t.bn(xx,xi);
		pv[i] = bn->nid();
	}
}

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
	cout << "xinfo: \n";
	for(size_t v=0;v!=xi.size();v++) {
		cout << "v: " << v << endl;
		for(size_t j=0;j!=xi[v].size();j++) cout << "j,xi[v][j]: " << j << ", " << xi[v][j] << endl;
	}
	cout << "\n\n";
}

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
	double xinc;

	//compute min and max for each x
	std::vector<double> minx(p,INFINITY);
	std::vector<double> maxx(p,-INFINITY);
	double xx;
	for(size_t i=0;i<p;i++) {
		for(size_t j=0;j<n;j++) {
			xx = *(x+p*j+i);
			if(xx < minx[i]) minx[i]=xx;
			if(xx > maxx[i]) maxx[i]=xx;
		}
	}
	//make grid of nc cutpoints between min and max for each x.
	xi.resize(p);
	for(size_t i=0;i<p;i++) {
		xinc = (maxx[i]-minx[i])/(nc+1.0);
		xi[i].resize(nc);
		for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
	}
}
// get min/max needed to make cutpoints
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx)
{
	double xx;

	for(size_t i=0;i<p;i++) {
		for(size_t j=0;j<n;j++) {
			xx = *(x+p*j+i);
			if(xx < minx[i]) minx[i]=xx;
			if(xx > maxx[i]) maxx[i]=xx;
		}
	}
}
//make xinfo = cutpoints give the minx and maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx)
{
	double xinc;
	//make grid of nc cutpoints between min and max for each x.
	xi.resize(p);
	for(size_t i=0;i<p;i++) {
		xinc = (maxx[i]-minx[i])/(nc+1.0);
		xi[i].resize(nc);
		for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
	}
}

// Check if a vector is sorted.  For checking z and zpred for causal funbart.
bool is_sort(arma::vec x) {
     int n=x.n_elem;
     for (int i=0; i<n-1; ++i)
         if (x[i] < x[i+1]) return false;
     return true;
}

//--------------------------------------------------
//log of the integrated likelihood
double lil(double n, double sy, double sy2, double sigma, double tau)
{
   double yb,yb2,S,sig2,d;
   double sum, rv;

   yb = sy/n;
   yb2 = yb*yb;
   S = sy2 - (n*yb2);
   sig2 = sigma*sigma;
   d = n*tau*tau + sig2;
   sum = S/sig2 + (n*yb2)/d;
   rv = -(n*LTPI/2.0) - (n-1)*log(sigma) -log(d)/2.0;
   rv = rv -sum/2.0;
   return rv;
}

double lilhet(double n, double sy, double sy2, double sigma, double tau)
{
   double d = 1/(tau*tau) + n;// n is \sum phi_i for het

   double out = -log(tau) - 0.5*log(d);
   out += 0.5*sy*sy/d - 0.5*sy2;
   return out;
}


//get sufficients stats for all bottom nodes (sy, sy2)
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
{
   // Bottom nodes are written to bnv.
   // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
   // Initialize data structures
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y

   bnv.clear();      // Clear the bnv variable if any value is already saved there.
   x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

   // Not sure what this part here is doing.
   typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
   // or long int, etc.  We don't have to know that ahead of time by using this notation.
   bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
   sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

   // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
   //map looks like
   // bottom node 1 ------ 1
   // bottom node 2 ------ 2

   // Sum the y values (sy) and the y^2 values (sy2) for each node and store in sv.
   // Loop through each observation.  Push each obs x down the tree and find its bottom node,
   // then index into the suff stat for the bottom node corresponding to that obs.
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
      y=di.y[i];           // Resolves to r.

      tbn = x.bn(xx,xi); // Find bottom node for this observation.
      ni = bnmap[tbn];   // Map bottom node to integer index

      // Update the sufficient stats for the bottom node to which that obs belongs.
      ++(sv[ni].n);
      sv[ni].sy += y;
      sv[ni].sy2 += y*y;
   }
}


//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         y = di.y[i];
         if(xx[v] < xi[v][c]) {
            sl.n++;
            sl.sy += y;
            sl.sy2 += y*y;
         } else {
            sr.n++;
            sr.sy += y;
            sr.sy2 += y*y;
         }
      }
   }
}

//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      if(bn==nl) {
         y = di.y[i];
         sl.n++;
         sl.sy += y;
         sl.sy2 += y*y;
      }
      if(bn==nr) {
         y = di.y[i];
         sr.n++;
         sr.sy += y;
         sr.sy2 += y*y;
      }
   }
}

void fit_linear(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate)
{
  double *xx;
  double *omega_i_tmp;
  tree::tree_cp bn;
  size_t var;

  //Getting the node pointers for every observation
  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    if(populate) {
      bn = t.bn(xx,xi);
      node_pointers[i] = bn;
    } else {
      bn = node_pointers[i];
    }

    omega_i_tmp = di.omega + i*(di.p+1);


    if(bn->nid()==1){// If is root node
      //Then the the rule is on the same node
      var = bn->getv()+1;
    }else{
      //Otherwise, is on the parent
      var =  (*(bn->getp())).getv() + 1;
    }

    //The fit only happens in the intercept and in the slope
    fv[i] = 0.0;
    fv[i] = omega_i_tmp[0]*bn->mu(0) + omega_i_tmp[var]*bn->mu(1);
    //Rcpp::Rcout << bn->mu(0)  << " mu0 -" << bn->mu(1) << "mu1" << endl;

  }
}

double fit_i_linear(size_t& i, std::vector<tree>& t, xinfo& xi, dinfo& di)
{
  double *xx;
  double fv = 0.0;
  double *omega_i_tmp;
  tree::tree_cp bn;
  size_t var;

  //Getting the node pointers for every observation
  for(size_t j=0; j<t.size(); ++j){
    xx = di.x + i*di.p;

    bn = t[j].bn(xx,xi);

    omega_i_tmp = di.omega + i*(di.p+1);

    if(bn->nid()==1){// If is root node
      //Then the the rule is on the same node
      var = bn->getv()+1;
    }else{
      //Otherwise, is on the parent
      var =  (*(bn->getp())).getv() + 1;
    }

    //The fit only happens in the intercept and in the slope
    fv += omega_i_tmp[0]*bn->mu(0) + omega_i_tmp[var]*bn->mu(1);
  }

  return fv;
}


void drmu_linear(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{

  bool debug = false;
  tree::npv bnv;
  std::vector<sinfo> sv; //will be resized in allsuff
  tree::npv bnv0; std::vector<sinfo> sv0;
  //debug is broken, would need to eat node_pointers. suff stat code works anyhow
  //  if(debug) allsuff_basis(t,xi,di,bnv0,sv0,node_pointers);

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  t.getbots(bnv);   // Save bottom nodes for x to bnv variable.
  tree::npv::size_type nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.
  arma::vec v = zeros(nb);
  
  double prior_prec = pi.Prec0(0,0);
  double prior_mean = pi.mu0(0);

  //Get the suff stats on every node
  for(tree::npv::size_type i=0; i<nb; i++) {
    sv[i] = bnv[i]->s;
    if(bnv[i]->nid()==1){
      v[i] = bnv[i]->getv();
    }else{
      v[i] = (*(bnv[i]->getp())).getv();
    }

    if(debug) {
      if(abs(sv[i].sy - sv0[i].sy)>1e-8) {

        Rcout << "depth "<< bnv[i]->depth() << endl;
        Rcout << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
        Rcout << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
        stop("shit");
      }
      else {
        //Rcout << "depth "<< bnv[i]->depth() << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec <<" new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << endl;
        //Rcoutt << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
        //Rcoutt << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
      }
      sv[i] = sv0[i];
    }

  }
  

  mat Phi;
  vec m;
  vec beta_draw;
  mat Prec_mi=eye(2,2); //Will be updated with m_i

  double s2 = (pi.sigma*pi.sigma);
  for(tree::npv::size_type i=0;i!=bnv.size();i++) {

    if(pi.gprior==1){
      Prec_mi(0,0) = sv[i].WtW(0,0)*2*pi.ntree/(sv[i].n);
      Prec_mi(1,0) = sv[i].WtW(1,0)*2*pi.ntree/(sv[i].n);
      Prec_mi(0,1) = sv[i].WtW(0,1)*2*pi.ntree/(sv[i].n);
      Prec_mi(1,1) = sv[i].WtW(1,1)*2*pi.ntree/(sv[i].n);
      // Prec_mi(0,0) = sv[i].WtW(0,0)*2*pi.ntree/(sv[i].n*pi.Prec0(1,1));
      // Prec_mi(1,0) = sv[i].WtW(1,0)*2*pi.ntree/(sv[i].n*pi.Prec0(1,1));
      // Prec_mi(0,1) = sv[i].WtW(0,1)*2*pi.ntree/(sv[i].n*pi.Prec0(1,1));
      // Prec_mi(1,1) = sv[i].WtW(1,1)*2*pi.ntree/(sv[i].n*pi.Prec0(1,1));
      // // Prec_mi(0,0) = sv[i].WtW(0,0)*2*pi.ntree/(sv[i].n*s2*pi.Prec0(1,1));
      // Prec_mi(1,0) = sv[i].WtW(1,0)*2*pi.ntree/(sv[i].n*s2*pi.Prec0(1,1));
      // Prec_mi(0,1) = sv[i].WtW(0,1)*2*pi.ntree/(sv[i].n*s2*pi.Prec0(1,1));
      // Prec_mi(1,1) = sv[i].WtW(1,1)*2*pi.ntree/(sv[i].n*s2*pi.Prec0(1,1));
    }else{
      //This is not correct, does not make sense, but it is not used anyway
      Prec_mi(1,1)=Prec_mi(1,1)*pi.Prec0(1,1);
      Prec_mi(0,0)=Prec_mi(0,0)*pi.Prec0(0,0);
    }

    //Prec_mi(1,1)=Prec_mi(1,1)*pi.Prec0(1,1);
    //Prec_mi(0,0)=Prec_mi(0,0)*pi.Prec0(0,0); //Receive the prior

    //Then update the variance prior using m_i
    //Prec_mi(1,1)=Prec_mi(1,1)*(pi.mi[(v[i])]/pi.ntree);
    //Removing the mi

    vec tt = Prec_mi*pi.mu0;

    Phi = sv[i].WtW/s2 + Prec_mi;
    m = tt + sv[i].sy_vec/s2;

    //Rcpp::Rcout<<"before draw" <<endl;
    //Rcpp::Rcout<<"WtW"<<sv[i].WtW <<endl;
    //Rcpp::Rcout<<"Phi"<<Phi <<endl;
    //Rcpp::Rcout<<"g"<< pi.ntree/(sv[i].n*s2) <<endl;
    //Rcpp::Rcout<<"n"<< sv[i].n <<endl;
    //Rcpp::Rcout<<"s2"<< s2 <<endl;
    //Rcpp::Rcout<<"Wty"<< sv[i].sy_vec <<endl;
    
    // Rcpp::Rcout<< "Beta_draw" << sv[i].WtW <<endl;
    // Rcpp::Rcout<< "Determinant " << arma::det(sv[i].WtW) << endl;
    // Rcpp::Rcout<< "sv[i].WtW " << sv[i].WtW << endl;
    
    
    
    //Solve the case for categorical binary variables (multicollinearity)
    // Rcpp::Rcout<< "sv[i].WtW " << sv[i].WtW << endl;
    // Rcpp::Rcout<< "arma::det(sv[i].WtW) " << arma::det(sv[i].WtW) << endl;

    //Since in these cases sum of x and sum of x square are the same
    //and that the positions (0,1) and (1,0) are the same
    //then I only need to check those two positions to know if it is fine
    //This will only happen when all x's are either 1 or 0,
    //therefore, it captures the binary case
    //if((sv[i].WtW(0,1)-sv[i].WtW(1,1)) < 1e-3){
    if((abs(sv[i].WtW(1,1)-sv[i].WtW(0,1)*sv[i].WtW(1,0)/sv[i].WtW(0,0))/sv[i].WtW(0,0) < 1e-8)){
      double tt_a = prior_prec*prior_mean;
      //Rcpp::Rcout<< "prior_prec 1: " << prior_prec << endl;
      double Phi_a = sv[i].WtW(0,0)/s2 + prior_prec;
      //Rcpp::Rcout<< "sv[i].n0 3: " << sv[i].WtW(0,0) << endl;
      double m_a = tt_a + sv[i].sy_vec(0)/s2;
      //Rcpp::Rcout<< "sy_vec 1: " << sv[i].sy_vec << endl;
      beta_draw.set_size(2);
      // Rcpp::Rcout<< "Is the solve here? 1" << endl;
      
      beta_draw(0) = m_a/Phi_a + gen.normal(0,1)/sqrt(Phi_a);
      
      // Rcpp::Rcout<< "prior_prec " << prior_prec << endl;
      // Rcpp::Rcout<< "prior_mean " << prior_mean << endl;
      // Rcpp::Rcout<< "mu_a " << m_a/Phi_a << endl;
      // Rcpp::Rcout<< "Var_a " << (1/Phi_a) << endl;
      // Rcpp::Rcout<< "s2 " << s2 << endl;
      // Rcpp::Rcout<< "sv[i].WtW(0,0) " << sv[i].WtW << endl;
      // Rcpp::Rcout<< "n " << sv[i].n << endl;
      // Rcpp::Rcout<< "n0 " << sv[i].n0 << endl;
      // Rcpp::Rcout<< "sv[i].sy_vec(0) " << sv[i].sy_vec << endl;
      // Rcpp::Rcout<< "sv[i].sy " << sv[i].sy << endl;

      
      beta_draw(1) = 0;
      // Rcpp::Rcout<< "not here" << endl;
      //Rcpp::Rcout<< "beta_draw 4: " << beta_draw << endl;
    }else{
      // Rcpp::Rcout<< "here 2?" << endl;
      // Rcpp::Rcout<< "mu " << m << endl;
      // Rcpp::Rcout<< "Var " << Phi << endl;
      beta_draw = rmvnorm_post(m, Phi);
      // Rcpp::Rcout<< "not here" << endl;
    }
    
    //Rcpp::Rcout<<"After draw" <<endl;

    // Assign bottom node values to new mu draw.
    bnv[i] -> setm(beta_draw);

    // Check for NA result.
    if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
      Rcpp::stop("drmu failed");
    }
  }

}


void allsuff_linear_death(tree& x,
                          tree::tree_cp nx,
                          xinfo& xi,
                          dinfo& di,
                          tree::npv& bnv,
                          std::vector<sinfo>& sv,
                          sinfo& st,
                          std::vector<tree::tree_cp>& node_pointers,
                          size_t newvar)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].n = 0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  double *omega_i_tmp;
  arma::vec omega_i(di.basis_dim);
  arma::vec omega_i_death(di.basis_dim);
  bool in_candidate_nog, prune;
  size_t split_var;
  size_t death_var;


  //Getting the variable for new suff stats after prune
  if(nx->nid()==1){//If root node, then
    death_var = newvar+1; //Variable will be stored on the root node
  }else{
    //Otherwise, is on the parent
    death_var =  (*(nx->getp())).getv() + 1;
  }


  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    y = di.y[i];

    // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    //Find the split rule for each node
    if(tbn->nid()==1){
      split_var = tbn->getv()+1;
    }else{
      split_var =  (*(tbn->getp())).getv() + 1;
    }

    //Will decide if the variable will be in the new leaf
    prune = false;

    tree::tree_cp lptr = nx->getl();
    tree::tree_cp rptr = nx->getr();

    in_candidate_nog = ((tbn==lptr) | (tbn==rptr));

    //For the proposed split we know the variable is v so when we compute
    //summary stats for sl ad sr, use variable v

    prune  = in_candidate_nog;

    omega_i_tmp = di.omega + i*(di.p+1);


    //get design vector
    omega_i[0] = omega_i_tmp[0];
    omega_i[1] = omega_i_tmp[split_var];

    omega_i_death[0] = omega_i_tmp[0];
    omega_i_death[1] = omega_i_tmp[death_var];


    //////////////////////////
    double b = 0;
    
    ++(sv[ni].n0);
      
    //If it is in the leaf, then sum
    sv[ni].n += omega_i_death[0];
    if(prune) st.n += omega_i_death[0];

    sv[ni].sy_vec += y*omega_i;
    if(prune) st.sy_vec += y*omega_i_death;
    
    sv[ni].sy += y;
    if(prune) st.sy += y;
    

    for(size_t j=0; j<di.basis_dim; ++j) {
      double a = omega_i[j]*omega_i[j];
      if(prune) b = omega_i_death[j]*omega_i_death[j];
      sv[ni].WtW.at(j,j) += a;
      if(prune) st.WtW.at(j,j) += b;
      for(size_t g=0; g<j; ++g) {
        double a = omega_i[j]*omega_i[g];
        if(prune) b = omega_i_death[j]*omega_i_death[g];
        sv[ni].WtW.at(g,j) += a;
        if(prune) st.WtW.at(g,j) += b;
      }
    }
  }

  //Fill the other half of the matrix
  for(size_t q=0; q<sv.size(); ++q) {
    for(size_t j=0; j<di.basis_dim; ++j) {
      for(size_t g=0; g<j; ++g) {
        sv[q].WtW(j,g) = sv[q].WtW(g,j);
        st.WtW(j,g) = st.WtW(g,j);
      }
    }
  }
}


void allsuff_linear_birth(tree& x, tree::tree_cp nx,
                          size_t v,
                          size_t c,
                          xinfo& xi,
                          dinfo& di,
                          tree::npv& bnv,
                          std::vector<sinfo>& sv,
                          sinfo& sl,
                          sinfo& sr,
                          std::vector<tree::tree_cp>& node_pointers)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].n = 0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  double *omega_i_tmp;
  arma::vec omega_i(di.basis_dim);
  arma::vec omega_i_birth(di.basis_dim);

  size_t birth_var = v+1;

  bool in_candidate_nog, left, right;
  size_t split_var;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    y = di.y[i];
    
    // Rcpp::Rcout<< "y[i] " << di.y[i] << endl;

    // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    //FInd the rule for each observation leaf node
    if(tbn->nid()==1){
      split_var = tbn->getv()+1;
    }else{
      split_var =  (*(tbn->getp())).getv() + 1;
    }

    //Get if the variable is in one of the new nodes
    left = false; right = false;

    in_candidate_nog = (tbn == nx);

    //For the proposed split we know the variable is v so when we compute
    //summary stats for sl ad sr, use variable v

    left  = in_candidate_nog & (xx[v] < xi[v][c]);
    right = in_candidate_nog & !(xx[v] < xi[v][c]);
    
    
    
    //Rcpp::Rcout<< "Posicao: " << xx[v] << endl;

    ++(sv[ni].n0);
    if(left) sl.n0 += 1;
    if(right) sr.n0 += 1;

    omega_i_tmp = di.omega + i*(di.p+1);
    //get design vector
    omega_i[0] = omega_i_tmp[0];
    omega_i[1] = omega_i_tmp[split_var];

    omega_i_birth[0] = omega_i_tmp[0];
    omega_i_birth[1] = omega_i_tmp[birth_var];
    
    sv[ni].n += omega_i_birth[0];
    if(left) sl.n += omega_i_birth[0];
    if(right) sr.n += omega_i_birth[0];
    
    //if(left) Rcpp::Rcout<< "Omega left: " << omega_i_birth << endl;
    //if(right) Rcpp::Rcout<< "Omega right: " << omega_i_birth << endl;

    //Computing the suff stats

    sv[ni].sy_vec += y*omega_i;
    if(left) sl.sy_vec += y*omega_i_birth;
    if(right) sr.sy_vec += y*omega_i_birth;
    
    // Rcpp::Rcout<< "sy_vec[i] " << sv[ni].sy_vec << endl;
    
    sv[ni].sy += y;
    if(left) sl.sy += y;
    if(right) sr.sy += y;

    // Rcpp::Rcout<< "sy[i] " << sv[ni].sy << endl;
    
    for(size_t j=0; j<di.basis_dim; ++j) {
      double a = omega_i[j]*omega_i[j];
      double b = omega_i_birth[j]*omega_i_birth[j];
      sv[ni].WtW.at(j,j) += a;
      if(left) sl.WtW.at(j,j) += b;
      if(right) sr.WtW.at(j,j) += b;
      for(size_t g=0; g<j; ++g) {
        double a = omega_i[j]*omega_i[g];
        double b = omega_i_birth[j]*omega_i_birth[g];
        sv[ni].WtW.at(g,j) += a;
        if(left) sl.WtW.at(g,j) += b;
        if(right) sr.WtW.at(g,j) += b;
      }
    }
  }

  //Fill the other half of the matrix
  for(size_t q=0; q<sv.size(); ++q) {
    for(size_t j=0; j<di.basis_dim; ++j) {
      for(size_t g=0; g<j; ++g) {
        sv[q].WtW(j,g) = sv[q].WtW(g,j);
        sl.WtW(j,g) = sl.WtW(g,j);
        sr.WtW(j,g) = sr.WtW(g,j);
      }
    }
  }

}

double lil_linear(sinfo& s, pinfo& pi, size_t v){

  double ll=0;
  double tt = 0;

  double s2 = (pi.sigma)*(pi.sigma);

  mat Prec_mi=eye(2,2);

  //Update prior with m_i
  if(pi.gprior==1){
    Prec_mi(0,0) = s.WtW(0,0)*2*pi.ntree/(s.n);
    Prec_mi(1,0) = s.WtW(1,0)*2*pi.ntree/(s.n);
    Prec_mi(0,1) = s.WtW(0,1)*2*pi.ntree/(s.n);
    Prec_mi(1,1) = s.WtW(1,1)*2*pi.ntree/(s.n);
    // Prec_mi(0,0) = s.WtW(0,0)*2*pi.ntree/(s.n*pi.Prec0(1,1));
    // Prec_mi(1,0) = s.WtW(1,0)*2*pi.ntree/(s.n*pi.Prec0(1,1));
    // Prec_mi(0,1) = s.WtW(0,1)*2*pi.ntree/(s.n*pi.Prec0(1,1));
    // Prec_mi(1,1) = s.WtW(1,1)*2*pi.ntree/(s.n*pi.Prec0(1,1));
    // Prec_mi(0,0) = s.WtW(0,0)*2*pi.ntree/(s.n*s2*pi.Prec0(1,1));
    // Prec_mi(1,0) = s.WtW(1,0)*2*pi.ntree/(s.n*s2*pi.Prec0(1,1));
    // Prec_mi(0,1) = s.WtW(0,1)*2*pi.ntree/(s.n*s2*pi.Prec0(1,1));
    // Prec_mi(1,1) = s.WtW(1,1)*2*pi.ntree/(s.n*s2*pi.Prec0(1,1));
  }else{
    Prec_mi(1,1)=Prec_mi(1,1)*pi.Prec0(1,1);
    Prec_mi(0,0)=Prec_mi(0,0)*pi.Prec0(0,0);
  }

  //Prec_mi(1,1)=Prec_mi(1,1)*pi.Prec0(1,1);
  //Prec_mi(0,0)=Prec_mi(0,0)*pi.Prec0(0,0);

  //Prec_mi(1,1)=Prec_mi(1,1)*(pi.mi[(v)]/pi.ntree);
  //removing the mi
  
  // Rcpp::Rcout<< "s.WtW 1: " << s.WtW(0,0) << endl;
  // Rcpp::Rcout<< "s.WtW 2: " << s.WtW(0,1) << endl;
  // Rcpp::Rcout<< "s.WtW 3: " << s.WtW(1,0) << endl;
  // Rcpp::Rcout<< "s.WtW 4: " << s.WtW(1,1) << endl;
  // Rcpp::Rcout<< "Determinant: " << (abs(s.WtW(1,1)-s.WtW(0,1)*s.WtW(1,0))/s.WtW(0,0)) << endl;
  // Rcpp::Rcout<< "Determinant: " << (abs(s.WtW(1,1)-s.WtW(0,1)*s.WtW(1,0))) << endl;
  // //if((s.WtW(0,1)-s.WtW(1,1)) < 1e-3){
  
  
  
  if((abs(s.WtW(1,1)-s.WtW(0,1)*s.WtW(1,0)/s.WtW(0,0))/s.WtW(0,0) < 1e-8)){
    // if(rcond(s.WtW) < 1e-17){
    // Rcpp::Rcout<< "Inside" << endl;
    // Rcpp::Rcout<< "s.WtW 1: " << s.WtW(0,0) << endl;
    // Rcpp::Rcout<< "s.WtW 2: " << s.WtW(0,1) << endl;
    // Rcpp::Rcout<< "s.WtW 3: " << s.WtW(1,0) << endl;
    // Rcpp::Rcout<< "s.WtW 4: " << s.WtW(1,1) << endl;


    // }
    // double tt_a = prior_prec*prior_mean;
    // double Phi_a = sv[i].WtW(0,0)/s2 + prior_prec;
    // double m_a = tt_a + sv[i].sy_vec(0)/s2;
    // beta_draw.set_size(2);
    // beta_draw(0) = m_a/Phi_a + gen.normal(0,1)/sqrt(Phi_a);
    // 
    // 
    // mat precpred = pi.Prec0(0,0) + s.WtW/s2;
    // 
    // double dotp = dot(s.sy_vec/s2, solve(precpred, s.sy_vec/s2));
    // double tt = 0;
    // 
    // tt = -0.5*log(det(precpred)) - 0.5*pi.logdetSigma0;
    // 
    // ll = -0.5*s.n*log(s2) + tt + 0.5*dotp;
    
    
    //Rcpp::Rcout<< "s.WtW 1: " << s.WtW(0,0) << endl;
    
    double precpred = pi.Prec0(0,0) + s.WtW(0,0)/s2;
    
    //double dotp = dot(s.sy_vec/s2, solve(precpred, s.sy_vec/s2));
    double dotp = (s.sy_vec(0)/s2)*(s.sy_vec(0)/s2)/precpred;
    double tt = 0;
    
    tt = -0.5*log(precpred) - 0.5*pi.logdetSigma0;
    
    ll = -0.5*s.n*log(s2) + tt + 0.5*dotp;
    
    // Rcpp::Rcout<< "solve " << precpred*s.sy_vec(0)/s2 << endl;
    // Rcpp::Rcout<< "s.sy_vec " << s.sy_vec(0) << endl;
    // Rcpp::Rcout<< "s.WtW " << s.WtW(0,0) << endl;
    // Rcpp::Rcout<< "pi.Prec0 " << Prec_mi(0,0) << endl;
    // Rcpp::Rcout<< "det(precpred) " << det(precpred) << endl;
    // Rcpp::Rcout<< "logdetSigma0 " << pi.logdetSigma0 << endl;
    // Rcpp::Rcout<< "s2 " << s2 << endl;
    // Rcpp::Rcout<< "dotp " << dotp << endl;
    // Rcpp::Rcout<< "ll " << ll << endl;
    
    
    return(ll);
  }else{
  // Rcpp::Rcout<< "No zero inside, matrix is good " << endl;
  // // // if(rcond(s.WtW) < 1e-15){
  // Rcpp::Rcout<< "Outside" << endl;
  // Rcpp::Rcout<< "s.WtW 1: " << s.WtW(0,0) << endl;
  // Rcpp::Rcout<< "s.WtW 2: " << s.WtW(0,1) << endl;
  // Rcpp::Rcout<< "s.WtW 3: " << s.WtW(1,0) << endl;
  // Rcpp::Rcout<< "s.WtW 4: " << s.WtW(1,1) << endl;


  // }
  
    
  mat precpred = Prec_mi + s.WtW/s2;

  // Rcpp::Rcout<< "precpred: " << precpred << endl;
  // Rcpp::Rcout<< "sy_vec: " << s.sy_vec << endl;
  // Rcpp::Rcout<< "s2: " << s2 << endl;
  // Rcpp::Rcout<< "s.WtW 1: " << s.WtW(0,0) << endl;
  // Rcpp::Rcout<< "s.WtW 2: " << s.WtW(0,1) << endl;
  // Rcpp::Rcout<< "s.WtW 3: " << s.WtW(1,0) << endl;
  // Rcpp::Rcout<< "s.WtW 4: " << s.WtW(1,1) << endl;
  // 
  // Rcpp::Rcout<<"before solve" <<endl;
  // Rcpp::Rcout<< "here 3?" << endl;
  double dotp = dot(s.sy_vec/s2, solve(precpred, s.sy_vec/s2));
  // Rcpp::Rcout<< "not here" << endl;
  // Rcpp::Rcout<<"after solve" <<endl;



  double val1;
  double sign1;
  log_det(val1, sign1, precpred);

  double val2;
  double sign2;
  //Rcpp::Rcout<<"before inv" <<endl;
  //Rcpp::Rcout<<s.WtW <<endl;
  //Rcpp::Rcout<<(pi.ntree/(s.n*s2)) <<endl;
  //Rcpp::Rcout<<pi.ntree <<endl;
  //Rcpp::Rcout<<s.n <<endl;
  //Rcpp::Rcout<<s2 <<endl;
  //Rcpp::Rcout<<Prec_mi <<endl;
  // Rcpp::Rcout<< "here? 4" << endl;
  log_det(val2, sign2, Prec_mi.i());
  // Rcpp::Rcout<< "not here" << endl;

  //Rcpp::Rcout<<"after inv" <<endl;


  //tt = -0.5*log(det(precpred)) - 0.5*log(det(Prec_mi.i()));
  tt = -0.5*val1 - 0.5*val2;

  ll = -0.5*s.n*log(s2) + tt + 0.5*dotp;

  return(ll);
  }
}

mat coef_linear(tree& t, xinfo& xi, dinfo& di)
{
  double *xx;

  tree::tree_cp bn;

  mat out(di.basis_dim,di.n);

  for(size_t i=0; i<di.n; i++) {
    xx = di.x + i*di.p;
    bn = t.bn(xx,xi);
    out.col(i) = bn->getm();

  }
  return(out);
}

