#include <iostream>

#include "info.h"
#include "tree.h"
#include "bd.h"
#include "funs.h"

using std::cout;
using std::endl;

/*
notation: (as in old code): going from state x to state y (eg, incoming tree is x).

note: rather than have x and making a tree y
we just figure out what we need from x, the drawn bottom node,the drawn (v,c).
note sure what the right thing to do is.
Could make y (using a birth) and figure stuff out from y.
That is how the old code works.
*/

//
// double bd_basis(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
// {
//   tree::npv goodbots;                    //nodes we could birth at (split on)
//   double PBx = getpb(x,xi,pi,goodbots);  //prob of a birth at x
//
//   //Rcpp::Rcout << "PBx: " << PBx << endl;
//
//   // If statement for selecting birth or death proposal.
//   if(gen.uniform() < PBx) {
//
// //    Rcpp::Rcout << "birth" << endl;
//
//     //--------------------------------------------------
//     // BIRTH PROPOSAL
//     //--------------------------------------------------
//
//     //--------------------------------------------------
//     //draw proposal
//
//     //draw bottom node, choose node index ni from list in goodbots
//     size_t ni = floor(gen.uniform()*goodbots.size());
//     tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at
//
//     //draw v,  the variable
//     std::vector<size_t> goodvars; //variables nx can split on
//     getgoodvars(nx,xi,goodvars);
//     size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
//     size_t v = goodvars[vi];
//
//     //draw c, the cutpoint
//     int L,U;
//     L=0; U = xi[v].size()-1;
//     nx->rg(v,&L,&U);
//     size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points
//
//     //--------------------------------------------------
//     //compute things needed for metropolis ratio
//
//     double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
//     size_t dnx = nx->depth();
//     double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx
//
//     double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
//     if(goodvars.size()>1) { //know there are variables we could split l and r on
//       PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
//       PGry = PGly;
//     } else { //only had v to work with, if it is exhausted at either child need PG=0
//       if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
//         PGly = 0.0;
//       } else {
//         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
//       }
//       if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
//         PGry = 0.0;
//       } else {
//         PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
//       }
//     }
//
//     double PDy; //prob of proposing death at y
//     if(goodbots.size()>1) { //can birth at y because splittable nodes left
//       PDy = 1.0 - pi.pb;
//     } else { //nx was the only node you could split on
//       if((PGry==0) && (PGly==0)) { //cannot birth at y
//         PDy=1.0;
//       } else { //y can birth at either l or r
//         PDy = 1.0 - pi.pb;
//       }
//     }
//
//     double Pnogy; //death prob of choosing the nog node at y
//     size_t nnogs = x.nnogs();
//     tree::tree_cp nxp = nx->getp();
//     if(nxp==0) { //no parent, nx is the top and only node
//       Pnogy=1.0;
//     } else {
//       //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
//       if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
//         Pnogy = 1.0/nnogs;
//       } else { //if parent is not a nog, y has one more nog.
//         Pnogy = 1.0/(nnogs+1.0);
//       }
//     }
//
//     //--------------------------------------------------
//     //compute sufficient statistics
//     sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
//     getsuff_basis(x,nx,v,c,xi,di,sl,sr);
//
//     //--------------------------------------------------
//     //compute alpha
//
//     double alpha=0.0,alpha1=0.0,alpha2=0.0;
//     double lill=0.0,lilr=0.0,lilt=0.0;
//
//     if((sl.n>=5) && (sr.n>=5)) { //cludge?
//
//       //lil_basis(double n, vec sy_vec, mat WtW, double sigma, vec mu0, mat Prec0, mat Sigma0)
//
//       // lill = lil_basis(sl.n, sl.sy_vec, sl.WtW, pi);                        // USED lil(...)
//       // lilr = lil_basis(sr.n, sr.sy_vec, sr.WtW, pi);                        // USED lil(...)
//       // lilt = lil_basis(sl.n + sr.n, sl.sy_vec + sr.sy_vec, sl.WtW + sr.WtW, pi);  // USED lil(...)
//       //
//       lill = lil_basis(sl, pi);                        // USD lil(...)
//       lilr = lil_basis(sr, pi);                        // USED lil(...)
//       lilt = lil_basis(sl,sr,pi);  // USED lil(...)
//
//
//       alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
//       alpha2 = alpha1*exp(lill+lilr-lilt);
//       alpha = std::min(1.0,alpha2);
//
//     } else {
//       alpha=0.0;
//     }
//
//     /*
//     cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
//     cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
//     cout << "L,U: " << L << "," << U << endl;
//     cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
//       PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
//         ", " << Pnogy << "," << Pbotx << endl;
//     cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
//     cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
//     cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
//     cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
//     */
//
//     //--------------------------------------------------
//     //finally ready to try metrop
//     //--------------------------------------------------
//
//     vec mul, mur;
//     mul = zeros(di.basis_dim);
//     mur = zeros(di.basis_dim);
//     if(gen.uniform() < alpha) {
//
//       //--------------------------------------------------
//       // do birth:
//       // Set mul and mur to zero vectors, since we will immediately
//       // fill them in the MCMC by using drmu.  This saves computation cost.
//
//       Rcpp::Rcout << "Prebirth... " << sr.sy << " " << sr.sy_vec << endl;
//       Rcpp::Rcout << "Prebirth... " << sl.sy << " " << sl.sy_vec << endl;
//       x.birth(nx->nid(),v,c,mul,mur,sl,sr);
//       return alpha;
//     } else {
//       sinfo st = sl;
//       st.n += sr.n;
//       st.n0 += sr.n0;
//       st.sy += sr.sy;
//       st.sy_vec += sr.sy_vec;
//       st.WtW += sr.WtW;
//
//       nx->s = st;
//
//       return alpha+10;
//     }
//   } else {
// //    Rcpp::Rcout << "death" << endl;
//     //--------------------------------------------------
//     // DEATH PROPOSAL
//     //--------------------------------------------------
//
//     //--------------------------------------------------
//     //draw proposal
//
//     //draw nog node, any nog node is a possibility
//     tree::npv nognds; //nog nodes
//     x.getnogs(nognds);
//     size_t ni = floor(gen.uniform()*nognds.size());
//     tree::tree_p nx = nognds[ni]; //the nog node we might kill children at
//
//     //--------------------------------------------------
//     //compute things needed for metropolis ratio
//
//     double PGny; //prob the nog node grows
//     size_t dny = nx->depth();
//     PGny = pi.alpha/pow(1.0+dny,pi.beta);
//
//     //better way to code these two?
//     double PGlx = pgrow(nx->getl(),xi,pi);
//     double PGrx = pgrow(nx->getr(),xi,pi);
//
//     double PBy;  //prob of birth move at y
//     //if(nx->ntype()=='t') { //is the nog node nx the top node
//     if(!(nx->p)) { //is the nog node nx the top node
//       PBy = 1.0;
//     } else {
//       PBy = pi.pb;
//     }
//
//     double Pboty;  //prob of choosing the nog as bot to split on when y
//     int ngood = goodbots.size();
//     if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
//     if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
//     ++ngood;  //know you can split at nx
//     Pboty=1.0/ngood;
//
//     double PDx = 1.0-PBx; //prob of a death step at x
//     double Pnogx = 1.0/nognds.size();
//
//     //--------------------------------------------------
//     //compute sufficient statistics
//     sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
//
//     // Test if problem is getl and getr.
//     getsuff_basis(x,nx->getl(),nx->getr(),xi,di,sl,sr);
//
//     //--------------------------------------------------
//     //compute alpha
//     /*
//     double lill = lil_ts(sl.n_vec, sl.sy_vec, sl.sy2, pi.sigma, pi.mu0, pi.Prec0);
//     double lilr = lil_ts(sr.n_vec, sr.sy_vec, sr.sy2, pi.sigma, pi.mu0, pi.Prec0);
//     double lilt = lil_ts(sl.n_vec + sr.n_vec, sl.sy_vec + sr.sy_vec, sl.sy2 + sr.sy2, pi.sigma, pi.mu0, pi.Prec0);
//     */
//     // double lill = lil_basis(sl.n, sl.sy_vec, sl.WtW, pi);                        // USED lil(...)
//     // double lilr = lil_basis(sr.n, sr.sy_vec, sr.WtW, pi);                        // USED lil(...)
//     // double lilt = lil_basis(sl.n + sr.n, sl.sy_vec + sr.sy_vec, sl.WtW + sr.WtW, pi);  // USED lil(...)
//     //
//     //
//
//     double lill = lil_basis(sl, pi);                        // USED lil(...)
//     double lilr = lil_basis(sr, pi);                        // USED lil(...)
//     double lilt = lil_basis(sl,sr, pi);  // USED lil(...)
//
//
//     double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
//     double alpha2 = alpha1*exp(lilt - lill - lilr);
//     double alpha = std::min(1.0,alpha2);
//
//     /*
//     cout << "death prop: " << nx->nid() << endl;
//     cout << "nognds.size(), ni, nx: " << nognds.size() << ", " << ni << ", " << nx << endl;
//     cout << "depth of nog node: " << dny << endl;
//     cout << "PGny: " << PGny << endl;
//     cout << "PGlx: " << PGlx << endl;
//     cout << "PGrx: " << PGrx << endl;
//     cout << "PBy: " << PBy << endl;
//     cout << "Pboty: " << Pboty << endl;
//     cout << "PDx: " << PDx << endl;
//     cout << "Pnogx: " << Pnogx << endl;
//     cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
//     cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
//     cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
//     cout << "sigma: " << pi.sigma << endl;
//     cout << "tau: " << pi.tau << endl;
//     cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
//     */
//
//     //--------------------------------------------------
//     //finally ready to try metrop
//
//     // All values wrong, but ok bc get reset immediately in the draw mu.
//     if(gen.uniform()<alpha) {
//       //draw mu for nog (which will be bot)
//       vec mu;
//       mu = zeros(di.basis_dim);
//
//       sinfo st = sl;
//       st.n += sr.n;
//       st.n0 += sr.n0;
//       st.sy += sr.sy;
//       st.sy_vec += sr.sy_vec;
//       st.WtW += sr.WtW;
//
//       //do death
//       x.death(nx->nid(),mu, st);
//
//       return -alpha;
//     } else {
//
//       nx->getl()->s = sl;
//       nx->getr()->s = sr;
//
//       return -alpha-10;
//     }
//   }
// }



double bd_basis(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen, std::vector<tree::tree_cp>& node_pointers)
{
  tree::npv goodbots;                    //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots);  //prob of a birth at x

  //Rcpp::Rcout << "PBx: " << PBx << endl;

  // If statement for selecting birth or death proposal.
  if(gen.uniform() < PBx) {

    //    Rcpp::Rcout << "birth" << endl;

    //--------------------------------------------------
    // BIRTH PROPOSAL
    //--------------------------------------------------

    //--------------------------------------------------
    //draw proposal

    //draw bottom node, choose node index ni from list in goodbots
    size_t ni = floor(gen.uniform()*goodbots.size());
    tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on
    getgoodvars(nx,xi,goodvars);
    //size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable

    // DART
    size_t vi;
    if(pi.dart) {
      std::vector<double> wt(goodvars.size());
      for(size_t j=0; j<goodvars.size(); ++j) {
        wt[j] = log(pi.var_probs[goodvars[j]]);
      }
      //vi = goodvars[rdisc_log_inplace(wt)]; //sample variable using the log-weights
      vi = rdisc_log_inplace(wt); //sample variable using the log-weights
    } else {
      vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
    }

    size_t v = goodvars[vi];

    //draw c, the cutpoint
    int L,U;
    L=0; U = xi[v].size()-1;
    nx->rg(v,&L,&U);
    size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
    size_t dnx = nx->depth();
    double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

    double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
    if(goodvars.size()>1) { //know there are variables we could split l and r on
      PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
      PGry = PGly;
    } else { //only had v to work with, if it is exhausted at either child need PG=0
      if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else {
        PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
      if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
        PGry = 0.0;
      } else {
        PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
    }

    double PDy; //prob of proposing death at y
    if(goodbots.size()>1) { //can birth at y because splittable nodes left
      PDy = 1.0 - pi.pb;
    } else { //nx was the only node you could split on
      if((PGry==0) && (PGly==0)) { //cannot birth at y
        PDy=1.0;
      } else { //y can birth at either l or r
        PDy = 1.0 - pi.pb;
      }
    }

    double Pnogy; //death prob of choosing the nog node at y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();
    if(nxp==0) { //no parent, nx is the top and only node
      Pnogy=1.0;
    } else {
      //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
      if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else { //if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs+1.0);
      }
    }

    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))

    sl.sy = 0;
    sl.n0 = 0.0;
    sl.sy_vec.zeros(di.basis_dim);
    sl.WtW.zeros(di.basis_dim, di.basis_dim);
    sr.sy = 0;
    sr.n0 = 0.0;
    sr.sy_vec.zeros(di.basis_dim);
    sr.WtW.zeros(di.basis_dim, di.basis_dim);

    // void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, sinfo& sl, sinfo& sr)
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff
    //allsuff_basis(t,xi,di,bnv0,sv0);
    allsuff_basis_birth(x,nx,v,c,xi,di,bnv,sv,sl,sr,node_pointers);

    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }

    //--------------------------------------------------
    //compute alpha

    double alpha=0.0,alpha1=0.0,alpha2=0.0;
    double lill=0.0,lilr=0.0,lilt=0.0;

    if((sl.n>=5) && (sr.n>=5)) { //cludge?

      //lil_basis(double n, vec sy_vec, mat WtW, double sigma, vec mu0, mat Prec0, mat Sigma0)

      // lill = lil_basis(sl.n, sl.sy_vec, sl.WtW, pi);                        // USED lil(...)
      // lilr = lil_basis(sr.n, sr.sy_vec, sr.WtW, pi);                        // USED lil(...)
      // lilt = lil_basis(sl.n + sr.n, sl.sy_vec + sr.sy_vec, sl.WtW + sr.WtW, pi);  // USED lil(...)
      //
      lill = lil_basis(sl, pi);                        // USD lil(...)
      lilr = lil_basis(sr, pi);                        // USED lil(...)
      lilt = lil_basis(sl,sr,pi);  // USED lil(...)
      // Rcpp::Rcout<< "Basis regular "<< endl;
      // Rcpp::Rcout<< "v: " << v << endl;
      // Rcpp::Rcout<< "lill: " << lill << endl;
      // Rcpp::Rcout<< "lilr: " << lilr << endl;
      // Rcpp::Rcout<< "lilt: " << lilt << endl;


      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1*exp(lill+lilr-lilt);
      alpha = std::min(1.0,alpha2);

    } else {
      alpha=0.0;
    }

    /*
    cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
    cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
    cout << "L,U: " << L << "," << U << endl;
    cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
      PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
        ", " << Pnogy << "," << Pbotx << endl;
    cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
    cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
    cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
    cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
    */

    //--------------------------------------------------
    //finally ready to try metrop
    //--------------------------------------------------

    vec mul, mur;
    mul = zeros(di.basis_dim);
    mur = zeros(di.basis_dim);
    if(gen.uniform() < alpha) {

      //--------------------------------------------------
      // do birth:
      // Set mul and mur to zero vectors, since we will immediately
      // fill them in the MCMC by using drmu.  This saves computation cost.

      //Rcpp::Rcout << "Prebirth... " << sr.sy << " " << sr.sy_vec << endl;
      //Rcpp::Rcout << "Prebirth... " << sl.sy << " " << sl.sy_vec << endl;

      //Rcpp::Rcout << "Birthing";


      //Rcpp::Rcout << " before birth" << endl;
      //Rcpp::Rcout << "L" << endl << sl.WtW << endl;
      //Rcpp::Rcout << "R" << endl << sr.WtW << endl;
      x.birth(nx->nid(),v,c,mul,mur,sl,sr);

      double *xx;        //current x
      for(size_t i=0;i<di.n;i++) {
        xx = di.x + i*di.p;
        bool in_candidate_nog, left, right;
        in_candidate_nog = (node_pointers[i] == nx);
        left  = in_candidate_nog & (xx[v] < xi[v][c]);
        right = in_candidate_nog & !(xx[v] < xi[v][c]);
        if(left)  node_pointers[i] = nx->getl();
        if(right) node_pointers[i] = nx->getr();
      }

      //Rcpp::Rcout << " after birth" << endl;
      //Rcpp::Rcout << "L" << endl << sl.WtW << endl;
      //Rcpp::Rcout << "R" << endl << sr.WtW << endl;
      //Rcpp::Rcout << " Birthed";
      return alpha;
    } else {
      /*
      sinfo st = sl;
      st.n += sr.n;
      st.n0 += sr.n0;
      st.sy += sr.sy;
      st.sy_vec += sr.sy_vec;
      st.WtW += sr.WtW;

      nx->s = st;
       */

      return alpha+10;
    }
  } else {
    //    Rcpp::Rcout << "death" << endl;
    //--------------------------------------------------
    // DEATH PROPOSAL
    //--------------------------------------------------

    //--------------------------------------------------
    //draw proposal

    //draw nog node, any nog node is a possibility
    tree::npv nognds; //nog nodes
    x.getnogs(nognds);
    size_t ni = floor(gen.uniform()*nognds.size());
    tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double PGny; //prob the nog node grows
    size_t dny = nx->depth();
    PGny = pi.alpha/pow(1.0+dny,pi.beta);

    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,pi);
    double PGrx = pgrow(nx->getr(),xi,pi);

    double PBy;  //prob of birth move at y
    //if(nx->ntype()=='t') { //is the nog node nx the top node
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = pi.pb;
    }

    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;

    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();

    //--------------------------------------------------
    //compute sufficient statistics
    //sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))

    // Test if problem is getl and getr.
    //getsuff_basis(x,nx->getl(),nx->getr(),xi,di,sl,sr);

    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))

    // void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, sinfo& sl, sinfo& sr)
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff
    //allsuff_basis(t,xi,di,bnv0,sv0);
    allsuff_basis(x,xi,di,bnv,sv,node_pointers);

    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }

    sl = nx->getl()->s;
    sr = nx->getr()->s;

    // call new function merge_suffstats with arguments like allsuff_basis but also
    // sl->getl() and sl->getr(). This function will populate another suff stat object
    // called st (for s total) corresponding to the node nx after death

    double lill = lil_basis(sl, pi);                        // USED lil(...)
    double lilr = lil_basis(sr, pi);                        // USED lil(...)
    double lilt = lil_basis(sl,sr, pi);  // this will become lilt = lil_basis(st, pi);


    double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    double alpha2 = alpha1*exp(lilt - lill - lilr);
    double alpha = std::min(1.0,alpha2);

    //--------------------------------------------------
    //finally ready to try metrop

    // All values wrong, but ok bc get reset immediately in the draw mu.
    if(gen.uniform()<alpha) {

      tree::tree_cp lptr = nx->getl();
      tree::tree_cp rptr = nx->getr();

      for(size_t i=0;i<di.n;i++) {
        if((node_pointers[i]==lptr) | (node_pointers[i]==rptr)) {
          node_pointers[i] = nx;
        }
      }

      //Rcpp::Rcout << "Deathing";

      //draw mu for nog (which will be bot)
      vec mu;
      mu = zeros(di.basis_dim);

      sinfo st = sl;
      st.n += sr.n;
      st.n0 += sr.n0;
      st.sy += sr.sy;
      st.sy_vec += sr.sy_vec;
      st.WtW += sr.WtW;

      //do death
      x.death(nx->nid(),mu, st);

      return -alpha;
    } else {

      //nx->getl()->s = sl;
      //nx->getr()->s = sr;

      return -alpha-10;
    }
  }
}


double bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{
   tree::npv goodbots;                    //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots);  //prob of a birth at x

   //Rcpp::Rcout << "PBx: " << PBx << endl;

   // If statement for selecting birth or death proposal.
   if(gen.uniform() < PBx) {
      //--------------------------------------------------
      // BIRTH PROPOSAL
      //--------------------------------------------------

      //--------------------------------------------------
      //draw proposal

      //draw bottom node, choose node index ni from list in goodbots
      size_t ni = floor(gen.uniform()*goodbots.size());
      tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

      //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      size_t v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_cp nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
         if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
           Pnogy = 1.0/(nnogs+1.0);
         }
      }

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuff_ts(x,nx,v,c,xi,di,sl,sr, di.tlen);

      //--------------------------------------------------
      //compute alpha

      double alpha=0.0,alpha1=0.0,alpha2=0.0;
      double lill=0.0,lilr=0.0,lilt=0.0;

      //Rcpp::Rcout << "birth alpha sl.n: " << sl.n << ", sr.n: " << sr.n << endl;
      if((sl.n0>=5) && (sr.n0>=5)) { //cludge?
      //if((sl.n>=1) && (sr.n>=1)) { //cludge?

         lill = lil_ts(sl.n_vec, sl.sy_vec, sl.sy2, pi.sigma, pi.mu0, pi.Prec0);                        // USED lil(...)
         lilr = lil_ts(sr.n_vec, sr.sy_vec, sr.sy2, pi.sigma, pi.mu0, pi.Prec0);                        // USED lil(...)
         lilt = lil_ts(sl.n_vec + sr.n_vec, sl.sy_vec + sr.sy_vec, sl.sy2 + sr.sy2, pi.sigma, pi.mu0, pi.Prec0);  // USED lil(...)

         alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
         alpha2 = alpha1*exp(lill+lilr-lilt);
         alpha = std::min(1.0,alpha2);

       } else {
         alpha=0.0;
      }

      /*
      cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
      cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
      cout << "L,U: " << L << "," << U << endl;
      cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
         PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
         ", " << Pnogy << "," << Pbotx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop
      //--------------------------------------------------

      vec mul, mur;
      mul = zeros(di.tlen);
      mur = zeros(di.tlen);
      if(gen.uniform() < alpha) {

         //--------------------------------------------------
         // do birth:
         // Set mul and mur to zero vectors, since we will immediately
         // fill them in the MCMC by using drmu.  This saves computation cost.

			   x.birth(nx->nid(),v,c,mul,mur);
         return alpha;
      } else {
         return alpha+10;
      }
   } else {
      //--------------------------------------------------
      // DEATH PROPOSAL
      //--------------------------------------------------

       //--------------------------------------------------
      //draw proposal

      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size());
      tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.beta);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,pi);
      double PGrx = pgrow(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
      //if(nx->ntype()=='t') { //is the nog node nx the top node
      if(!(nx->p)) { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl(di.tlen),sr(di.tlen); //sl for left from nx and sr for right from nx (using rule (v,c))

      // Test if problem is getl and getr.
      getsuff_ts(x,nx->getl(),nx->getr(),xi,di,sl,sr,di.tlen);

      //--------------------------------------------------
      //compute alpha

      double lill = lil_ts(sl.n_vec, sl.sy_vec, sl.sy2, pi.sigma, pi.mu0, pi.Prec0);
      double lilr = lil_ts(sr.n_vec, sr.sy_vec, sr.sy2, pi.sigma, pi.mu0, pi.Prec0);
      double lilt = lil_ts(sl.n_vec + sr.n_vec, sl.sy_vec + sr.sy_vec, sl.sy2 + sr.sy2, pi.sigma, pi.mu0, pi.Prec0);

      double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
      double alpha2 = alpha1*exp(lilt - lill - lilr);
      double alpha = std::min(1.0,alpha2);

      /*
      cout << "death prop: " << nx->nid() << endl;
      cout << "nognds.size(), ni, nx: " << nognds.size() << ", " << ni << ", " << nx << endl;
      cout << "depth of nog node: " << dny << endl;
      cout << "PGny: " << PGny << endl;
      cout << "PGlx: " << PGlx << endl;
      cout << "PGrx: " << PGrx << endl;
      cout << "PBy: " << PBy << endl;
      cout << "Pboty: " << Pboty << endl;
      cout << "PDx: " << PDx << endl;
      cout << "Pnogx: " << Pnogx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "sigma: " << pi.sigma << endl;
      cout << "tau: " << pi.tau << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop

      // All values wrong, but ok bc get reset immediately in the draw mu.
      if(gen.uniform()<alpha) {
         //draw mu for nog (which will be bot)
         vec mu;
         mu = zeros(di.tlen);

         //do death
			x.death(nx->nid(),mu);
         return -alpha;
      } else {
         return -alpha-10;
      }
   }
}

// Het variances.
double bdhet(tree& x, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen)
{
   tree::npv goodbots;                    //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots);  //prob of a birth at x
   //Rcpp::Rcout << "PBx: " << PBx << endl;

   // If statement for selecting birth or death proposal.
   if(gen.uniform() < PBx) {

      //--------------------------------------------------
      // BIRTH PROPOSAL
      //--------------------------------------------------

       //--------------------------------------------------
      //draw proposal

      //draw bottom node, choose node index ni from list in goodbots
      size_t ni = floor(gen.uniform()*goodbots.size());
      tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

       //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      size_t v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_cp nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
         if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
            Pnogy = 1.0/(nnogs+1.0);
         }
      }

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuffhet_ts(x,nx,v,c,xi,di,phi,sl,sr,di.tlen);

      //--------------------------------------------------
      //compute alpha

      double alpha=0.0,alpha1=0.0,alpha2=0.0;
      double lill=0.0,lilr=0.0,lilt=0.0;

      //if((sl.n>=5) && (sr.n>=5)) {
      if((abs(sl.n)>=2) && (abs(sr.n)>=2)) { //cludge?

         lill = lilhet_ts(sl.n0, sl.n, sl.n_vec, sl.sy_vec, sl.sy2, pi.mu0, pi.Prec0);
         lilr = lilhet_ts(sr.n0, sr.n, sr.n_vec, sr.sy_vec, sr.sy2, pi.mu0, pi.Prec0);
         lilt = lilhet_ts(sl.n0+sr.n0, sl.n+sr.n, sl.n_vec+sr.n_vec,
                          sl.sy_vec+sr.sy_vec,
                          sl.sy2+sr.sy2, pi.mu0, pi.Prec0);

         alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
         alpha2 = alpha1*exp(lill+lilr-lilt);
         alpha = std::min(1.0,alpha2);
      } else {
         alpha=0.0;
      }

      /*
      cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
      cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
      cout << "L,U: " << L << "," << U << endl;
      cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
         PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
            ", " << Pnogy << "," << Pbotx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop
      //--------------------------------------------------

      vec mul, mur;
      mul = zeros(di.tlen);
      mur = zeros(di.tlen);
      if(gen.uniform() < alpha) {

         //--------------------------------------------------
         // do birth:
         // Set mul and mur to zero vectors, since we will immediately
         // fill them in the MCMC by using drmu.  This saves computation cost.

         x.birth(nx->nid(),v,c,mul,mur);
         return alpha;
      } else {
         return alpha+10;
      }
   } else {
      //--------------------------------------------------
      // DEATH PROPOSAL
      //--------------------------------------------------

      //--------------------------------------------------
      //draw proposal

      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size());
      tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.beta);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,pi);
      double PGrx = pgrow(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
      //if(nx->ntype()=='t') { //is the nog node nx the top node
      if(!(nx->p)) { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl(di.tlen),sr(di.tlen); //sl for left from nx and sr for right from nx (using rule (v,c))

      // Test if problem is getl and getr.
      getsuffhet_ts(x,nx->getl(),nx->getr(),xi,di,phi,sl,sr,di.tlen);

      //--------------------------------------------------
      //compute alpha

      double lill = lilhet_ts(sl.n0, sl.n, sl.n_vec, sl.sy_vec, sl.sy2, pi.mu0, pi.Prec0);
      double lilr = lilhet_ts(sr.n0, sr.n, sr.n_vec, sr.sy_vec, sr.sy2, pi.mu0, pi.Prec0);
      double lilt = lilhet_ts(sl.n0+sr.n0, sl.n+sr.n, sl.n_vec+sr.n_vec,
                              sl.sy_vec+sr.sy_vec,
                              sl.sy2+sr.sy2, pi.mu0, pi.Prec0);

      double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
      double alpha2 = alpha1*exp(lilt - lill - lilr);
      double alpha = std::min(1.0,alpha2);

      /*
      cout << "death prop: " << nx->nid() << endl;
      cout << "nognds.size(), ni, nx: " << nognds.size() << ", " << ni << ", " << nx << endl;
      cout << "depth of nog node: " << dny << endl;
      cout << "PGny: " << PGny << endl;
      cout << "PGlx: " << PGlx << endl;
      cout << "PGrx: " << PGrx << endl;
      cout << "PBy: " << PBy << endl;
      cout << "Pboty: " << Pboty << endl;
      cout << "PDx: " << PDx << endl;
      cout << "Pnogx: " << Pnogx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "sigma: " << pi.sigma << endl;
      cout << "tau: " << pi.tau << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop

      // All values wrong, but ok bc get reset immediately in the draw mu.
      if(gen.uniform()<alpha) {
         //draw mu for nog (which will be bot)
         vec mu;
         mu = zeros(di.tlen);

         //do death
         x.death(nx->nid(),mu);

         //return true;
         return -alpha;
      } else {
         return -alpha-10;
      }
   }
}

double bd_linear(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen, std::vector<tree::tree_cp>& node_pointers)
{
  tree::npv goodbots;                    //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots);  //prob of a birth at x
  
  // Rcpp::Rcout<< "PBx: " << PBx << endl;
  
  // If statement for selecting birth or death proposal.
  if(gen.uniform() < PBx) {
    // Rcpp::Rcout << "Inside Birth" << endl;
    //--------------------------------------------------
    //--------------------------------------------------
    // BIRTH PROPOSAL
    //--------------------------------------------------
    //--------------------------------------------------
    //draw proposal

    //draw bottom node, choose node index ni from list in goodbots
    size_t ni = floor(gen.uniform()*goodbots.size());
    tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on
    getgoodvars(nx,xi,goodvars);
    //size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable

    // DART
    size_t vi;
    if(pi.dart) {
      std::vector<double> wt(goodvars.size());
      for(size_t j=0; j<goodvars.size(); ++j) {
        wt[j] = log(pi.var_probs[goodvars[j]]);
      }
      //vi = goodvars[rdisc_log_inplace(wt)]; //sample variable using the log-weights
      vi = rdisc_log_inplace(wt); //sample variable using the log-weights
    } else {
      vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
    }

    size_t v = goodvars[vi];

    //draw c, the cutpoint
    int L,U;
    L=0; U = xi[v].size()-1;
    nx->rg(v,&L,&U);
    size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
    size_t dnx = nx->depth();
    double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

    double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
    if(goodvars.size()>1) { //know there are variables we could split l and r on
      PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
      PGry = PGly;
    } else { //only had v to work with, if it is exhausted at either child need PG=0
      if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else {
        PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
      if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
        PGry = 0.0;
      } else {
        PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
    }
    double PDy; //prob of proposing death at y
    if(goodbots.size()>1) { //can birth at y because splittable nodes left
      PDy = 1.0 - pi.pb;
    } else { //nx was the only node you could split on
      if((PGry==0) && (PGly==0)) { //cannot birth at y
        PDy=1.0;
      } else { //y can birth at either l or r
        PDy = 1.0 - pi.pb;
      }
    }

    double Pnogy; //death prob of choosing the nog node at y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();
    if(nxp==0) { //no parent, nx is the top and only node
      Pnogy=1.0;
    } else {
      //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
      if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else { //if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs+1.0);
      }
    }

    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl(di.basis_dim),sr(di.basis_dim),st(di.basis_dim);; //sl for left from nx and sr for right from nx (using rule (v,c))

    sl.sy = 0;
    sl.n0 = 0.0;
    sl.sy_vec.zeros(di.basis_dim);
    sl.WtW.zeros(di.basis_dim, di.basis_dim);
    sr.sy = 0;
    sr.n0 = 0.0;
    sr.sy_vec.zeros(di.basis_dim);
    sr.WtW.zeros(di.basis_dim, di.basis_dim);

    // void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, sinfo& sl, sinfo& sr)
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff
    //allsuff_basis(t,xi,di,bnv0,sv0);
    allsuff_linear_birth(x,nx,v,c,xi,di,bnv,sv,sl,sr,node_pointers);
    // Rcpp::Rcout<< "After allsuff" << endl;
    // Rcpp::Rcout<< "sr: " << sr.WtW << endl;
    // Rcpp::Rcout<< "sl: " << sl.WtW << endl;
    // Rcpp::Rcout<< "sr.n: " << sr.n << endl;
    // Rcpp::Rcout<< "sl.n: " << sl.n << endl;
    
    
    
    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }

    //The variable to update the suff. stats
    st = nx->s;
    size_t vp;
    if(nx->nid()==1){
      vp = nx->getv();
    }else{
      vp = (*(nx->getp())).getv();
    }

    //--------------------------------------------------
    //compute alpha

    double alpha=0.0,alpha1=0.0,alpha2=0.0;
    double lill=0.0,lilr=0.0,lilt=0.0;

    //////
    double ldet_right;
    double ldet_left;
    double sign;
    log_det(ldet_right,sign,sr.WtW);
    log_det(ldet_left,sign,sl.WtW);
    //////
    // Rcpp::Rcout<< "ldet_left: " << ldet_left << endl;
    // Rcpp::Rcout<< "ldet_right: " << ldet_right << endl;

    //if(((sl.n>=5) && (sr.n>=5))&&
       //((ldet_left>-5)&&(ldet_right>-5))) { //cludge?

    if(((sl.n>=5) && (sr.n>=5))){
       
      //The likelihood must be updated at each draw
      //because the prior m_i is always changing
      // Rcpp::Rcout<< "v: " << v << endl;
      // Rcpp::Rcout<< "vp: " << vp << endl;
      // 
      
      lill = lil_linear(sl, pi, v);
      // Rcpp::Rcout<< "lill: " << lill << endl;
      lilr = lil_linear(sr, pi, v);
      // Rcpp::Rcout<< "lilr: " << lilr << endl;
      lilt = lil_linear(st, pi, vp);
      // Rcpp::Rcout<< "lilt: " << lilt << endl;


      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1*exp(lill+lilr-lilt);
      alpha = std::min(1.0,alpha2);

    } else {
      alpha=0.0;
    }

    //--------------------------------------------------
    //finally ready to try metrop
    //--------------------------------------------------

    vec mul, mur;
    mul = zeros(di.basis_dim);
    mur = zeros(di.basis_dim);
    if(gen.uniform() < alpha) {
      //Here Birth Proposal has been accepted


      ////////////////////////////////////////////////
      //This section will update the m_is when death occurs
      ////////////////////////////////////////////////

      //Verify if the rule already exists on the tree
      //arma::vec maux;
      //arma::vec nv;
      //size_t vv;
      //nv = zeros(2);
      //maux = zeros(di.p);

      //Get rule for every single node on the tree
      //for(size_t kk=0; kk<bnv.size(); ++kk) {
      //  if(bnv[kk]->nid()==1){
      //    vv = bnv[kk]->getv();
      //  }else{
      //    vv = (*(bnv[kk]->getp())).getv();
      //  }
      //  maux[vv] += 1;
      //}

      //Getting the new set of rules
      //if(nx->nid()==1){
      //  nv[0] = nx->getv();
      //  nv[1] = v;
      //}else{
      //  nv[0] = (*(nx->getp())).getv();
      //  nv[1] = v;
      //}

      //One node had this rule. Now it has two children
      //maux[(nv[0])] -= 1; //Update the usage of the variable
      //if(maux[(nv[0])] == 0){
      //  //Update m_i
      //  pi.mi[(nv[0])] -= 1;
      //}
      //Two new nodes have this new rule, then add 2
      //maux[(nv[1])] += 2;

      //If the rule has only be used twice, then is is a new rule
      //Update m_i
      //if(maux[(nv[1])] == 2){
      //  pi.mi[(nv[1])] += 1;
      //}
      ////////////////////////////////////////////////



      //--------------------------------------------------
      // do birth:
      // Set mul and mur to zero vectors, since we will immediately
      // fill them in the MCMC by using drmu.  This saves computation cost.
      x.birth(nx->nid(),v,c,mul,mur,sl,sr);


      //Update the pointers to final nodes
      double *xx;        //current x
      for(size_t i=0;i<di.n;i++) {
        xx = di.x + i*di.p;
        bool in_candidate_nog, left, right;
        in_candidate_nog = (node_pointers[i] == nx);
        left  = in_candidate_nog & (xx[v] < xi[v][c]);
        right = in_candidate_nog & !(xx[v] < xi[v][c]);
        if(left)  node_pointers[i] = nx->getl();
        if(right) node_pointers[i] = nx->getr();
      }

      return alpha;
    } else {

      return alpha+10;
    }
  } else {
    // Rcpp::Rcout << "Inside Death" << endl;
    //--------------------------------------------------
    //--------------------------------------------------
    // DEATH PROPOSAL
    //--------------------------------------------------
    //--------------------------------------------------
    //draw proposal

    //draw nog node, any nog node is a possibility
    tree::npv nognds; //nog nodes
    x.getnogs(nognds);
    size_t ni = floor(gen.uniform()*nognds.size());
    tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double PGny; //prob the nog node grows
    size_t dny = nx->depth();
    PGny = pi.alpha/pow(1.0+dny,pi.beta);

    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,pi);
    double PGrx = pgrow(nx->getr(),xi,pi);

    double PBy;  //prob of birth move at y
    //if(nx->ntype()=='t') { //is the nog node nx the top node
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = pi.pb;
    }

    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;

    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();

    //--------------------------------------------------
    //Objects that will receive the suff. stats
    sinfo sl(di.basis_dim),sr(di.basis_dim),st(di.basis_dim);

    //This object will hold the bottom nodes
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff

    //Variable used in case nx is root.
    size_t newvar = floor(gen.uniform()*(di.p));

    //Function to calculate sufficient statistics
    allsuff_linear_death(x,nx,xi,di,bnv,sv,st,node_pointers, newvar);

    //Getting sufficient statistics for every node
    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }

    //Getting the sufficient statistics for the leaves
    sl = nx->getl()->s;
    sr = nx->getr()->s;

    //Get the variables that will be used on the likelihood function
    size_t v = nx->getv(); //The old nodes rules
    size_t vp;
    //The rule of the new node
    if(nx->nid()==1){
      vp = newvar; //Use this draw if it is a root node
    }else{
      vp = (*(nx->getp())).getv(); //Otherwise, use parent
    }

    //The likelihood must be updated at each draw
    //because the prior m_i is always changing
    double lill = lil_linear(sl, pi, v);
    double lilr = lil_linear(sr, pi, v);
    double lilt = lil_linear(st, pi, vp);
    
    // Rcpp::Rcout<< "v: " << v << endl;
    // Rcpp::Rcout<< "lill: " << lill << endl;
    // Rcpp::Rcout<< "lilr: " << lilr << endl;
    // Rcpp::Rcout<< "lilt: " << lilt << endl;


    double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    double alpha2 = alpha1*exp(lilt - lill - lilr);
    double alpha = std::min(1.0,alpha2);

    //--------------------------------------------------
    //finally ready to try metrop

    // All values wrong, but ok bc get reset immediately in the draw mu.
    if(gen.uniform()<alpha) {
      //Here Death Proposal has been accepted


      ////////////////////////////////////////////////
      //This section will update the m_is when death occurs
      ////////////////////////////////////////////////

      //Verify if the rule already exists on the tree
      //arma::vec maux; //Vector that will count the nodes rules on this tree
      //arma::vec nv; //Vector os new variables
      //size_t vv; //Iterator
      //nv = zeros(2); // Initializing nv
      //maux = zeros(di.p); //Initializing maux

      //Since the object bnv already exists, I will use it to
      //find the variables
      //for(size_t kk=0; kk<bnv.size(); ++kk) {
      //  if(bnv[kk]->nid()==1){
      //    //Special case: If the node is root, then the variable
      //    //is stored in the same node
      //    vv = bnv[kk]->getv();
      //  }else{
      //    //Otherwise, is on the parent
      //    vv = (*(bnv[kk]->getp())).getv();
      //  }
      //  //Count how many times each variable is used on the tree
      //  maux[vv] += 1;
      //}

      //Get the new variables
      //if(nx->nid()==1){
      //  nv[0] = nx->getv();
      //  nv[1] = newvar;
      //}else{
      //  nv[0] = nx->getv();
      //  nv[1] = (*(nx->getp())).getv();
      //}

      //Two node had the same rule, then -2 because they collapsed
      //maux[(nv[0])] -= 2;

      //If there are no more nodes with that rule, update m_i
      //if(maux[(nv[0])] == 0){
      //  pi.mi[(nv[0])] -= 1;
      //}

      //One node has this new rule, then +1
      //maux[(nv[1])] += 1;

      //If there is only one node with that rule, then it is a new rule
      //Therefore, m_i must be updated
      //if(maux[(nv[1])] == 1){
      //  pi.mi[(nv[1])] += 1;
      //}
      ////////////////////////////////////////////////



      tree::tree_cp lptr = nx->getl();
      tree::tree_cp rptr = nx->getr();
      //Update the pointers on root nodes
      for(size_t i=0;i<di.n;i++) {
        if((node_pointers[i]==lptr) | (node_pointers[i]==rptr)) {
          node_pointers[i] = nx;
        }
      }

      //Rcpp::Rcout << "Deathing";

      //draw mu for nog (which will be bot)
      vec mu;
      mu = zeros(di.basis_dim);


      //do death
      x.death(nx->nid(),mu, st);

      //After the nodes are collapsed, if it is a root node,
      //then change the variable at the node
      if(nx->nid()==1){
        nx->setv(newvar);
      }

      return -alpha;
    } else {

      return -alpha-10;
    }
  }
}
