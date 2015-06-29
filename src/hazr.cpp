#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector hazr(double t, double m, double bm, double bt, NumericMatrix dat){
  // MAKE SURE dat IS SORTED BY TIME OR THIS FUNCTION WILL NOT WORK!!!
  // bm should be on marker scale
  // bt should be on time scale
  // ---------------------------------------------------------------------------
  NumericVector out(2);
  
  // Collect subset of observations in marker trip in mstrip
  // ---------------------------------------------------------------------------
  int j = 0; // number of values in marker strip
  NumericMatrix mstrip(dat.nrow(), 2); // put obsevations in marker strip here
  for (int i = 0; i < dat.nrow(); i++){
    // populate mstrip with observations
    if ((dat(i, 2) < m+bm) && (dat(i,2) > m-bm)){
      mstrip(j, 0) = dat(i, 0);
      mstrip(j, 1) = dat(i, 1);
      j += 1;
    }
  }
  if (j == 0) return(out);
  
  // - Do the local-in-marker hazard and N-A calculations
  // ---------------------------------------------------------------------------
  NumericVector na(j); // create vector for local N-A estimates
  double nac = 0;      // holder for N-A contributions
  na[0] = mstrip(0, 1)/j; // N-A estimate at earliest time in strip
  if (fabs(t - mstrip(0, 0)) <= bt){
    // if first value is close enough to t, include in hazard calculation
    out[0] = na[0];
    out[1] = na[0]; // set N-A estimator
  }
  for (int i = 1; i < j; i++){
    // include rest of qualified observations into hazard calculation
    // calculate remaining N-A estimates
    nac = (mstrip(i, 1)/(j - i)); // censoring indicat. divided by risk set size
    na[i] = na[i-1] + nac;
    if (fabs(t - mstrip(i, 0)) <= bt){
      // if close enough to t, include in hazard calculation
      out[0] += nac;
    }
    if (t >= mstrip(i, 0)){
      // if we have reached t, return this for N-A estimate
      out[1] = na[i];
    }
  }
  out[0] = out[0]/(2*bt);
  return(out);
}

// [[Rcpp::export]]
NumericVector nnhazr(double t, double m, double bm, double bt, NumericMatrix dat){
  // MAKE SURE dat IS SORTED BY TIME OR THIS FUNCTION WILL NOT WORK!!!
  // bm should be on the marker scale
  // bt should be number of nearest neighbors for time
  // The implied bandwidth this function uses is unconventional: independently
  //   calculates a "left" and "right" bandwidth, and then takes their sum
  // ---------------------------------------------------------------------------
  NumericVector out(6); // first two slots are conditional hazard rate and 
                        // N-A estimator; others are just diagnostic output
  //  if ((bt < 0) || (bt>1)){
  //    return(out);
  //  }
  
  // Collect subset of observations in marker trip in mstrip
  // ---------------------------------------------------------------------------
  int j = 0; // number of values in marker strip
  NumericMatrix mstrip(dat.nrow(), 2); // put obsevations in marker strip here
  for (int i = 0; i < dat.nrow(); i++){
    // populate mstrip with observations
    if ((dat(i, 2) < m+bm) && (dat(i,2) > m-bm)){
      mstrip(j, 0) = dat(i, 0);
      mstrip(j, 1) = dat(i, 1);
      j += 1;
    }
  }
  bt = bt/j;
  // Calculate K-M estimates and put in km
  // Needed for nearest neighbors calculation of hazard function
  // ---------------------------------------------------------------------------
  NumericVector km(j); // create vector for local K-M estimates
  km[0] = (j - mstrip(0, 1))/j; // K-M estimate at earliest time in strip
  double tp = 1.0; // find percentile for t within strip (for nearest neighbors)
  if (t >= mstrip(0, 0)){
    tp = km[0];
  }
  for (int i = 1; i < j; i++){
    // calculate remaining K-M estimates
    // meanwhile find out percentile for t within strip
    km[i] = km[i-1] * (j - i - mstrip(i, 1))/(j - i);
    if (t >= mstrip(i, 0)){
      tp = km[i];
    }
  }
  
  // - Do the local-in-marker hazard and N-A calculations
  // ---------------------------------------------------------------------------
  NumericVector na(j); // create vector for local N-A estimates
  double nac = 0;      // holder for N-A contributions
  double hrc = 0;      // counter for number of observations in hazard estimator
  double btl = 0;      // implied left bw from NN
  double btr = 0;      // implied right bw from nn
  na[0] = mstrip(0, 1)/j; // N-A estimate at earliest time in strip
  if (fabs(tp - km[0]) <= bt){
    // if first value is close enough to t, include in hazard calculation
    out[0] = na[0];
    if ((tp - km[0]) <= 0){
      // if t percentile is greater than or equal to the percentile of the first
      // time value in the strip then update left bandwidth value
      btl = t - mstrip(0, 0);
    }
    else{
      // otherwise update the right bandwidth value
      btr = mstrip(0, 0) - t;
    }
    //bti = fabs(t - mstrip(0, 0));
    out[1] = na[0]; // set N-A estimator
  }
  for (int i = 1; i < j; i++){
    // include rest of qualified observations into hazard calculation
    // calculate remaining N-A estimates
    nac = (mstrip(i, 1)/(j - i)); // censoring indicat. divided by risk set size
    na[i] = na[i-1] + nac;
    if (fabs(tp - km[i]) <= bt){
      // if close enough to t, include in hazard calculation
      out[0] += nac;
      // if (fabs(t-mstrip(i, 0)) > bti){
      // also check how far on time-scale it is for eventual bandwidth choice
      //bti = fabs(t - mstrip(i, 0));
      if (((tp - km[i]) < 0) && ((t-mstrip(i, 0)) > btl)){
        btl = t-mstrip(i, 0);
      }
      else if (((tp - km[i]) >= 0) && ((mstrip(i, 0)-t) > btr)){
        btr = mstrip(i, 0)-t;
      }
      //      }
    }
    if (tp == km[i]){
      // if we have reached t, return this for N-A estimate
      out[1] = na[i];
    }
  }
  out[0] = out[0]/(btl+btr);
  //  if (btl > btr){out[0] = out[0]/(btl*2);}
  //  else{out[0] = out[0]/(btr*2);}
  if ((btl+btr) == 0){
    out[0] = 0;
  }
  out[2] = tp;
  out[3] = km[0];
  out[4] = btl;
  out[5] = btr;
  return(out);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#hazr(42)
*/
