#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
SEXP  residualcokurtosisSF(SEXP sstockM2, SEXP sstockM4, SEXP mfactorM2, SEXP bbeta  ){ 
  using namespace Rcpp ;
  NumericVector stockM2(sstockM2);
  NumericVector stockM4(sstockM4);
  double factorM2 = as<double>(mfactorM2);
  NumericVector beta(bbeta);
  int N = beta.size();
  
  NumericMatrix M4SF(N, pow(N,3));
  double kijkl;
  
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      for(int k=0; k<N; k++) {
        for(int l=0; l<N; l++) {
          // in the general case: i, j , k and l are all different
          kijkl = 0;
          // if at least one of them are equal: 
            if( (i==j) || (i==k) || (i==l) || (j==k) || (j==l) || (k==l) ) {
              if( (i==j) && (i==k) && (i==l) ) { 
                // These are the kurtosis estimates of the individual assets: E[u^4]
                kijkl = 6*pow(beta[i],2)*factorM2*stockM2[i]+stockM4[i]; 
              } else {
                if( ((i==j) && (i==k)) || ((i==j) && (i==l)) || ((i==k) && (i==l)) || ((j==k) && (j==l)) ) {
                  // kiij E[ U[,i]^3*U[,j] ] = r3*sqrt( vm6[i]*vm2[j] )
                  if( (i==j) && (i==k) ) {
                    kijkl = 3*beta[i]*beta[l]*factorM2*stockM2[i];
                  } else
                    if( (i==j) && (i==l) ) {
                      kijkl = 3*beta[i]*beta[k]*factorM2*stockM2[i];
                    } else
                      if( (i==k) && (i==l) ) {
                        kijkl = 3*beta[i]*beta[j]*factorM2*stockM2[i];
                      } else
                        if( (j==k) && (j==l) ) {
                          kijkl = 3*beta[j]*beta[i]*factorM2*stockM2[j];
                        }
                } else {
                  if( ((i==j) && (k==l)) || ((i==k) && (j==l)) || ((i==l) && (j==k)) ) { 
                    // kiijj = E[ U[,i]^2*U[,j]^2 ] = r5*sqrt( vm4[i]*vm4[j]  )
                    if( (i==j) && (k==l) ) {
                      kijkl = pow(beta[i],2)*factorM2*stockM2[k] + pow(beta[k],2)*factorM2*stockM2[i]+stockM2[i]*stockM2[k];
                    } else
                      if( (i==k) && (j==l) ) {
                        kijkl = pow(beta[i],2)*factorM2*stockM2[j] + pow(beta[j],2)*factorM2*stockM2[i]+stockM2[i]*stockM2[j];
                      } else
                        if( (i==l) && (j==k) ) {
                          kijkl = pow(beta[i],2)*factorM2*stockM2[j] + pow(beta[j],2)*factorM2*stockM2[i]+stockM2[i]*stockM2[j];
                        } 
                  } else {
                    // kiijk = E[ U[,i]^2*U[,j]*U[,k] ] = r6*sqrt( vm4[i]*r5*sqrt( vm4[j]*vm4[k] ) )
                    if( (i==j) ) {
                      kijkl = beta[k]*beta[l]*factorM2*stockM2[i];
                    } else
                      if( (i==k) ) {
                        kijkl = beta[j]*beta[l]*factorM2*stockM2[i];
                      } else
                        if( (i==l) ) {
                          kijkl = beta[j]*beta[k]*factorM2*stockM2[i];
                        } else
                          if( (j==k) ) {
                            kijkl = beta[i]*beta[l]*factorM2*stockM2[j];
                          } else
                            if( (j==l) ) {
                              kijkl = beta[i]*beta[k]*factorM2*stockM2[j];
                            } else
                              if( (k==l) ) {
                                kijkl = beta[i]*beta[j]*factorM2*stockM2[k];
                              }
                  }
                }
              }  // no kurtosis estimates of individual assets
            }  // at least one of them are not equal
          // i denotes the subcomoment matrix
          M4SF(k,i*pow(N,2)+j*N+l) = kijkl;
        } // loop l
      } // loop k
    } // loop j
  } // loop i
  return M4SF;
}
