// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP coskewkurtosis_rcpp_hello_world() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = rcpp_hello_world();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// residualcokurtosisMF
SEXP residualcokurtosisMF(SEXP NN, SEXP sstockM2, SEXP sstockM4, SEXP bbetacov);
RcppExport SEXP coskewkurtosis_residualcokurtosisMF(SEXP NNSEXP, SEXP sstockM2SEXP, SEXP sstockM4SEXP, SEXP bbetacovSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type NN(NNSEXP );
        Rcpp::traits::input_parameter< SEXP >::type sstockM2(sstockM2SEXP );
        Rcpp::traits::input_parameter< SEXP >::type sstockM4(sstockM4SEXP );
        Rcpp::traits::input_parameter< SEXP >::type bbetacov(bbetacovSEXP );
        SEXP __result = residualcokurtosisMF(NN, sstockM2, sstockM4, bbetacov);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// residualcokurtosisSF
SEXP residualcokurtosisSF(SEXP sstockM2, SEXP sstockM4, SEXP mfactorM2, SEXP bbeta);
RcppExport SEXP coskewkurtosis_residualcokurtosisSF(SEXP sstockM2SEXP, SEXP sstockM4SEXP, SEXP mfactorM2SEXP, SEXP bbetaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type sstockM2(sstockM2SEXP );
        Rcpp::traits::input_parameter< SEXP >::type sstockM4(sstockM4SEXP );
        Rcpp::traits::input_parameter< SEXP >::type mfactorM2(mfactorM2SEXP );
        Rcpp::traits::input_parameter< SEXP >::type bbeta(bbetaSEXP );
        SEXP __result = residualcokurtosisSF(sstockM2, sstockM4, mfactorM2, bbeta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
