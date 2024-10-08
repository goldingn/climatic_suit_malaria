// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// DLStempindex
NumericVector DLStempindex(NumericVector myDD, NumericVector myovi, NumericMatrix myP, NumericVector Zout, NumericVector Zoutovi);
RcppExport SEXP _tempsuitcalc_DLStempindex(SEXP myDDSEXP, SEXP myoviSEXP, SEXP myPSEXP, SEXP ZoutSEXP, SEXP ZoutoviSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type myDD(myDDSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type myovi(myoviSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type myP(myPSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Zout(ZoutSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Zoutovi(ZoutoviSEXP);
    rcpp_result_gen = Rcpp::wrap(DLStempindex(myDD, myovi, myP, Zout, Zoutovi));
    return rcpp_result_gen;
END_RCPP
}
// DLStempindexnew
NumericVector DLStempindexnew(NumericVector myDD, NumericMatrix myP, NumericVector Zout, NumericVector Aoutmov);
RcppExport SEXP _tempsuitcalc_DLStempindexnew(SEXP myDDSEXP, SEXP myPSEXP, SEXP ZoutSEXP, SEXP AoutmovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type myDD(myDDSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type myP(myPSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Zout(ZoutSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Aoutmov(AoutmovSEXP);
    rcpp_result_gen = Rcpp::wrap(DLStempindexnew(myDD, myP, Zout, Aoutmov));
    return rcpp_result_gen;
END_RCPP
}
// WaterVolumeSim
NumericVector WaterVolumeSim(NumericVector temp, NumericVector dewp, NumericVector rain, double alt, double lat, NumericVector Vout);
RcppExport SEXP _tempsuitcalc_WaterVolumeSim(SEXP tempSEXP, SEXP dewpSEXP, SEXP rainSEXP, SEXP altSEXP, SEXP latSEXP, SEXP VoutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dewp(dewpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rain(rainSEXP);
    Rcpp::traits::input_parameter< double >::type alt(altSEXP);
    Rcpp::traits::input_parameter< double >::type lat(latSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Vout(VoutSEXP);
    rcpp_result_gen = Rcpp::wrap(WaterVolumeSim(temp, dewp, rain, alt, lat, Vout));
    return rcpp_result_gen;
END_RCPP
}
// WaterVolumeSimWithEvap
NumericVector WaterVolumeSimWithEvap(NumericVector temp, NumericVector rain, NumericVector evap, NumericVector Vout);
RcppExport SEXP _tempsuitcalc_WaterVolumeSimWithEvap(SEXP tempSEXP, SEXP rainSEXP, SEXP evapSEXP, SEXP VoutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rain(rainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type evap(evapSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Vout(VoutSEXP);
    rcpp_result_gen = Rcpp::wrap(WaterVolumeSimWithEvap(temp, rain, evap, Vout));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _tempsuitcalc_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tempsuitcalc_DLStempindex", (DL_FUNC) &_tempsuitcalc_DLStempindex, 5},
    {"_tempsuitcalc_DLStempindexnew", (DL_FUNC) &_tempsuitcalc_DLStempindexnew, 4},
    {"_tempsuitcalc_WaterVolumeSim", (DL_FUNC) &_tempsuitcalc_WaterVolumeSim, 6},
    {"_tempsuitcalc_WaterVolumeSimWithEvap", (DL_FUNC) &_tempsuitcalc_WaterVolumeSimWithEvap, 4},
    {"_tempsuitcalc_rcpp_hello_world", (DL_FUNC) &_tempsuitcalc_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_tempsuitcalc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
