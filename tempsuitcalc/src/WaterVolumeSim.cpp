#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector WaterVolumeSim (NumericVector temp, NumericVector dewp, NumericVector rain, double alt, double lat, NumericVector Vout)
{
  const int L = 10416 / 2;
  const double A = 1587.401;
  double rainfall = 0;
  double evapnum = 0;
  
  Vout[0] = 5947.08;   
  // This is equal to half of maximal volume, where maximal volume is A^(3/2)/(3*sqrt(pi))
  // At this volume, the area is approximately 1000
  
  for(int i=0; i < L-1; i++){
    
    rainfall = rain[i] * A;
    
    evapnum = 15*(temp[i]-dewp[i]) + 700*((temp[i]+0.006*alt)/(100-lat));
    
    Vout[i+1] = std::min(std::max(Vout[i] + rainfall - evapnum/(960 - 12*temp[i])*(pow(5.317362*Vout[i],2.0/3.0)) , 0.0) , 11894.16);
  }
  
  return(1);
}
