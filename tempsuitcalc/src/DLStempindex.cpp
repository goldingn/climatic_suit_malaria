#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector DLStempindex (NumericVector myDD, NumericVector myovi, NumericMatrix myP, NumericVector Zout, NumericVector Zoutovi)
{
  const int L = 10416 / 2;
  const int L2 = 792 / 2;
  const int MnB = 47 / 2;
  const double Kappa = 0.01;
  double Y[L2];
  double s[L2];
  double t[L2];
    
  for(int i=0; i < (L - L2); i++){

    // start new cohort
    double M = 100.0;
   
    for(int v=0; v < L2; v++){
      // set s and Y to 0
      s[v] = 0.0;
      t[v] = 0.0;
      Y[v] = 0.0;
    }
    
    for(int j=0; j < MnB; j++){
      // apply normal mortality whilst not biting
      M = M * myP(j,i);
    }
        
    for(int k = MnB; k < L2; k++){
      // loop through biting age mosquitoes

      // infect mosquitoes
      Y[k] = Kappa * M;
      
      for(int r=MnB; r <= k; r++){
        // evaluate degree-day accumulation
	      s[r] += myDD[(i + r)];
        t[r] += myovi[(i + r)];
        // apply mortality to infected mosquitoes
	      Y[r] = Y[r] * myP(k, (i + r));
        // if they have reached virus incubation
	      if(s[r] >= 1){
          // add them to infectious population
	        Zout[(i + r)] += Y[r];
        }
        // if they have reached oviposition
        if(t[r] >= 1){
          // add them to parous population
	        Zoutovi[(i + r)] += Y[r];
        }
    	}

      // apply mortality to uninfected mosquitoes
      M = M * myP(k,i);
    }
	}
  return(1);
}
