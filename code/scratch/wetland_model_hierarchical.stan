functions {
   #include wetland_model.stanfunctions
}
data {
   // Need to define 365 x K matrice for each wetlandModel input
   // all wetland model inputs
   // sites
   // weights
   // obs_sd for daily chage in WL
   int D; // number of days in year
   int K; // number of sites
   matrix[K,D] pet;
   matrix[K,D] rain;
   matrix[K,D] melt;
   matrix[K,D] wghts;
   matrix[K,D] y; // PET
   matrix[K,4] esyParams; // minESY, esyA, esyB, esyC
   vector[K] maxWL;
   vector[K] ySD;
}
parameters {
  real<lower=0, upper = 3> bPETPop;
  real<lower=0, upper = 3> bRainPop;
  real<lower=0, upper = 3> bMeltPop;
  real<lower=0, upper = 3> bQPop;
  real<lower=0> tauPETPop;
  real<lower=0> tauRainPop;
  real<lower=0> tauMeltPop;
  real<lower=0> tauQPop;
  vector<lower=0, upper = 3>[K] bPET;
  vector<lower=0, upper = 3>[K] bRain;
  vector<lower=0, upper = 3>[K] bMelt;
  vector<lower=0, upper = 3>[K] bQ;
  vector<lower=0>[K] sigma;
}
model {
   matrix[K,D] yHat;
   bPETPop ~ normal(1, 0.5);
   bRainPop ~ normal(1.5, 0.75);
   bMeltPop ~ normal(1, 0.5);
   bQPop ~ normal(0.5, 0.25);
   tauPETPop ~ gamma(1, 0.1);
   tauRainPop ~ gamma(1, 0.1);
   tauMeltPop ~ gamma(1, 0.1);
   tauQPop ~ gamma(1, 0.1);

   for(k in 1:K) {
      sigma[k] ~ normal(ySD, 2);
      bPET[k] ~ normal(bPETPop, tauPETPop);
      bRain[k] ~ normal(bRainPop, tauRainPop);
      bMelt[k] ~ normal(bMeltPop, tauMeltPop);
      bQ[k] ~ normal(bQPop, tauQPop);
      yHat[k] = wetlandModel(D, maxWL[k], y[k], melt[k], bMelt[k], pet[k], bPET[k], rain[k], bRain[k], bQ[k], esyParams[k,1], esyParams[k,2], esyParams[k,3], esyParams[k,4]);
      for(i in 1:D){
         target +=  normal_lpdf(y[k,i] | yHat[k,i], sigma[k]) * wghts[k,i];
      }
   }
}
// generated quantities{
//    // real bMelt = 1; // ~ normal(1, 0.5);
//    real bQ = 0.5; // ~ normal(0.5, 0.5);
//    array[K,D] real wl;
//    for(k in 1:K) {
//           wl[k] = normal_rng(wetlandModel(D, maxWL[k], y[k], melt[k], bMelt, pet[k], bPET, rain[k], bRain, bQ, esyParams[k,1], esyParams[k,2], esyParams[k,3], esyParams[k,4]), sigma[k]);
//    }
// }
