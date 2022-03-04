functions {
   #include wetland_model.stanfunctions
}
data {
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
   // vector<lower = 0>[4] params;
   real<lower=0, upper=2> bPET;
   real<lower=0, upper=2> bRain;
   real<lower=0, upper=2> bMelt;
   real<lower=0, upper=2> bQ;
   real<lower=0> sigma;
}
model {
   matrix[K,D] yHat;
   // real bPET;
   // real bRain;
   // real bMelt;
   // real bQ;
   // row_vector[4] priors;
   // priors = [1, 1.5, 1, 0.5];
   // Sigma ~ lkj_corr(2);
   // params ~ multi_normal(priors, Sigma);
   // bPET = params[1];
   // bRain = params[2];
   // bMelt = params[3];
   // bQ = params[4];
   bPET ~ normal(1, 0.5);
   bRain ~ normal(1.5, 0.75);
   bMelt ~ normal(1, 0.5);
   bQ ~ normal(0.5, 0.25);
   sigma ~ normal(4, 0.05);

   for(k in 1:K) {
      // sigma[k] ~ normal(ySD[k], 2);
      yHat[k] = wetlandModel(D, maxWL[k], y[k], melt[k], bMelt, pet[k], bPET, rain[k], bRain, bQ, esyParams[k,1], esyParams[k,2], esyParams[k,3], esyParams[k,4]);
      for(i in 1:D){
         target +=  normal_lpdf(y[k,i] | yHat[k,i], sigma) * wghts[k,i];
      }
   }
}
// generated quantities{
//    array[K,D] real wl;
//    for(k in 1:K) {
//           wl[k] = normal_rng(wetlandModel(D, maxWL[k], y[k], melt[k], bMelt, pet[k], bPET, rain[k], bRain, bQ, esyParams[k,1], esyParams[k,2], esyParams[k,3], esyParams[k,4]), sigma[k]);
//    }
// }
