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
   vector[K] ySD; // SD of daily water level change by site
}
parameters {
   real bPET;
   real bRain;
   real bMelt;
   real bQ;
   vector<lower=0>[K] sigma;
}
model {
   matrix[K,D] yHat;
   bPET ~ normal(1, 0.5);
   bRain ~ normal(1.5, 0.75);
   bMelt ~ normal(1, 0.5);
   bQ ~ normal(0.5, 0.25);

   for(k in 1:K) {
      sigma[k] ~ normal(ySD[k], 1);
      yHat[k] = wetlandModel(D, maxWL[k], y[k], melt[k], bMelt, pet[k], bPET, rain[k], bRain, bQ, esyParams[k,1], esyParams[k,2], esyParams[k,3], esyParams[k,4]);
      for(i in 1:D){
         target +=  normal_lpdf(y[k,i] | yHat[k,i], sigma[k]) * wghts[k,i];
      }
   }
}
// generated quantities{
//    array[K,D] real wl;
//    for(k in 1:K) {
//           wl[k] = normal_rng(wetlandModel(D, maxWL[k], y[k], melt[k], bMelt, pet[k], bPET, rain[k], bRain, bQ, esyParams[k,1], esyParams[k,2], esyParams[k,3], esyParams[k,4]), sigma[k]);
//    }
// }
