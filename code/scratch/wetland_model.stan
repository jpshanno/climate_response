functions {
   real predEsy(real wl, real esyA, real esyB, real esyC){
      return esyA - (esyA - esyB) * exp(esyC * wl);
   }

   row_vector wetlandModel(int D, real maxWL, row_vector obsWL, row_vector melt, real bMelt, row_vector pet, real bPET, row_vector rain, real bRain, real bQ, real minESY, real esyA, real esyB, real esyC){
      /*
      if(is.na(future.forest.change)){
            pet_fun <- 
            function(wl, max.wl, future.forest.change){ 1 }
      } else {
            pet_fun <- 
            function(wl, max.wl, future.forest.change){
            (1 - future.forest.change) + future.forest.change * pmin(1, abs(1/(1.45077 - 0.05869 * (wl - max.wl))))
            }
      }
      */
      
      // P <- as.numeric(filter(data$rain_cm, filter = phiP, method = "rec", sides = 1))
      // M <- as.numeric(filter(data$melt_cm, filter = phiM, method = "rec", sides = 1))
      // Create empty vectors
      row_vector[D] wlHat;
      vector[D] gradient;
      vector[D] qHat;
      vector[D] mHat;
      vector[D] pHat;
      vector[D] petHat;
      
      // Initialize model at full water level
      wlHat[1] = maxWL;
      
      // Loop through weather data
      for(t in 2:D){

            // Esy
            // Calculate gradient of drawdown 
            if(wlHat[t-1] > maxWL) {
            gradient[t] = minESY;
            } else {
            gradient[t] = predEsy(wlHat[t-1], esyA, esyB, esyC);
            }

            // PET or P times Esy
            // Use net input to determine if water level increases or decreases
            // Assuming AET is negligible on days where P >= PET
            if((rain[t] + pet[t]) <= 0){
            // Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
            // less than interception (not necessarily true, but works as a
            // simplifying assumption)
            petHat[t] = bPET * pet[t] * gradient[t];
                  // pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
            wlHat[t] = wlHat[t-1] + petHat[t];
            } else {
            pHat[t] = bRain * rain[t] * gradient[t];
            wlHat[t] = wlHat[t-1] + pHat[t];
            }

            // Snowmelt
            mHat[t] = bMelt * melt[t] * gradient[t];
            wlHat[t] = wlHat[t] + mHat[t];

            // Streamflow
            // If WL is above spill point threshold then lose some to streamflow. 
            // This could probably be improved using the morphology models to determine
            // streamflow
            if(wlHat[t-1] > maxWL){
            qHat[t] = bQ * (wlHat[t-1] - maxWL);
            if(qHat[t] > wlHat[t] - maxWL){
                  qHat[t] = wlHat[t] - maxWL;
            }
            wlHat[t] = wlHat[t] - qHat[t];
            }

      }
      return wlHat;
   }
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
   vector[K] obs_sigma; // mean of daily water level change by site
}
parameters {
   real<lower = 0> bPET;
   real<lower = 0> bRain;
   real<lower = 0> bMelt;
   real<lower = 0> bQ;
   real<lower = 0> sigma;
   // vector<lower = 0, upper = 2>[4] params;
   // cholesky_factor_corr[4] Sigma;
}
model {
   matrix[K,D] yHat;
   // bPET ~ normal(1, 0.5);
   target += normal_lpdf(bPET | 1, 0.25);
   // bRain ~ normal(1.5, 0.75);
   target += normal_lpdf(bRain | 1.2, 0.25);
   // bMelt ~ normal(1, 0.5);
   target += normal_lpdf(bMelt | 1.2, 0.25);
   // bQ ~ normal(0.5, 0.25);
   target += normal_lpdf(bQ | 0.6, 0.1);
   target += normal_lpdf(sigma | obs_sigma, obs_sigma / 4);
   // sigma ~ std_normal();

   // Sigma ~ lkj_corr_cholesky(0.5);
   // params ~ multi_normal_cholesky([1, 1.5, 1, 0.5], Sigma);
   // real bPET;
   // real bRain;
   // real bMelt;
   // real bQ;
   // bPET = params[1];
   // bRain = params[2];
   // bMelt = params[3];
   // bQ = params[4];

   for(k in 1:K) {
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
