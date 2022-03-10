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
            // Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
            // less than interception (not necessarily true, but works as a
            // simplifying assumption)
            petHat[t] = bPET * pet[t] * gradient[t];
                  // pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
            wlHat[t] = wlHat[t-1] + petHat[t];
            pHat[t] = bRain * rain[t] * gradient[t];
            wlHat[t] = wlHat[t-1] + pHat[t];
            

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
   vector<lower = 0>[4] bPop;
   matrix<lower = 0>[K, 4] bGroup;
   real<lower = 0> sigma;
   vector<lower = 0>[4] tau;
}
model {
   matrix[K,D] yHat;

   // Population Estimates
   target += lognormal_lpdf(bPop[1] | 1, 1);
   target += lognormal_lpdf(bPop[2] | 2, 1);
   target += lognormal_lpdf(bPop[3] | 1, 1);
   target += lognormal_lpdf(bPop[4] | 0.6, 0.25);
   target += std_normal_lpdf(sigma);

   // Group Effects
   for(p in 1:4) {
      target += normal_lpdf(tau[p] | 0, 0.5);
   }

   for(k in 1:K) {
      for(p in 1:4) {
         target += lognormal_lpdf(bGroup[k, p] | bPop[p], tau[p]);
      }
      yHat[k] = wetlandModel(D, maxWL[k], y[k], melt[k], bGroup[k, 1], pet[k], bGroup[k, 2], rain[k], bGroup[k, 3], bGroup[k, 4], esyParams[k,1], esyParams[k,2], esyParams[k,3], esyParams[k,4]);
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
