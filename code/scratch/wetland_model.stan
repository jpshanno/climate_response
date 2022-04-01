functions {

   // row_vector fitEsy(
   //    row_vector waterBalance,
   //    row_vector waterLevel,
   //    real asympA,
   //    real asympB,
   //    real asympC
   // ){
   //    row_vector drawdown;
   //    row_vector esyHat;
   //    row_vector esyCoefs;
   //    drawdown = a - (a - b) * exp(-c * waterBalance);
   //    esyHat = (a - b) * (exp(-c * drawdown) * c);
   //    return esyHat;
   // }

   row_vector wetlandModel(
     real bPET,
     real bRain,
     real bMelt,
     real bQ,
     real phiRain,
     real phiMelt,
     row_vector pet,
     row_vector rain,
     row_vector melt,
     int D,
     real maxWL,
     real esyint,
     real esyslope,
     real minObsWL){
            
      // Create empty vectors
      row_vector[D] wlHat;
      vector[D] Rain;
      vector[D] Melt;
      vector[D] gradient;
      vector[D] qHat;
      vector[D] mHat;
      vector[D] pHat;
      vector[D] petHat;
      
      // Initialize model at full water level
      wlHat[1] = maxWL;
      
      // Loop through weather data
      for(t in 2:D){
            Rain[t] = rain[t] + rain[t-1] * phiRain;
            Melt[t] = melt[t] + melt[t-1] * phiMelt;
            wlHat[t] = wlHat[t - 1];
            // Esy
            // Calculate gradient of drawdown 
            gradient[t] = fmax(esyint - esyslope * wlHat[t], 1);
            
            // PET or P times Esy
            // Use net input to determine if water level increases or decreases
            // Assuming AET is negligible on days where P >= PET
            // Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
            // less than interception (not necessarily true, but works as a
            // simplifying assumption)
            petHat[t] = bPET * pet[t] * gradient[t];
                  // pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
            wlHat[t] = wlHat[t] + petHat[t] * (pet[t] > rain[t]);
            pHat[t] = bRain * Rain[t] * gradient[t];
            wlHat[t] = wlHat[t] + pHat[t] * (pet[t] <= rain[t]);
            

            // Snowmelt
            mHat[t] = bMelt * Melt[t] * gradient[t];
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
   // vector[K] nEsy; // number of days used in fitting Esy
   matrix[K,D] pet;
   matrix[K,D] rain;
   matrix[K,D] melt;
   matrix[K,D] wghts;
   matrix[K,D] y; // PET
   matrix[K,6] esyParams; // minESY, esyA, esyB, esyC
   vector[K] maxWL;
   vector[K] obs_sigma; // mean of daily water level change by site
}
parameters {
   real<lower = 0> bPET;
   real<lower = 0> bRain;
   real<lower = 0> bMelt;
   real<lower = 0> bQ;
   real<lower = 0> bphiRain;
   real<lower = 0> bphiMelt;
   real<lower = 0> bEsyInt;
   real<lower = 0> bEsySlope;
   real<lower = 0> taubPET;
   real<lower = 0> taubRain;
   real<lower = 0> taubMelt;
   real<lower = 0> taubQ;
   real<lower = 0> tauphiRain;
   real<lower = 0> tauphiMelt;
   real<lower = 0> tauEsyInt;
   real<lower = 0> tauEsySlope;
   row_vector<offset = bPET, multiplier = taubPET>[K] gPET;
   row_vector<offset = bRain, multiplier = taubRain>[K] gRain;
   row_vector<offset = bMelt, multiplier = taubMelt>[K] gMelt;
   row_vector<offset = bQ, multiplier = taubQ>[K] gQ;
   row_vector<offset = bphiRain, multiplier = tauphiRain>[K] gphiRain;
   row_vector<offset = bphiMelt, multiplier = tauphiMelt>[K] gphiMelt;
   row_vector<offset = bEsyInt, multiplier = tauEsyInt>[K] gEsyInt;
   row_vector<offset = bEsySlope, multiplier = tauEsySlope>[K] gEsySlope;
   real<lower = 0> sigma;
   real<lower = 0> sigmaEsy;
}
model {
   matrix[K,D] yHat;
   real minObsWL;
   
   // Population Estimates
   target += gamma_lpdf(bPET | 10, 10);
   target += gamma_lpdf(bRain | 11, 10);
   target += gamma_lpdf(bMelt | 11, 10);
   target += gamma_lpdf(bQ | 3, 4);
   target += exponential_lpdf(bphiRain | 10);
   target += exponential_lpdf(bphiMelt | 10);
   target += normal_lpdf(bEsyInt | 2, 1);
   target += gamma_lpdf(bEsySlope | 0.5, 10);
   target += std_normal_lpdf(sigma);
   target += std_normal_lpdf(sigmaEsy);

   // Standard Deviation of Group Effects around Population Effects
   target += normal_lpdf(taubPET | 0, 0.05);
   target += normal_lpdf(taubRain | 0, 0.05);
   target += normal_lpdf(taubMelt | 0, 0.05);
   target += normal_lpdf(taubQ | 0, 0.05);
   target += normal_lpdf(tauphiRain | 0, 0.05);
   target += normal_lpdf(tauphiMelt | 0, 0.05);
   target += std_normal_lpdf(tauEsyInt);
   target += normal_lpdf(tauEsySlope | 0, 0.05);
   
   for(k in 1:K) {
      target += normal_lpdf(gPET[k] | bPET, taubPET);
      target += normal_lpdf(gRain[k] | bRain, taubRain);
      target += normal_lpdf(gMelt[k] | bMelt, taubMelt);
      target += normal_lpdf(gQ[k] | bQ, taubQ);
      target += normal_lpdf(gphiRain[k] | bphiRain, tauphiRain);
      target += normal_lpdf(gphiMelt[k] | bphiMelt, tauphiMelt);
      target += normal_lpdf(gEsyInt[k] | bEsyInt, tauEsyInt);
      target += normal_lpdf(gEsySlope[k] | bEsySlope, tauEsySlope);
      // vector[N[k]] esyHat;
      // esyHat = fitEsy(waterBalance, y, gasympA[k], gasympB[k], gasympC[k]);
      // target += normal_lpdf()
      minObsWL = max(y[k, ]);
      yHat[k] = wetlandModel(
        gPET[k],
        gRain[k],
        gMelt[k],
        gQ[k],
        gphiRain[k],
        gphiMelt[k],
        pet[k],
        rain[k],
        melt[k],
        D,
        maxWL[k],
        gEsyInt[k],
        gEsySlope[k],
        minObsWL);
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
