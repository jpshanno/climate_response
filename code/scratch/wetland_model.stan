functions {

   row_vector wetlandModel(
     real bPET,
     real bRain,
     real bMelt,
     real phiRain,
     real bQ,
     row_vector pet,
     row_vector rain,
     row_vector melt,
     int D,
     real maxWL,
     real esyA,
     real esyB,
     real esyC,
     real esymin){
            
      // Create empty vectors
      row_vector[D] wlHat;
      vector[D] Rain;
      // vector[D] Melt;
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
            // Melt[t] = melt[t] + melt[t-1] * phiMelt;
            wlHat[t] = wlHat[t - 1];
            // Esy
            // Calculate gradient of drawdown 
            // gradient[t] = fmax(esyint - esyslope * wlHat[t], esymin);
            gradient[t] = fmax(
               esymin,
               esyA - (esyA - esyB) * exp(esyC * wlHat[t])
            );
            // gradient[t] = fmin(gradient[t], esymax);
            // PET or P times Esy
            // Use net input to determine if water level increases or decreases
            // Assuming AET is negligible on days where P >= PET
            // Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
            // less than interception (not necessarily true, but works as a
            // simplifying assumption)
            petHat[t] = bPET * pet[t] * gradient[t];
                  // pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
            wlHat[t] = wlHat[t] + petHat[t];
            pHat[t] = bRain * Rain[t] * gradient[t];
            wlHat[t] = wlHat[t] + pHat[t];
            

            // Snowmelt
            mHat[t] = bMelt * melt[t];
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
   // vector[K] obs_sigma; // mean of daily water level change by site
}
parameters {
   real<lower = 0> bPET;
   real<lower = 0> bRain;
   real<lower = 0> bMelt;
   real<lower = 0> bQ;
   real<lower = 0> bphiRain;
   // real<lower = 0> bphiMelt;
   // real<lower = 0> bEsyInt;
   // real<lower = 0> bEsySlope;
   // real<lower = 0> bEsyMin;
   // real<lower = 0> taubPET;
   // real<lower = 0> taubRain;
   // real<lower = 0> taubMelt;
   // real<lower = 0> taubQ;
   // real<lower = 0> tauphiRain;
   // real<lower = 0> tauphiMelt;
   // real<lower = 0> taubEsyInt;
   // real<lower = 0> taubEsySlope;
   // real<lower = 0> taubEsyMin;
   // row_vector<offset = bPET, multiplier = taubPET>[K] gPET;
   // row_vector<offset = bRain, multiplier = taubRain>[K] gRain;
   // row_vector<offset = bMelt, multiplier = taubMelt>[K] gMelt;
   // row_vector<offset = bQ, multiplier = taubQ>[K] gQ;
   // row_vector<offset = bphiRain, multiplier = tauphiRain>[K] gphiRain;
   // row_vector<offset = bphiMelt, multiplier = tauphiMelt>[K] gphiMelt;
   // row_vector<offset = bEsyInt, multiplier = taubEsyInt>[K] gEsyInt;
   // row_vector<offset = bEsySlope, multiplier = taubEsySlope>[K] gEsySlope;
   // row_vector<offset = bEsyMin, multiplier = taubEsyMin>[K] gEsyMin;
   // real<lower = 0> omega;
   // real<upper = 0> alpha;
   // real<lower = 0> bEsyMax;
   real<lower = 0> sigma;
}
model {
   matrix[K,D] yHat;
   
   // Population Estimates
   // target += gamma_lpdf(bPET | 10, 10);
   // target += gamma_lpdf(bRain | 11, 10);
   // target += gamma_lpdf(bMelt | 11, 10);
   // target += gamma_lpdf(bQ | 3, 4);
   // target += exponential_lpdf(bphiRain | 10);
   // target += exponential_lpdf(bphiMelt | 10);
   bPET ~ normal(1, 0.25);
   bRain ~ normal(1, 0.25);
   bMelt ~ normal(8, 1);
   bphiRain ~ normal(0.2, 0.05);
   bQ ~ normal(0.2, 0.05);
   // bEsyInt ~ normal(2, 1);
   // bEsySlope ~ normal(0.02, 0.005);
   // bEsyMin ~ normal(0.75, 0.2);
   // target += std_normal_lpdf(omega);
   // target += std_normal_lpdf(alpha);
   sigma ~ std_normal();
   // bEsyMax ~ gamma(20, 10);

   // Standard Deviation of Group Effects around Population Effects
   // target += normal_lpdf(taubPET | 0, 0.05);
   // target += normal_lpdf(taubRain | 0, 0.05);
   // target += normal_lpdf(taubMelt | 0, 0.05);
   // target += normal_lpdf(taubQ | 0, 0.05);
   // target += normal_lpdf(tauphiRain | 0, 0.05);
   // target += normal_lpdf(tauphiMelt | 0, 0.05);
   // target += std_normal_lpdf(taubEsyInt);
   // target += normal_lpdf(taubEsySlope | 0, 0.005);
   // target += normal_lpdf(taubEsyMin | 0, 0.05);
   
   for(k in 1:K) {
      // target += normal_lpdf(gPET[k] | bPET, taubPET);
      // target += normal_lpdf(gRain[k] | bRain, taubRain);
      // target += normal_lpdf(gMelt[k] | bMelt, taubMelt);
      // target += normal_lpdf(gQ[k] | bQ, taubQ);
      // target += normal_lpdf(gphiRain[k] | bphiRain, tauphiRain);
      // target += normal_lpdf(gphiMelt[k] | bphiMelt, tauphiMelt);
      // target += normal_lpdf(gEsyInt[k] | bEsyInt, taubEsyInt);
      // target += normal_lpdf(gEsySlope[k] | bEsySlope, taubEsySlope);
      // target += normal_lpdf(gEsyMin[k] | bEsyMin, taubEsyMin);
      // vector[N[k]] esyHat;
      // esyHat = fitEsy(waterBalance, y, gasympA[k], gasympB[k], gasympC[k]);
      // target += normal_lpdf()
      // maxEsy = bEsyInt - bEsySlope * min(y[k, ]);
      yHat[k] = wetlandModel(
      //   bParams[1],
      //   bParams[2],
        bPET,
        bRain,
        bMelt,
        bphiRain,
        bQ,
        pet[k],
        rain[k],
        melt[k],
        D,
        maxWL[k],
        esyParams[k, 2],
        esyParams[k, 3],
        esyParams[k, 4],
        esyParams[k, 1]
        );
      for(i in 1:D){
         target +=  normal_lpdf(y[k,i] | yHat[k,i], sigma) * wghts[k,i];
      }
   }
}
// generated quantities{
//    matrix[K,D] fitted;
//    for(k in 1:K) {
//       // target += normal_lpdf(gPET[k] | bPET, taubPET);
//       // target += normal_lpdf(gRain[k] | bRain, taubRain);
//       // target += normal_lpdf(gMelt[k] | bMelt, taubMelt);
//       // target += normal_lpdf(gQ[k] | bQ, taubQ);
//       // target += normal_lpdf(gphiRain[k] | bphiRain, tauphiRain);
//       // target += normal_lpdf(gphiMelt[k] | bphiMelt, tauphiMelt);
//       // target += normal_lpdf(gEsyInt[k] | bEsyInt, taubEsyInt);
//       // target += normal_lpdf(gEsySlope[k] | bEsySlope, taubEsySlope);
//       // target += normal_lpdf(gEsyMin[k] | bEsyMin, taubEsyMin);
//       // vector[N[k]] esyHat;
//       // esyHat = fitEsy(waterBalance, y, gasympA[k], gasympB[k], gasympC[k]);
//       // target += normal_lpdf()
//       fitted[k] = to_row_vector(normal_rng(
//          wetlandModel(
//             pet[k],
//             rain[k],
//             melt[k],
//             D,
//             maxWL[k],
//             bEsyInt,
//             bEsySlope,
//             bEsyMin
//             ),
//          sigma
//       ));
//    }
// }
