functions {

   row_vector wetlandModel(
     real bPET,
     real bPETTreat,
     int T,
     real bRain,
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
            wlHat[t] = wlHat[t - 1];
            
            // Esy
            // Calculate gradient of drawdown 
            gradient[t] = fmax(
               esymin,
               esyA - (esyA - esyB) * exp(esyC * wlHat[t])
            );

            // PET or P times Esy
            // Use net input to determine if water level increases or decreases
            // Assuming AET is negligible on days where P >= PET
            // Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
            // less than interception (not necessarily true, but works as a
            // simplifying assumption)
            petHat[t] = (bPET * (T == 1) + bPETTreat * (T == 2)) * pet[t] * gradient[t];
                  // pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
            wlHat[t] = wlHat[t] + petHat[t];
            pHat[t] = bRain * Rain[t] * gradient[t];
            wlHat[t] = wlHat[t] + pHat[t];
            

            // Snowmelt
            mHat[t] = bRain * melt[t] * gradient[t];
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
   int T; // number of site statuses (2 for control & treatment)
   array[2] matrix[K,D] pet;
   array[2] matrix[K,D] rain;
   array[2] matrix[K,D] melt;
   array[2] matrix[K,D] wghts;
   array[2] matrix[K,D] y; // PET
   matrix[K,6] esyParams; // minESY, esyA, esyB, esyC
   vector[K] maxWL;
}
parameters {
   real<lower = 0> bPET;
   real<lower = 0> bRain;
   real<lower = 0> bQ;
   real<lower = 0> bphiRain;
   real<lower = 0> taubPET;
   // real<lower = 0> taubPETControl;
   // real<lower = 0> taubPETTreated;
   real<lower = 0> taubRain;
   real<lower = 0> tauphiRain;
   row_vector<lower = 0>[T] taubPETStatus;
   // real<lower = 0> bPETControl;
   // real<lower = 0> bPETTreated;
   row_vector<offset = bRain, multiplier = taubRain>[K] gRain;
   row_vector<offset = bPET, multiplier = taubPET>[T] bPETStatus;
   // row_vector<offset = bPET, multiplier = taubPET>[K] gPET;
   row_vector<offset = bPETStatus[1], multiplier = taubPETStatus[1]>[K] gPETControl;
   row_vector<offset = bPETStatus[2], multiplier = taubPETStatus[2]>[K] gPETTreated;
   row_vector<offset = bphiRain, multiplier = tauphiRain>[K] gphiRain;
   real<lower = 0> sigma;
}
// TODO: Test using real gQ in model
transformed parameters{
   row_vector[K] gQ;
   gQ = rep_row_vector(bQ, K);
}
model {
   array[2] matrix[K,D] yHat;
   
   // Population Estimates
   bPET ~ normal(1, 0.25);
   bRain ~ normal(1, 0.25);
   bphiRain ~ normal(0.2, 0.05);
   bQ ~ normal(0.2, 0.05);
   sigma ~ std_normal();

   // Standard Deviation of Group Effects around Population Effects
   taubPET ~ normal(0, 0.05);
   // taubPETControl ~ normal(0, 0.05);
   // taubPETTreated ~ normal(0, 0.05);
   taubPETStatus ~ normal(0, 0.05);
   taubRain ~ normal(0, 0.05);
   tauphiRain ~ normal(0, 0.05);

   // Treatment PET
   bPETStatus ~ normal(bPET, taubPET);
   // bPETControl ~ normal(bPET, tauPET);
   // bPETTreated ~ normal(bPET, tauPET);

   for(k in 1:K) {
      for(t in 1:T) {
         // Sample Group Level Parameters
         gPETControl[k] ~ normal(bPETStatus[1], taubPETStatus[1]);
         gPETTreated[k] ~ normal(bPETStatus[2], taubPETStatus[2]);
         gRain[k] ~ normal(bRain, taubRain);
         gphiRain[k] ~ normal(bphiRain, tauphiRain);
         
         yHat[t, k] = wetlandModel(
            gPETControl[k],
            gPETTreated[k],
            t,
            gRain[k],
            gphiRain[k],
            gQ[k],
            pet[t, k],
            rain[t, k],
            melt[t, k],
            D,
            maxWL[k],
            esyParams[k, 2],
            esyParams[k, 3],
            esyParams[k, 4],
            esyParams[k, 1]
            );
         for(i in 1:D){
            target +=  normal_lpdf(y[t,k,i] | yHat[t,k,i], sigma) * wghts[t,k,i];
         }
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
