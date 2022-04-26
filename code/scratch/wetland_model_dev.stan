functions {

   row_vector wetlandModel(
     real bPET,
     real bTreat,
     int T,
     real bRain,
     real bInt,
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
            // Create lagged rain
            Rain[t] = rain[t] + rain[t-1] * phiRain;
            wlHat[t] = wlHat[t - 1];
            // Esy
            // Calculate gradient of drawdown 
            gradient[t] = fmax(
               esymin,
               esyA - (esyA - esyB) * exp(esyC * wlHat[t])
            );
            // PET or P times Esy
            petHat[t] = (bPET * (1 - bTreat * (T - 1))) * pet[t] * gradient[t];
                  // pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
            wlHat[t] = wlHat[t] + petHat[t];
            pHat[t] = bRain * Rain[t] * gradient[t];
            wlHat[t] = wlHat[t] + pHat[t];
            
            // Snowmelt
            mHat[t] = (bRain * (1 + bInt * (T-1))) * melt[t] * gradient[t];
            wlHat[t] = wlHat[t] + mHat[t];

            // Streamflow
            // If WL is above spill point threshold then lose some to streamflow. 
            // This could probably be improved using the morphology models to determine
            // streamflow
            if(wlHat[t-1] > maxWL){
               qHat[t] = bQ * (wlHat[t-1] - maxWL);
               qHat[t] = fmin(qHat[t], wlHat[t] - maxWL);
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
   matrix[K,4] esyParams; // minESY, esyA, esyB, esyC
   matrix[K, T] maxWL;
}
parameters {
   real<lower = 0> bPET;
   real<lower = 0, upper = 1> bTreat;
   real<lower = 0, upper = 1> bInt;
   real<lower = 0> bRain;
   real<lower = 0> bQ;
   real<lower = 0> bphiRain;

   real<lower = 0> taubPET;
   real<lower = 0> taubRain;
   real<lower = 0> tauphiRain;

   row_vector<offset = bPET, multiplier = taubPET>[K] gPET;
   row_vector<offset = bRain, multiplier = taubRain>[K] gRain;
   row_vector<offset = bphiRain, multiplier = tauphiRain>[K] gphiRain;
   
   real<lower = 0> sigma;
}
transformed parameters{
   row_vector[K] gQ;
   row_vector[K] gTreat;
   row_vector[K] gInt;

   gQ = rep_row_vector(bQ, K); // repeating for easy coef extraction of pop and site level estimates
   gTreat = rep_row_vector(bTreat, K); // repeating for easy coef extraction of pop and site level estimates
   gInt = rep_row_vector(bInt, K); // repeating for easy coef extraction of pop and site level estimates
}
model {
   array[2] matrix[K,D] yHat;
   
   // Population Estimates
   bPET ~ normal(1.75, 0.25);
   bTreat ~ normal(0.5, 0.1);
   bRain ~ normal(2, 0.25);
   bInt ~ normal(0.2, 0.1);
   bphiRain ~ normal(0.2, 0.05);
   bQ ~ normal(0.5, 0.1);
   sigma ~ std_normal();

   // Standard Deviation of Group Effects around Population Effects
   taubPET ~ normal(0, 0.05);
   taubRain ~ normal(0, 0.05);
   tauphiRain ~ normal(0, 0.05);
   
   for(k in 1:K) {
      for(t in 1:T) {

         gPET[k] ~ normal(bPET, taubPET);
         gRain[k] ~ normal(bRain, taubRain);
         gphiRain[k] ~ normal(bphiRain, tauphiRain);
         
         yHat[t, k] = wetlandModel(
            gPET[k],
            gTreat[k],
            t,
            gRain[k],
            gInt[k],
            gphiRain[k],
            gQ[k],
            pet[t, k],
            rain[t, k],
            melt[t, k],
            D,
            maxWL[k, t],
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

