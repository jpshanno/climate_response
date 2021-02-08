library(targets)
source("code/load_project.R")

# see time_series_models.R for prewhitening and ccf. This or the approach from 
# Hyndman should probably be done for each site to impove models and highlight
# potential differences based on site morphology or site status

tr_dat <- 
  tar_read(training_data)

tr_dat[, Ds_raw := Ds_cm / sy]

tr_dat[, `:=`(lag1_wl_cm = shift(wl_initial_cm, 1),
              lag1_precip_cm = shift(best_precip_cm, 1),
              lag2_precip_cm = shift(best_precip_cm, 2),
              lag3_precip_cm = shift(best_precip_cm, 3),
              lag1_pet_cm = shift(pet_cm, 1),
              lag2_pet_cm = shift(pet_cm, 2),
              lag3_pet_cm = shift(pet_cm, 3)),
       by = .(site, sample_year)]

# tr_dat[, 
#        .(pet_mod = list(forecast::auto.arima(pet_cm, max.q = 0, allowdrift = FALSE, allowmean = FALSE)),
#          precip_mod = list(forecast::auto.arima(best_precip_cm, max.q = 0, allowdrift = FALSE, allowmean = FALSE)),
#          wl_ts = list(wl_initial_cm),
#          pet_ts = list(pet_cm),
#          precip_ts = list(best_precip_cm)),
#        keyby = .(site, sample_year)
# ][, `:=`(pet_coefs = map(pet_mod, coef), 
#          precip_coefs = map(precip_mod, coef))
# ][, prewhite_wl := pmap(list(wl_ts, pet_coefs, precip_coefs),
#                         function(WL, PET, P){
#                           filter(
#                             x = filter(x = WL,
#                                        filter = PET,
#                                        method = "recursive",
#                                        sides = 1),
#                             filter = P,
#                             method = "recursive",
#                             sides = 1)
#                         }
# )]

ds_mmods <-
  tr_dat[, .(mod = list(lmrob(Ds_raw ~ 
                               best_precip_cm + 
                               lag1_precip_cm + 
                               lag2_precip_cm + 
                               lag3_precip_cm + 
                               pet_cm +
                               lag1_pet_cm + 
                               lag2_pet_cm + 
                               lag3_pet_cm,
                             data = .SD,
                             setting = "KS2014"))),
               keyby = .(site, site_status)]

# wl_mmods <-
#   tr_dat[, .(mod = list(lmer(wl_initial_cm ~ 0 + 
#                                lag1_wl + 
#                                best_precip_cm + 
#                                lag1_precip_cm + 
#                                lag2_precip_cm + 
#                                pet_cm +
#                                (0 + 
#                                   lag1_wl + 
#                                   best_precip_cm + 
#                                   lag1_precip_cm + 
#                                   lag2_precip_cm + 
#                                   pet_cm || site),
#                         data = .SD))),
#          keyby = .(site_status)]

# wl_coefs <- 
#   CJ(site = unique(tr_dat$site),
#      site_status = unique(tr_dat$site_status))
# 
# wl_coefs[, pop_coefs := list(list(lme4::fixef(ds_mmods[.BY[[1]], mod[[1]]]))),
#          by = .(site_status)]
# 
# wl_coefs[, group_coefs := list(list(as.numeric((lme4::ranef(ds_mmods[.BY[[1]], mod[[1]]])$site[.BY[[2]],])))),
#          by = .(site_status, site)]
# 
# wl_coefs[, group_coefs := map(group_coefs,
#                               ~ifelse(is.na(first(.x)),
#                                       rep(0, length(wl_coefs[1, group_coefs[[1]]])),
#                                       .x))]
# 
# na_plus <- function(x, y){nafill(x, "const", 0) + nafill(y, "const", 0)}
# 
# wl_coefs[, full_coefs := map2(pop_coefs,
#                               group_coefs,
#                               na_plus)]
# 
# wl_coefs[, full_coefs := map2(pop_coefs, full_coefs,
#                               ~set_names(.y, names(.x)))]
# wl_coefs[, f_predict := map(full_coefs,
#                             ~as.function(
#                               list(water.level = NULL,
#                                    pet = NULL,
#                                    precip = NULL,
#                                    bquote(
#                                      {precip1 <- 
#                                        shift(precip, 1, 0)
#                                      
#                                      precip2 <- 
#                                        shift(precip, 2, 0)
#                                      
#                                      precip3 <- 
#                                        shift(precip, 3, 0)
#                                      
#                                      pet1 <- 
#                                        shift(pet, 1, 0)
#                                      
#                                      pet2 <- 
#                                        shift(pet, 2, 0)
#                                      
#                                      pet3 <- 
#                                        shift(pet, 3, 0)
#                                      
#                                      .(wl_m) * water.level +
#                                        .(p_m) * precip +
#                                        .(p1_m) * precip1 +
#                                        .(p2_m) * precip2 +
#                                        .(p3_m) * precip3 +
#                                        .(pet_m) * pet +
#                                        .(pet1_m) * pet1 +
#                                        .(pet2_m) * pet2 +
#                                        .(pet3_m) * pet3
#                                      },
#                                      where = list(wl_m = .x[["wl_initial_cm"]],
#                                                   p_m = .x[["best_precip_cm"]],
#                                                   p1_m = .x[["lag1_precip_cm"]],
#                                                   p2_m = .x[["lag2_precip_cm"]],
#                                                   p3_m = .x[["lag3_precip_cm"]],
#                                                   pet_m = .x[["pet_cm"]],
#                                                   pet1_m = .x[["lag1_pet_cm"]],
#                                                   pet2_m = .x[["lag2_pet_cm"]],
#                                                   pet3_m = .x[["lag3_pet_cm"]])
#                                    )
#                               )))]


ds_mmods[, coefs := map(mod, coef)]

ds_mmods[, f_predict := map(coefs,
                            ~as.function(
                              list(pet = NULL,
                                   precip = NULL,
                                   bquote(
                                     {precip1 <- 
                                       shift(precip, 1, 0)
                                     
                                     precip2 <- 
                                       shift(precip, 2, 0)
                                     
                                     precip3 <- 
                                       shift(precip, 3, 0)
                                     
                                     pet1 <- 
                                       shift(pet, 1, 0)
                                     
                                     pet2 <- 
                                       shift(pet, 2, 0)
                                     
                                     pet3 <- 
                                       shift(pet, 3, 0)
                                     
                                     .(b) +
                                       .(p_m) * precip +
                                       .(p1_m) * precip1 +
                                       .(p2_m) * precip2 +
                                       .(p3_m) * precip3 +
                                       .(pet_m) * pet +
                                       .(pet1_m) * pet1 +
                                       .(pet2_m) * pet2 +
                                       .(pet3_m) * pet3
                                     },
                                     where = list(b = .x[["(Intercept)"]],
                                                  p_m = .x[["best_precip_cm"]],
                                                  p1_m = .x[["lag1_precip_cm"]],
                                                  p2_m = .x[["lag2_precip_cm"]],
                                                  p3_m = .x[["lag3_precip_cm"]],
                                                  pet_m = .x[["pet_cm"]],
                                                  pet1_m = .x[["lag1_pet_cm"]],
                                                  pet2_m = .x[["lag2_pet_cm"]],
                                                  pet3_m = .x[["lag3_pet_cm"]])
                                   )
                              )))]


tr_dat[, Ds_hat := ds_mmods[CJ(.BY[[1]], .BY[[2]]), f_predict[[1]]](pet_cm, best_precip_cm),
       by = .(site, site_status)]

ggplot(tr_dat,
       aes(x = Ds_cm,
           y = Ds_hat)) + 
  geom_point() + 
  geom_abline() +
  facet_wrap(~site,
             scales = "free")

ggplot(tr_dat,
       aes(x = Ds_cm,
           y = Ds_hat)) + 
  geom_point() + 
  geom_abline() +
  facet_wrap(~site) +
  coord_cartesian(xlim = c(-1.5, 1.5),
                  ylim = c(-1.5, 1.5))

tr_dat[, c("ytd_Ds_cm", "ytd_Ds_hat") := map(.SD, ~cumsum(nafill(.x, "const", 0))),
       by = .(site, sample_year),
       .SDcols = c("Ds_cm", "Ds_hat")]
tr_dat[, wl_hat := shift(wl_initial_cm + Ds_hat, 1),
       by = .(site, sample_year)]

ggplot(tr_dat[site == "077"],
       aes(x = doy)) +
  geom_line(aes(y = ytd_Ds_cm),
            color = "gray45") + 
  geom_line(aes(y = ytd_Ds_hat),
            color = "blue",
            linetype = "dashed") +
  facet_wrap(~sample_year,
             scales = "free_y") +
  theme_bw()

val_dat <- 
  tar_read(validation_data)

val_dat[, Ds_hat := wl_coefs[CJ(.BY[[1]], .BY[[2]]), 
                            f_predict[[1]]](wl_initial_cm, pet_cm, best_precip_cm),
       by = .(site, site_status)]
val_dat[, Ds_con := wl_coefs[CJ(.BY[[1]], "Control"),
                             f_predict[[1]]](wl_initial_cm, pet_cm, best_precip_cm),
        by = .(site)]

val_dat[, wl_hat := shift(wl_initial_cm + Ds_hat, 1),
       by = .(site, sample_year)]

val_dat[, wl_hat_con := shift(wl_initial_cm + Ds_con, 1),
        by = .(site, sample_year)]

ggplot(val_dat[site == "077"],
       aes(x = doy)) +
  geom_line(aes(y = wl_initial_cm),
            color = "gray45") + 
  geom_line(aes(y = wl_hat),
            color = "blue",
            linetype = "dashed") +
  geom_line(aes(y = wl_hat_con),
            color = "red",
            linetype = "dotted") +
  facet_wrap(~sample_year,
             scales = "free_y") +
  theme_bw()

val_dat[ds_mmods, 
        f_predict := i.f_predict,
        on = c("site", "site_status")]


predict_wl <- 
  function(initial.wl,
           pet,
           precip,
           wl_function){
    
    if(inherits(wl_function, "list")){
      wl_function <- 
        wl_function[[1]]
    }
    
    if(length(initial.wl) > 1){
      first_wl <- 
        min(which(!is.na(initial.wl)))
      
      initial.wl <- 
        initial.wl[first_wl]
    }
    
    wl.out <- 
      rep(NA, length(pet))
    
    ds.out <- 
      rep(NA, length(pet))
    
    wl.out[first_wl] <- 
      initial.wl
    
    for(i in (first_wl + 1):length(wl.out)){
      ds.out[i] <- 
        wl_function(wl.out[i-1], pet[i-1], precip[i-1])
      
      wl.out[i] <- 
        wl.out[i-1] + ds.out[i]
    }
    
    data.table(ds.out, wl.out)
  }

val_dat2 <- 
  val_dat[doy >= 138]

val_dat2[, 
        c("ds_hat_funct", "wl_hat_func") := 
          predict_wl(initial.wl = wl_initial_cm,
                     pet = pet_cm,
                     precip = best_precip_cm,
                     wl_function = first(f_predict)),
        by = .(site, sample_year)]

ggplot(val_dat2[site == "135"],
       aes(x = doy)) +
  geom_line(aes(y = wl_initial_cm),
            color = "gray45") + 
  geom_line(aes(y = wl_hat_func),
            color = "red",
            linetype = "dotted") +
  facet_wrap(~sample_year,
             scales = "free_y") +
  theme_bw()

