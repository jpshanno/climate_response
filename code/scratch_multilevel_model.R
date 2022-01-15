source("code/load_project.R")
dat <- tar_read(training_data)[["control"]]
tar_load(esy_functions)
dat[esy_functions, max_wl := i.max_wl, on = "site"]

dat[
  j = l1_wl_initial_cm := shift(wl_initial_cm, 1),
  by = .(site, water_year)
]

dat[esy_functions, esy := i.pred_fun[[1]](l1_wl_initial_cm, min.esy = 1), on = "site"]

trimmed_dat <- 
  dat[!is.na(wl_initial_cm + l1_wl_initial_cm + pet_cm + rain_cm + melt_cm)]

int_form <- brmsformula(
  wl_initial_cm ~ int,
  int ~ 1 + 1 | site,
  nl = TRUE
)

l1_form <- brmsformula(
  wl_initial_cm ~ int + mwl * l1_wl_initial_cm,
  int + mwl ~ 1 + 1 | site,
  nl = TRUE
)

preds_form <- brmsformula(
  wl_initial_cm ~ int + mwl * l1_wl_initial_cm + mwb * (pet_cm + rain_cm) + mm * melt_cm,
  int + mwl + mwb + mm ~ 1 + 1 | site,
  nl = TRUE
)

maxwl_form <- brmsformula(
  wl_initial_cm ~ int + mwl * l1_wl_initial_cm +
    (mwb * (pet_cm + rain_cm) + mm * melt_cm) +
    (l1_wl_initial_cm > max_wl) * mq * (l1_wl_initial_cm - max_wl),
  int + mwl + mwb + mm + mq ~ 1 + 1 | site,
  nl = TRUE
)

ex_esy_form <- brmsformula(
  wl_initial_cm ~ int + mwl * l1_wl_initial_cm +
    (esy * mwb * (pet_cm + rain_cm) + mm * melt_cm) +
    (l1_wl_initial_cm > max_wl) * mq * (l1_wl_initial_cm - max_wl),
  int + mwl + mwb + mm + mq ~ 1 + 1 | site,
  nl = TRUE
)

# Internal ESy might be too complicated because esy_emp needs to be fit to a
# subset of the data
in_esy_form <- brmsformula(
  # FIX FORMULA TO MATCH maxwl_form
  wl_initial_cm ~ (wlT >= l1_wl_initial_cm) * (int + mwl * l1_wl_initial_cm + esy * mwb * (pet_cm + rain_cm) + mm * melt_cm) + 
    (wlT < l1_wl_initial_cm) * mq * (l1_wl_initial_cm - wlT),
  int + mwl + mwb + mm + wlT + mq ~ 1 + 1 | site,
  # This stuff is just copied over from wetland_model_functions.R
  wl_initial_cm ~ quad(ytd_water_balance, b0 = b0, b1 = b1, b2 = b2),
  esy_emp ~ quad_prime(mod = mod, wa = .SD[, ytd_water_balance]),
  esy_emp ~ a - (a - b) * exp (c * wl_initial_cm),
  nl = TRUE
)

int_mod <-
  brm(
    formula = int_form,
    data = trimmed_dat
  )

l1_mod <-
  brm(
    formula = l1_form,
    data = trimmed_dat
  )

preds_mod <-
  brm(
    formula = preds_form,
    data = trimmed_dat,
    cores = 4
  )

maxwl_mod <-
  brm(
    formula = maxwl_form,
    data = trimmed_dat,
    cores = 4
  )

external_esy_mod <-
  brm(
    formula = ex_esy_form,
    data = trimmed_dat,
    cores = 4
  )

internal_esy_mod <-
  brm(
    formula = in_esy_form,
    data = trimmed_dat,
    cores = 4
  )

loo_compare(loo(preds_mod), loo(maxwl_mod), loo(external_esy_mod))

pp_check(preds_mod)
pp_check(maxwl_mod)
pp_check(external_esy_mod)
plot(esy_mod)
dat[, sample_date2 := sample_date]
first_rows <- dat[!is.na(wl_initial_cm) & month(sample_date) %in% 3:5, .(first_date = min(sample_date)), by = .(site)]
test_dat <- copy(dat[first_rows, on = c("site" = "site", "sample_date2 >= first_date")])
test_dat[esy_functions, pred_fun := i.pred_fun, on = "site"]

create_function <- function(x) {
  as.function(
    list(
      l1_wl = NULL,
      esy = NULL,
      pet = NULL,
      rain = NULL,
      melt = NULL,
      bquote(
        .(int_Intercept) + .(mwl_Intercept) * l1_wl +
          esy * (.(mwb_Intercept) * (pet + rain) + .(mm_Intercept) * melt) + 
          ((l1_wl > .(max_wl))) * .(mq_Intercept)  * (l1_wl - .(max_wl)),
        where = x)
      )
    )
}

create_coef_df <- function(x) {
  coefs <- coef(x)$site

  sites <- rownames(coefs)

  purrr::map_dfr(
    sites,
    ~ as.data.table(c(site = .x, as.list(coefs[.x, , ]["Estimate", ])))
  )
}

coef_df <- create_coef_df(maxwl_mod)[esy_functions, on = "site"]
pred_functions <- coef_df[, .(fun = list(create_function(.SD[1,]))), keyby = .(site)]

predict_water_levels <- function(initial_wl, site, functions, esy_funs, covars) {

  pred_fun <- functions[site, fun][[1]]
  esy_fun <- esy_funs[site, pred_fun][[1]]

  .fitted <- numeric(nrow(covars))
  .fitted[1] <- initial_wl

  for(i in 2:length(.fitted)) {
    .fitted[i] <- pred_fun(
      l1_wl = .fitted[i-1],
      esy = esy_fun(.fitted[i-1], 1),
      pet = covars[i, pet_cm],
      rain = covars[i, rain_cm],
      melt = covars[i, melt_cm]
    )
  }
  .fitted
}

test_dat[, .fitted := predict_water_levels(.SD[1, wl_initial_cm], .BY[[1]], pred_functions, esy_functions, .SD), by = .(site, water_year)]

ggplot(test_dat) +
  aes(x = sample_date) +
  geom_line(aes(y = wl_initial_cm), linetype = "dashed") +
  geom_line(aes(y = .fitted), color = "blue") +
  facet_wrap(~site, scales = "free")
