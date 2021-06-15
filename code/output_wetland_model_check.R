source("code/load_project.R")


dat <- 
  tar_read(test_data_fits)

ggplot(dat[!is.na(wl_initial_cm + wl_hat) & !(site %in% c("113", "119", "156", "157"))]) + 
  aes(x = wl_hat,
      y = wl_initial_cm) + 
  geom_point() +
  geom_abline(color = "red") + 
  facet_wrap(~site, 
             scales = "free")

ggplot(dat) + 
  aes(x = yday(sample_date)) + 
  geom_line(aes(y = wl_initial_cm),
            color = "blue") +
  geom_line(aes(y = wl_hat),
            color = "red",
            linetype = "dashed") +
  facet_wrap(~site + year(sample_date), 
             scales = "free")

res <- 
  dat[!is.na(wl_initial_cm + wl_hat)
      # & month(sample_date) %in% 5:9, 
      # lapply(.SD, mean, na.rm = TRUE),
      # .SDcols = c("wl_hat", "wl_initial_cm"),
      # by = .(site, water_year, sample_week = round_date(sample_date, "week"))][
        ,.(site_status = unique(site_status),
           rsr = hydroGOF::rsr(wl_hat, wl_initial_cm),
        nse = hydroGOF::NSE(wl_hat, wl_initial_cm),
        pbias = abs(hydroGOF::pbias(wl_hat, wl_initial_cm)),
        d = hydroGOF::d(wl_hat, wl_initial_cm),
        r2 = hydroGOF::rPearson(wl_hat, wl_initial_cm),
        br2 = hydroGOF::br2(wl_hat, wl_initial_cm),
        kge = hydroGOF::KGE(wl_hat, wl_initial_cm)),
      by = .(site, water_year)] %>% 
  melt(id.vars = c("site", "water_year", "site_status"),
       variable.name = "metric")

res_summary <- 
  res[, 
      .(metric, 
        site_status,
        value = fcase(metric == "rsr" & between(value, 0, 0.5), "Very Good", 
                      metric == "rsr" & between(value, 0.5, 0.6), "Good", 
                      metric == "rsr" & between(value, 0.6, 0.7), "Satisfactory", 
                      metric == "rsr" & value > 0.7, "Unsatisfactory",
                      metric == "nse" & between(value, 0.75, 1), "Very Good", 
                      metric == "nse" & between(value, 0.65, 0.75), "Good", 
                      metric == "nse" & between(value, 0.5, 0.65), "Satisfactory", 
                      metric == "nse" & value <= 0.5, "Unsatisfactory",
                      metric == "pbias" & between(value, 0, 10), "Very Good", 
                      metric == "pbias" & between(value, 10, 15), "Good", 
                      metric == "pbias" & between(value, 15, 25), "Satisfactory", 
                      metric == "pbias" & value > 25, "Unsatisfactory",
                      metric == "r2" & between(value, 0.85, 1), "Very Good", 
                      metric == "r2" & between(value, 0.75, 0.85), "Good", 
                      metric == "r2" & between(value, 0.6, 0.75), "Satisfactory", 
                      metric == "r2" & value <= 0.6, "Unsatisfactory",
                      metric == "d" & value > 0.9, "Very Good", 
                      metric == "d" & between(value, 0.85, 0.9), "Good", 
                      metric == "d" & between(value, 0.75, 0.85), "Satisfactory", 
                      metric == "d" & value <= 0.75, "Unsatisfactory")
        )
      ][, .N, by = .(metric, site_status, value)] %>% 
    dcast(value ~ metric + site_status, 
          value.var = "N") %>% 
  transform(value = factor(value,
                           levels = c("Very Good", "Good", "Satisfactory", "Unsatisfactory"),
                           ordered = TRUE)) %>% 
  .[order(value)]

res_summary
res[metric == "r2" & value < 0.6]

ggplot(res[metric == "r2"]) +
  aes(x = site,
      y = value) +
  geom_boxplot() +
  facet_wrap(~metric,
             ncol = 1,
             scales = "free")

ggplot(res) +
  aes(x = site,
      y = value) +
  geom_jitter() +
  geom_hline(data = data.frame(metric = c("rsr", "nse"), threshold = c()))
  facet_wrap(~metric,
             ncol = 1,
             scales = "free")
