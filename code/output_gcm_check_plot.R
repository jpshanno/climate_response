source('code/load_project.R')
library(ggdist)

green <- as.character(palette.colors()[4])
orange <- as.character(palette.colors()[7])

obs <- tar_read(swg_data)[station_name == "bergland_dam"]
# obs <- tar_read(swg_simulations_loca)[scenario == "historical"]
loca <- tar_read(loca_simulations)[station_name == "bergland_dam" & scenario == "historical"]


dat <- 
  rbind(obs[, .(type = 'Observed', sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)],
        loca[, .(type = 'LOCA', sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)])

obs[loca[gcm == "gfdl-cm3"], 
    `:=`(gfdl_tmin_c = i.tmin_c,
         gfdl_tmax_c = i.tmax_c,
         gfdl_precip_cm = i.precip_cm), 
    on = c("sample_date")]

obs[loca[gcm == "ccsm4"], 
    `:=`(ccsm_tmin_c = i.tmin_c,
         ccsm_tmax_c = i.tmax_c,
         ccsm_precip_cm = i.precip_cm), 
    on = c("sample_date")]

obs[, sample_season := as.climate_season(sample_month, TRUE)]

panel_out <- 
  {ggplot(obs[na.omit(obs)],
          aes(x = sample_season)) + 
      geom_jitter(aes(y = tmin_c,
                       fill = "Observed"),
                   width = 0.5) +
      geom_jitter(aes(y = gfdl_tmin_c,
                       fill = "LOCA"), 
                   width = 0.5, 
                   side = 'left')} +
  {ggplot(obs[na.omit(obs)],
          aes(x = sample_season)) + 
      stat_halfeye(aes(y = tmax_c,
                       fill = "Observed"),
                   width = 0.5) +
      stat_halfeye(aes(y = gfdl_tmax_c,
                       fill = "LOCA"), 
                   width = 0.5, 
                   side = 'left')} +
  {ggplot(obs[precip_cm > 0 & gfdl_precip_cm > 0],
          aes(x = sample_season)) + 
      stat_halfeye(aes(y = precip_cm,
                       fill = "Observed"),
                   width = 0.5) +
      stat_halfeye(aes(y = gfdl_precip_cm,
                       fill = "LOCA"), 
                   width = 0.5, 
                   side = 'left') +
      scale_y_log10()} +
  {ggplot(obs[na.omit(obs)],
          aes(x = sample_season)) + 
      stat_halfeye(aes(y = tmin_c,
                       fill = "Observed"),
                   width = 0.5) +
      stat_halfeye(aes(y = ccsm_tmin_c,
                       fill = "LOCA"), 
                   width = 0.5, 
                   side = 'left')} +
  {ggplot(obs[na.omit(obs)],
          aes(x = sample_season)) + 
      stat_halfeye(aes(y = tmax_c,
                       fill = "Observed"),
                   width = 0.5) +
      stat_halfeye(aes(y = ccsm_tmax_c,
                       fill = "LOCA"), 
                   width = 0.5, 
                   side = 'left')} +
  {ggplot(obs[precip_cm > 0 & ccsm_precip_cm > 0],
          aes(x = sample_season)) + 
      stat_halfeye(aes(y = precip_cm,
                       fill = "Observed"),
                   width = 0.5) +
      stat_halfeye(aes(y = ccsm_precip_cm,
                       fill = "LOCA"), 
                   width = 0.5, 
                   side = 'left') +
      scale_y_log10()} +
  plot_layout(ncol = 3,
              nrow = 2,
              guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  scale_fill_manual(values = c(Observed = green,
                               LOCA = orange)) &
  labs(x = NULL) &
  theme_bw(base_size = 16) &
  guides(fill = guide_legend(title = NULL)) &
  theme(legend.position = "bottom")

ggsave(plot = panel_out,
       filename = "output/figures/panel_slab_intervals_gcm_evaluations.png",
       width = 12,
       height = 8,
       dpi = 300)

dat[, sample_season := as.climate_season(sample_date, TRUE)]

ggplot(dat,
       aes(x = sample_season,
           y = tmin_c,
           fill = type)) + 
  geom_boxplot()
ggplot(dat,
       aes(x = sample_season,
           y = tmax_c,
           fill = type)) + 
  geom_boxplot()
ggplot(dat,
       aes(x = sample_season,
           y = precip_cm,
           fill = type)) + 
  geom_boxplot() +
  scale_y_log10()


ks.test(obs[as.climate_season(sample_month, FALSE) == "jja", tmin_c], loca[as.climate_season(sample_month, FALSE) == "jja", tmin_c])
plot(density(obs$tmin_c))
lines(density(loca$tmin_c), col = 'red')

curve(ecdf(obs[as.climate_season(sample_month, FALSE) == "jja", tmin_c])(x), form = -5, to = 25)
curve(ecdf(loca[as.climate_season(sample_month, FALSE) == "jja", tmin_c])(x), col = 'blue', add = TRUE)

goftest::ad.test(obs[as.climate_season(sample_month, FALSE) == "jja", tmin_c],
                 ecdf(loca[as.climate_season(sample_month, FALSE) == "jja", tmin_c]))

qqplot(obs[as.climate_season(sample_month, FALSE) == "jja", tmin_c],
       loca[as.climate_season(sample_month, FALSE) == "jja", tmin_c])
abline(0, 1)
