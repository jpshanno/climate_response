#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param data
#' @param esy_data
#' @param output.file
#' @param ...
#' @return
#' @author Joe Shannon
#' @export
create_esy_impact_plot <- function(data, esy_data, parameters, output.file, ...) {

  dat <- copy(data)
  dat[parameters, on = "site", min_esy := params[[1]][["min_esy"]]]
  dat[esy_data, on = "site", esy_hat := esy_func[[1]](wl_initial_cm, min_esy)]
  dat[, `:=`(
    adj_pet = pet_cm * esy_hat,
    adj_rain = rain_cm * esy_hat,
    adj_melt = melt_cm * esy_hat)
    ]
  
  mod <- quantreg::rq(
    Ds_cm ~ 0 + pet_cm + rain_cm + melt_cm,
    data = dat)
  
  adj_mod <- quantreg::rq(
    Ds_cm ~ 0 + adj_pet + adj_rain + adj_melt,
    data = dat)
  
  aug_raw <- setDT(broom::augment(mod))
  r2_raw <- pretty_round(cor(x = aug_raw$Ds_cm, y = aug_raw$.fitted), 2)
  setnames(aug_raw, c("Ds_cm", ".fitted"), c("Observed Daily Change in Water Level (cm)", "Predicted Daily Change in Water Level (cm)"))
  
  aug_adj <- setDT(broom::augment(adj_mod))
  r2_adj <- pretty_round(cor(x = aug_adj$Ds_cm, y = aug_adj$.fitted), 2)
  setnames(aug_adj, c("Ds_cm", ".fitted"), c("Observed Daily Change in Water Level (cm)", "Predicted Daily Change in Water Level (cm)"))
  
  adj_regression <- quantreg::rq(`Predicted Daily Change in Water Level (cm)` ~ `Observed Daily Change in Water Level (cm)`, data = aug_adj)
  adj_coefs <- pretty_round(coef(adj_regression), 2)

  raw_regression <- quantreg::rq(`Predicted Daily Change in Water Level (cm)` ~ `Observed Daily Change in Water Level (cm)`, data = aug_raw)
  raw_coefs <- pretty_round(coef(raw_regression), 2)
  
  raw_plot <- 
    ggplot(data = aug_raw) +
    ggtitle("A. Raw Inputs") +
    aes(x = `Observed Daily Change in Water Level (cm)`, y = `Predicted Daily Change in Water Level (cm)`) +
    geom_point(color = "gray30", alpha = 0.5, stroke = 0) +
    # Add invisible points to force equal ylims
    geom_point(data = aug_adj, color = NA) +
    geom_abline(aes(slope = 1, intercept = 0), size = 10/14, color = "black", linetype = "dashed") +
    geom_quantile(quantiles = 0.5, formula = y ~ x, color = 'black') +
    annotate(x = min(c(aug_adj[["Observed Daily Change in Water Level (cm)"]], aug_raw[["Observed Daily Change in Water Level (cm)"]])),
             y = max(c(aug_adj[["Predicted Daily Change in Water Level (cm)"]], aug_raw[["Predicted Daily Change in Water Level (cm)"]])),
             label = paste0("pseudo-R<sup>2</sup> = ", r2_raw, "<br> Regression: y = ", raw_coefs[1], " +", raw_coefs[2], " * x"),
             fill = "#ffffffe6",
             geom = ggtext::GeomRichText,
             size = 3,
             hjust = 0,
             vjust = 1) + 
    ggforce::facet_zoom(xlim = c(-5, 5), ylim = c(-5, 5), zoom.size = 1) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(fill = NA, color = "black"),
      panel.ontop = FALSE,
      strip.background = element_rect(fill = "gray85", color = NA),
      plot.title = ggtext::element_markdown())
  
  esy_plot <- 
    ggplot(data = aug_adj) +
    ggtitle("B. E<sub>Sy</sub>-Adjusted Inputs") +
    aes(x = `Observed Daily Change in Water Level (cm)`, y = `Predicted Daily Change in Water Level (cm)`) +
    geom_point(color = "gray30", alpha = 0.5, stroke = 0) +
    # Add invisible points to force equal ylims
    geom_point(data = aug_raw, color = NA) +
    geom_abline(aes(slope = 1, intercept = 0), size = 10/14, color = "black", linetype = "dashed") +
    geom_quantile(quantiles = 0.5, formula = y ~ x, color = 'black') +
    annotate(x = min(c(aug_adj[["Observed Daily Change in Water Level (cm)"]], aug_raw[["Observed Daily Change in Water Level (cm)"]])),
             y = max(c(aug_adj[["Predicted Daily Change in Water Level (cm)"]], aug_raw[["Predicted Daily Change in Water Level (cm)"]])),
             label = paste0("pseudo-R<sup>2</sup> = ", r2_adj, "<br> Regression: y = ", adj_coefs[1], " +", adj_coefs[2], " * x"),
             fill = "#ffffffe6",
             geom = ggtext::GeomRichText,
             size = 3,
             hjust = 0,
             vjust = 1) + 
    ggforce::facet_zoom(xlim = c(-5, 5), ylim = c(-5, 5), zoom.size = 1) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(fill = NA, color = "black"),
      panel.ontop = FALSE,
      strip.background = element_rect(fill = "gray85", color = NA),
      plot.title = ggtext::element_markdown())
  
  fig <- {wrap_elements(raw_plot)} / {wrap_elements(esy_plot)}
  
  ggsave(plot = fig,
         filename = output.file,
         ...)
  
  output.file

}
