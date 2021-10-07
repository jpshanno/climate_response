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
  dat[parameters, on = "site", minESY := params[[1]][["minESY"]]]
  dat[esy_data, on = "site", esy_hat := pred_fun[[1]](wl_initial_cm, minESY)]
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
  
  
  
  raw_plot <- 
    ggplot(data = aug_raw) +
    ggtitle("A. Raw Inputs") +
    aes(x = `Observed Daily Change in Water Level (cm)`, y = `Predicted Daily Change in Water Level (cm)`) +
    annotate(x = min(aug_raw[["Observed Daily Change in Water Level (cm)"]]),
             y = max(aug_raw[["Predicted Daily Change in Water Level (cm)"]]),
             label = paste0("pseudo-R<sup>2</sup> = ", r2_raw),
             fill = "#ffffff",
             geom = ggtext::GeomRichText,
             hjust = 0,
             vjust = 1) + 
    geom_point(color = "gray30", alpha = 0.5, stroke = 0) +
    geom_abline(aes(slope = 1, intercept = 0), size = 10/14, color = "black", linetype = "dashed") +
    ggforce::facet_zoom(xlim = c(-5, 5), ylim = c(-5, 5), zoom.size = 1) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "gray85", color = NA),
      plot.title = ggtext::element_markdown())
  
  esy_plot <- 
    ggplot(data = aug_adj) +
    ggtitle("B. E<sub>Sy</sub> Inputs") +
    aes(x = `Observed Daily Change in Water Level (cm)`, y = `Predicted Daily Change in Water Level (cm)`) +
    annotate(x = min(aug_adj[["Observed Daily Change in Water Level (cm)"]]),
             y = max(aug_adj[["Predicted Daily Change in Water Level (cm)"]]),
             label = paste0("pseudo-R<sup>2</sup> = ", r2_adj),
             fill = "#ffffff",
             geom = ggtext::GeomRichText,
             hjust = 0,
             vjust = 1) + 
    geom_point(color = "gray30", alpha = 0.5, stroke = 0) +
    geom_abline(aes(slope = 1, intercept = 0), size = 10/14, color = "black", linetype = "dashed") +
    ggforce::facet_zoom(xlim = c(-5, 5), ylim = c(-5, 5), zoom.size = 1) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "gray85", color = NA),
      plot.title = ggtext::element_markdown())
  
  fig <- {wrap_elements(raw_plot)} / {wrap_elements(esy_plot)}
  
  ggsave(plot = fig,
         filename = output.file,
         ...)
  
  output.file

}
