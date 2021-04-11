##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param out.dir
##' @return
##' @author Joe Shannon
##' @export
create_model_evaulation_plots <- 
  function(data, out.dir) {
    
    fs::dir_create(out.dir)

    for(i in unique(data$site)){
      dat <- 
        data[site == i]
      fig <- 
        ggplot(dat) +
        aes(x = sample_date) +
        geom_line(aes(y = wl_initial_cm,
                      color = 'Observed',
                      linetype = 'Observed')) +
        geom_line(aes(y = wl_hat,
                      color = 'Modeled')) +
        facet_wrap(~water_year,
                   scales = "free_x") +
        scale_color_manual(values = c(Observed = 'gray40',
                                      Modeled = 'black')) +
        scale_linetype_manual(values = c(Observed = 'dashed',
                                         Modeled = 'solid')) +
        labs(y = "Water Level Relative to Ground (cm)") +
        theme_minimal(base_size = 14) +
        theme(legend.position = 'bottom',
              axis.title.x = element_blank())
      
      ggsave(filename = fs::path(out.dir, paste0('wetland_evaluation_linegraph_', i), ext = "pdf"),
             plot = fig,
             width = 7.5,
             height = 10,
             units = "in")
      
    }
    
    fit_fig <- 
      ggplot(data[!is.na(wl_hat + wl_initial_cm), 
                  .(site_status = site_status[1],
                    md = hydroGOF::md(wl_hat, wl_initial_cm)),
                  by = .(site, water_year)]) +
      aes(x = site,
          y = md) +
      geom_boxplot() +
      labs(y = "Modified d") +
      facet_wrap(~site_status) +
      theme_minimal(base_size = 14) +
      theme(axis.title.x = element_blank())
    
    ggsave(filename = fs::path(out.dir, "wetland_evaluation_summary_boxplot", ext = "pdf"),
           plot = fit_fig,
           width = 8,
           height = 4,
           units = "in")
    
    scatter_fig <- 
      ggplot(data) + 
      aes(x = wl_initial_cm, 
          y = wl_hat, 
          color = site_status) + 
      geom_point(alpha = 0.5,
                 shape = 20) +
      geom_abline(color = 'black',
                  linetype = 'dashed') + 
      geom_smooth(method = 'lm',
                  formula = y ~ x, se = FALSE) +
      labs(x = "Observed (cm)",
           y = "Modeled (cm)") +
      facet_wrap(~site, scales = 'free') +
      theme_minimal(base_size = 14) +
      theme(legend.position = 'bottom')
    
    ggsave(filename = fs::path(out.dir, "wetland_evaluation_scatterplot", ext = "pdf"),
           plot = scatter_fig,
           width = 12,
           height = 8,
           units = 'in')
    
    fs::dir_ls(out.dir)
    
  }
