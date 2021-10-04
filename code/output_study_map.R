source("code/load_project.R")
library(rnaturalearth)
library(sf)
library(ggplot2)
library(emojifont)
library(patchwork)

plot_bbox <- st_bbox(c(xmin = -93.5, ymin = 37.5, xmax = -65, ymax = 50),
                     crs = 4326)
lakes_file <- "data/spatial/natural_earth_great_lakes.gpkg"
hypso_file <- "data/spatial/natural_earth_hypsotint.tif"

if(fs::file_exists(lakes_file)){
  great_lakes <- 
    read_sf(lakes_file)
} else {
  great_lakes <- 
    ne_download(scale = 10, 
                type = 'lakes', 
                category = 'physical',
                load = TRUE,
                returnclass = 'sf') %>% 
    subset(name_alt == "Great Lakes") %>% 
    subset(grepl("Lake", name)) %>% 
    write_sf(lakes_file)
}


boundaries <- 
  ne_states(country = c("United States of America", "Canada"),
            returnclass = 'sf')

if(fs::file_exists(hypso_file)){
  hypso <- 
    stars::read_stars(hypso_file)
} else {
  tempfile <- 
    tempfile(fileext = ".zip")
  
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/HYP_50M_SR.zip",
                tempfile)
  
  hypso_prox <- 
    stars::read_stars(paste0("/vsizip/", tempfile, "/HYP_50M_SR/HYP_50M_SR.tif"))
  
  hypso <- 
    hypso_prox[plot_bbox]
  
  stars::write_stars(hypso, hypso_file)
}

study_coords <- 
  rbind(unique(tar_read(water_budget)[, .(type = "Wetland", station_name = site, lat, lon)]),
        unique(tar_read(external_met)[, .(type = "GHCND", station_name, lat, lon)]),
        unique(tar_read(mesowest_met)[, .(type = "Mesowest", station_name, lat, lon)]),
        unique(tar_read(swg_data)[station_name == "bergland_dam", .(type = "Norms_Station", station_name, lat, lon)])) %>% 
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326)

study_bounds <- 
  st_as_sfc(st_bbox(study_coords))


small_hypso <- hypso[study_bounds]

# Patchwork doesn't like inseting a plot with geom_sf() so converting to UTM and
# plotting as points and lines

utm_proj4 <- 
  "+proj=utm +zone=16 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

study_utm <- 
  study_coords %>% 
  st_transform(utm_proj4) %>% 
  cbind(., st_coordinates(.))

study_utm$symbol <- 
  fcase(study_utm$type == "Wetland", emoji('seedling'),
            study_utm$type == "GHCND", emoji('droplet'),
            study_utm$type == "Mesowest", emoji('sunny'),
            study_utm$type == "Norms_Station", emoji('sun_behind_rain_cloud'))

boundaries_utm <- 
  boundaries %>% 
  st_transform(utm_proj4) %>% 
  st_intersection(st_buffer(st_as_sfc(st_bbox(study_utm)), 5000))

boundaries_utm_coords <- 
  boundaries_utm %>% 
  as.data.table() %>% 
  .[, as.data.table(st_coordinates(geometry)), by = .(name)]
  
blue <- 
  palette.colors(palette = "Okabe-Ito")[6]


locus_map <- 
  ggplot() +
  # geom_sf(data = great_lakes) +
  geom_sf(data = st_simplify(boundaries, dTolerance = 0.05),
          fill = 'tan',
          color = 'white') +
  geom_sf(data = study_bounds) +
  coord_sf(xlim = c(-93.5, -65),
           ylim = c(37.5, 50)) +
  theme(panel.background = element_rect(fill = blue),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "#FFFFFFE6",
                                       color = NA))

detail_map <- 
  ggplot(study_utm) +
  aes(x = X, 
      y = Y) + 
  geom_polygon(data = boundaries_utm_coords,
               aes(group = name),
               fill = 'tan', 
               color = 'white') +
  geom_point(aes(shape = type),
             size = rel(5),
             fill = 'black') +
  scale_shape_manual(name = NULL,
                     values = c(Wetland = 19,
                                GHCND = 6,
                                Mesowest = 25,
                                Norms_Station = 8)) +
  labs(x = NULL,
       y = NULL) +
  scale_y_continuous(expand = expansion(-0.01, 0),
                     position = "right") +
  scale_x_continuous(expand = expansion(-0.01, 0),
                     position = "bottom") +
  theme(legend.position = "bottom", 
        panel.background = element_rect(fill = blue),
        panel.grid = element_blank(),
        legend.key = element_rect(fill = NA))

pub_map <- 
  detail_map + 
  inset_element(wrap_elements(plot = locus_map), 
                left = 0,
                bottom = 0.66, 
                right = 0.34, 
                top = 1)

ggsave(plot = pub_map,
       filename = "output/figures/climate_study_map.png",
       width = 12,
       height = 8,
       units = "in",
       dpi = 72,
       device = 'png')
