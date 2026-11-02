library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggnewscale)
library(ggspatial)

source('R/00_setup.R')
study_proj <- proj4_aea
#read in HI shapefile
hi <- vect(wd$basemap) %>% 
    st_as_sf()

# let's restrict and shapefile to only oahu
bbox <- st_bbox(c(xmin = -158.4, xmax = -157.5, ymin = 21.2, ymax = 21.8))
oahu <- hi %>% 
    st_crop(bbox)
oahu <- st_transform(oahu, st_crs(study_proj))

#get species just on Oahu
occs <- data.table::fread(wd$occs) %>% 
    dplyr::mutate(decimalLongitude = Corrected_longitude,
                  decimalLatitude = Corrected_latitude) %>% 
    dplyr::filter(!is.na(decimalLongitude),
                  !is.na(decimalLatitude)) %>% 
    dplyr::mutate(bw = str_replace(Search_Taxon_Name, " ", "_"))
occs_sf <- st_as_sf(occs, 
                    coords = c("decimalLongitude", "decimalLatitude"),
                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 
occs_sf <- st_transform(occs_sf, crs = study_proj)

oahu_species <- occs_sf %>% st_intersection(oahu) %>% pull(bw) %>% unique()

# Load elevation data
elev2 <- geodata::elevation_3s( lon = -157.8, lat = 21.65, path = "data/elevation")
elev2 <- terra::project(elev2, study_proj)

elev_oahu <- crop(elev2, oahu)

#load rasters
SDM_rasters <- list.files(path = wd$out, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
SDM_rasters_oahu <- SDM_rasters[sapply(SDM_rasters, function(x) any(str_detect(x, oahu_species)))]

raster_df <- lapply(SDM_rasters_oahu, function(f) {
    r <- terra::rast(f)
    df <- terra::crop(r, vect(oahu)) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE)
    
    names(df)[3] <- "suitability"
    
    df %>%
        mutate(
            sppLab = basename(dirname(f)) %>%
                gsub("Auriculella_", "A. ", .)
        )
}) %>%
    bind_rows()

pp <- raster_df %>%
    mutate(sppLab = factor(sppLab, levels = sort(unique(sppLab)))) %>%
    ggplot() +
    geom_spatraster(data = elev_oahu, alpha = 0.5) +
    scale_fill_gradientn(
        colors = c("grey50", "grey90", "white"),
        na.value = NA,
        guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_sf(data = oahu, fill = NA) +
    geom_tile(aes(x = x, y = y, fill = suitability)) +
    facet_wrap(~sppLab, nrow = 6, ncol = 3) +
    scale_fill_viridis_c(
        name = "Habitat Suitability",
        na.value = NA
    ) +
    coord_sf() +
    theme_void() +
    theme(
        strip.text = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5)
    ) +
    guides(
        fill = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5
        )
    )

ggsave(pp, filename = file.path(wd$figures, "manuscript", "speciesDistributionMaps.jpg"), dpi = 600, width = 5, height = 7)
