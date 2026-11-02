library(sf)
library(terra)
library(ggplot2)
library(tidyterra)

source("R/00_setup.R")

# 1. Prepare the Oahu boundary
hi <- terra::vect(wd$basemap) %>%
  sf::st_as_sf() %>%
  sf::st_transform(4326)

bbox <- sf::st_bbox(
  c(xmin = -158.4, ymin = 21.2,
    xmax = -157.5, ymax = 21.8),
  crs = sf::st_crs(4326)
)

oahu <- sf::st_crop(hi, bbox) %>%
  sf::st_union()

# 2. Find all completed range polygons
range_files <- list.files(
  path = wd$out,
  pattern = "_range\\.shp$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(range_files) == 0) {
  stop("No _range.shp files found in wd$out.")
}

# Use an existing SDM raster to set the grid resolution and CRS
model_files <- list.files(
  path = wd$out,
  pattern = "_SDM\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(model_files) == 0) {
  stop("No _SDM.tif files found to define the output grid.")
}

template <- terra::rast(model_files[1])[[1]]

oahu <- sf::st_transform(oahu, terra::crs(template))
oahu_vect <- terra::vect(oahu)

# Extend the template grid across all of Oahu, retaining its alignment
template <- terra::extend(template, terra::ext(oahu_vect))
template <- terra::crop(template, oahu_vect)
terra::values(template) <- 0

# 3. Add each modeled range to the richness map
richness <- template

for (f in range_files) {
  
  message("Adding: ", basename(f))
  
  species_range <- terra::vect(f)
  species_range <- terra::project(
    species_range,
    terra::crs(template)
  )
  
  # Cells whose centers fall inside the range receive 1
  presence <- terra::rasterize(
    species_range,
    template,
    field = 1,
    background = 0
  )
  
  richness <- richness + presence
}

# Restrict the result to Oahu's land area
richness <- terra::mask(richness, oahu_vect)
names(richness) <- "richness"

# 4. Save the richness raster
result_dir <- file.path(wd$figures, "richness")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

terra::writeRaster(
  richness,
  filename = file.path(result_dir, "oahu_richness.tif"),
  overwrite = TRUE
)

# 5. Plot and save
p <- ggplot() +
  tidyterra::geom_spatraster(data = richness) +
  geom_sf(data = oahu, fill = NA, linewidth = 0.3) +
  scale_fill_viridis_c(
    name = "Species richness",
    na.value = "transparent",
    breaks = scales::breaks_pretty()
  ) +
  coord_sf() +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Predicted species richness on Oahu")

print(p)

ggsave(
  filename = file.path(result_dir, "oahu_richness.jpg"),
  plot = p,
  dpi = 600,
  width = 7,
  height = 6
)
