library(sf)
library(terra)
library(ggplot2)
library(tidyterra)

source("R/00_setup.R")

# 1. Prepare the Oahu boundary
hi <- terra::vect(wd$basemap) |>
  sf::st_as_sf() |>
  sf::st_transform(4326)

bbox <- sf::st_bbox(
  c(xmin = -158.4, ymin = 21.2,
    xmax = -157.5, ymax = 21.8),
  crs = sf::st_crs(4326)
)

oahu <- sf::st_crop(hi, bbox) |>
  sf::st_union()

# 2. Select species
wanted_species <- c(
  "Aidemedia_chascax",
  "Aidemedia_zanclops",
  "Akialoa_ellisiana",
  "Akialoa_upupirostris",
  "Asio_flammeus_sandwichensis",
  "Buteo_solitarius",
  "Carpodacus_sp.",
  "Chaetoptila_angustipluma",
  "Chasiempis_ibidis",
  "Chloridops_regiskongi",
  "Chloridops_wahi",
  "Chlorodrepanis_flava",
  "Circus_dossenus",
  "Corvus_impluviatus",
  "Corvus_viriosus",
  "Falco_peregrinus",
  "Fulica_alai",
  "Gallinula_galeata",
  "Grallistrix_orion",
  "Haliaeetus_albicilla",
  "Hemignathus_lucidus",
  "Himatione_sanguinea",
  "Loxioides_bailleui",
  "Loxops_wolstenholmei",
  "Melamprosops_phaeosoma",
  "Moho_apicalis",
  "Myadestes_obscurus",
  "Myadestes_woahensis",
  "Paroreomyza_maculata",
  "Psittirostra_psittacea",
  "Rhodacanthis_flaviceps",
  "Rhodacanthis_litotes",
  "Telespiza_persecutrix",
  "Telespiza_ultima",
  "Vestiaria_coccinea",
  "Xestospiza_fastigialis",
  "Zapornia_ralphorum",
  "Zapornia_ziegleri"
)

range_files <- list.files(
  path = wd$out,
  pattern = "_range\\.shp$",
  recursive = TRUE,
  full.names = TRUE
)

model_species <- sub("_range\\.shp$", "", basename(range_files))
model_species <- gsub("[[:space:]]+", "_", model_species)

# Match the existing peregrine subspecies model to the species
model_species[
  model_species == "Falco_peregrinus_tundrius"
] <- "Falco_peregrinus"

keep <- model_species %in% wanted_species
range_files <- range_files[keep]
model_species <- model_species[keep]

if (length(range_files) == 0) {
  stop("No matching range files found.")
}

# 3. Define a common raster grid
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

template <- terra::extend(template, terra::ext(oahu_vect))
template <- terra::crop(template, oahu_vect)
terra::values(template) <- 0

# 4. Stack binary ranges to calculate richness
richness <- template

for (spp in unique(model_species)) {
  
  message("Adding: ", spp)
  
  # Combine any matching files so each species counts only once
  species_presence <- template
  
  for (f in range_files[model_species == spp]) {
    
    species_range <- terra::vect(f)
    species_range <- terra::project(
      species_range,
      terra::crs(template)
    )
    
    presence <- terra::rasterize(
      species_range,
      template,
      field = 1,
      background = 0
    )
    
    species_presence <- (species_presence + presence) > 0
  }
  
  richness <- richness + species_presence
}

richness <- terra::mask(richness, oahu_vect)
names(richness) <- "richness"

# 5. Save the richness raster
result_dir <- file.path(wd$figures, "richness_native")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

terra::writeRaster(
  richness,
  filename = file.path(result_dir, "oahu_native_richness.tif"),
  overwrite = TRUE
)

# 6. Plot the richness map
p <- ggplot() +
  tidyterra::geom_spatraster(data = richness) +
  geom_sf(data = oahu, fill = NA, linewidth = 0.3) +
  scale_fill_viridis_c(
    name = "Species richness",
    na.value = "transparent",
    breaks = function(x) {
      b <- pretty(x)
      b[b >= 0 & b == floor(b)]
    }
  ) +
  coord_sf() +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Predicted native species richness on Oahu")

print(p)

ggsave(
  filename = file.path(result_dir, "oahu_native_richness.jpg"),
  plot = p,
  dpi = 600,
  width = 7,
  height = 6
)