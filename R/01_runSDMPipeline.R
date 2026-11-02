############################ Set up the pipeline #############################

## Step 1: Load in the necessary libraries ----
source('R/00_setup.R')
library(rnaturalearth)
#library(rgeos)
library(dismo)
library(ENMeval)
library(stringr)
library(dplyr)
library(terra)
library(sf)
library(rJava)
library(smoothr)
library(units)
library(ggplot2)

## Step 2: Load in helper functions ----
list.files(wd$fun, full.names = T) %>% 
  lapply(source)

## Step 3: Prepare the spatial data ----
### 3A. load occurrence records
occs <- data.table::fread(wd$occs) %>% 
  dplyr::mutate(decimalLongitude = Corrected_longitude,
                decimalLatitude = Corrected_latitude) %>% 
  dplyr::filter(!is.na(decimalLongitude),
                !is.na(decimalLatitude)) %>% 
  dplyr::mutate(bw = str_replace(Search_Taxon_Name, " ", "_")) #bw stands for "Binomial With _ between genus and specific epithet"

### 3B. Read in basemap for visualizing
species_params <- read.csv(wd$species_params)

### 3C. Read in basemap for visualizing
world <- vect(x = wd$basemap) %>% st_as_sf()

### 3D. Choose projection and project data if necessary
study_proj <- proj4_aea #aea projection
world <- st_transform(world, crs = study_proj) 

occs_sf <- st_as_sf(occs, 
                    coords = c("decimalLongitude", "decimalLatitude"),
                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 
occs_sf <- st_transform(occs_sf, crs = study_proj)
occs <- occs %>% 
  mutate(x = st_coordinates(occs_sf)[,1], y = st_coordinates(occs_sf)[,2])

oahu <- sf::st_read("data/shapefiles/Oahu.shp", quiet = TRUE) %>%
  sf::st_transform(crs = study_proj) %>%
  sf::st_union()

# Match the Spatial format returned by your existing function
oahu_area <- sf::as_Spatial(sf::st_as_sf(oahu))

############################ Execute the pipeline #############################

rerunPipeline <- FALSE

#### Step 1: Set up species list
spp_list <- unique(occs$bw)

#### Step 2: Create function for executing the pipeline

sdm_pipeline <- function(x){
  
  print("=========================================================")
  print(paste("============== Working on", spp_list[x], "=============="))  
  print("=========================================================")
  binomial <- spp_list[x]
  species_df_raw <-   dplyr::filter(as.data.frame(occs), bw == binomial) %>% 
    filter(!is.na(x),
           !is.na(y))
  stopifnot(nrow(species_df_raw) >= 1)

  ### 2A. Create accessible area - Clip all potential predictor variables based on all cleaned occurrence records from this species.
  # For this project set aa to all of Oahu since birds are mobile
  aa_shp <- oahu_area
  mod_vars <- rast(wd$vars)
  
  ## Clip environmental variable layers to the defined accessible area
  mymod_vars <- clip_variableLayers(rstack = mod_vars, accessibleArea = aa_shp)

  ### 2B. Thin points based on the accessible area.
  ## Prepare the coordinates for rarefaction
  coordinates(species_df_raw) <- ~ x + y
  
  ## Perform the rarefaction
  species_df <- thinPoints(
    spp_df = species_df_raw, 
    area_sqkm = area_sqkm, 
    bio = mymod_vars[[2]], 
    method = "simple",
    simpleMult = 4
  )
  spp_df <- species_df_raw
  
  #### Step 3: Test and fine-tune the model
  
  ### 3A. Select top performing variables to reduce colinearity
  ## First run a test model
  print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
  max_model <- maxent(x = mymod_vars, p = coordinates(spp_df), progress = "text", silent = TRUE) 
  
  ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
  print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
  predictors <- select_sdmVariables(pred_vars = mymod_vars, maxent_mod = max_model, maxVIF = 5)
  
  ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
  # It is possible to change the number of cores
  print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
  
  eval1 <- ENMeval::ENMevaluate(
    occ = coordinates(spp_df), env = predictors,
    method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
    fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
    parallel = TRUE, numCores = 10, algorithm = 'maxent.jar',
    quiet = FALSE, updateProgress = TRUE
  )
  
  #### Step 4: Create the output for the best performing model
  
  ### 4A. Prepare the output path and coordinates
  bw <- spp_df$bw[1] # make a speciesName string without a space for better saving
  bw <- str_replace(bw, " ", "_")
  
  resultDir <- file.path(wd$out, bw)
  
  if (!dir.exists(resultDir)) {
      dir.create(resultDir, recursive = TRUE)
  }
  
  ### 4B. Output the model
  ## Output the model evaluation
  #  save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
  ## Output the best performing model, including both the actual model and the presence-absence model
  print(">>>>>>>>> Saving best model <<<<<<<<<")
  
  sp_df2 <- as.data.frame(coordinates(species_df_raw)) %>% 
    dplyr::rename(x = 1, y = 2)
  
  save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = resultDir,
                   spp = bw, occ_df = sp_df2)
  
  #### Step 5: Create LPT threshold cowplot
  ## Load in the rasters
  r    <- terra::rast(file.path(resultDir, paste0( bw,"_SDM.tif")))
  
  ## Classify rasters by a bunch of different thresholding values
  
  lpt_occs <- sp_df2 %>% 
      dplyr::select(x, y)
  lpt <- terra::extract(x = r, y = lpt_occs) 
  
  ## LPT 0  
  lptVal <- 0
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_000 <- terra::classify(r, pa_mat)
  
  ## LPT 0.01
  lptVal <- 0.01
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_001 <- terra::classify(r, pa_mat)
  
  ## LPT 0.025 
  lptVal <- 0.025
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_0025 <- terra::classify(r, pa_mat)
  
  ## LPT 0.05 
  lptVal <- 0.05
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_005 <- terra::classify(r, pa_mat)
  
  ## LPT 0.1  
  lptVal <- 0.1
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_01 <- terra::classify(r, pa_mat)
  
  ## LPT 0.15 
  lptVal <- 0.15
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_015 <- terra::classify(r, pa_mat)
  
  tiledf <- as.data.frame(r, xy = T) %>% na.omit() %>% rename(ClogLog = 3)
  
  p0 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(as.data.frame(r_best_pa_000,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
                mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
                alpha = 0.8) +
      geom_point(sp_df2, mapping = aes(x = x, y = y),
                 size = 1, alpha = 0.33, shape = 1) +
      coord_sf(xlim = c(min(sp_df2$x) - 5000, max(sp_df2$x) + 5000),
               ylim = c(min(sp_df2$y) - 5000, max(sp_df2$y) + 5000)) +
      labs(fill = "Presence") +
      scale_fill_manual(values = c("grey90", "gold")) +
      ggtitle("LPT 0%") +
      theme_classic()
  
  p001 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(as.data.frame(r_best_pa_001,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
                mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
                alpha = 0.8) +
      geom_point(sp_df2, mapping = aes(x = x, y = y),
                 size = 1, alpha = 0.33, shape = 1) +
      coord_sf(xlim = c(min(sp_df2$x) - 5000, max(sp_df2$x) + 5000),
               ylim = c(min(sp_df2$y) - 5000, max(sp_df2$y) + 5000)) +
      labs(fill = "Presence") +
      scale_fill_manual(values = c("grey90", "gold")) +
      ggtitle("LPT 1%") +
      theme_classic()
  
  p0025 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(as.data.frame(r_best_pa_0025,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
                mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
                alpha = 0.8) +
      geom_point(sp_df2, mapping = aes(x = x, y = y),
                 size = 1, alpha = 0.33, shape = 1) +
      coord_sf(xlim = c(min(sp_df2$x) - 5000, max(sp_df2$x) + 5000),
               ylim = c(min(sp_df2$y) - 5000, max(sp_df2$y) + 5000)) +
      labs(fill = "Presence") +
      scale_fill_manual(values = c("grey90", "gold")) +
      ggtitle("LPT 2.5%") +
      theme_classic()
  
  p005 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(as.data.frame(r_best_pa_005,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
                mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
                alpha = 0.8) +
      geom_point(sp_df2, mapping = aes(x = x, y = y),
                 size = 1, alpha = 0.33, shape = 1) +
      coord_sf(xlim = c(min(sp_df2$x) - 5000, max(sp_df2$x) + 5000),
               ylim = c(min(sp_df2$y) - 5000, max(sp_df2$y) + 5000)) +
      labs(fill = "Presence") +
      scale_fill_manual(values = c("grey90", "gold")) +
      ggtitle("LPT 5%") +
      theme_classic()
  
  p01 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(as.data.frame(r_best_pa_01,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
                mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
                alpha = 0.8) +
      geom_point(sp_df2, mapping = aes(x = x, y = y),
                 size = 1, alpha = 0.33, shape = 1) +
      coord_sf(xlim = c(min(sp_df2$x) - 5000, max(sp_df2$x) + 5000),
               ylim = c(min(sp_df2$y) - 5000, max(sp_df2$y) + 5000)) +
      labs(fill = "Presence") +
      scale_fill_manual(values = c("grey90", "gold")) +
      ggtitle("LPT 10%") +
      theme_classic()
  
  p015 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(as.data.frame(r_best_pa_015,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
                mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
                alpha = 0.8) +
      geom_point(sp_df2, mapping = aes(x = x, y = y),
                 size = 1, alpha = 0.33, shape = 1) +
      coord_sf(xlim = c(min(sp_df2$x) - 5000, max(sp_df2$x) + 5000),
               ylim = c(min(sp_df2$y) - 5000, max(sp_df2$y) + 5000)) +
      labs(fill = "Presence") +
      scale_fill_manual(values = c("grey90", "gold")) +
      ggtitle("LPT 15%") +
      theme_classic()
  
  cp <- cowplot::plot_grid(p0, p001, p0025, p005, p01, p015,
                           nrow = 3,
                           labels = c("A", "B", "C",
                                      "D", "E", "F"))
  
  ggsave(plot = cp, filename = file.path(resultDir, paste0(bw, "_lpt_thresholds.jpg")),
         width = 8, height = 6)
  
  #### Step 6: Create polygons from rasters
  print(">>>>>>>>> Creating range polygons <<<<<<<<<")
  
  lpt_val <- species_params %>%
      filter(species == binomial) %>%
      pull(lpt)
  
  rasterworking <- switch(
      as.character(lpt_val),
      "0"     = r_best_pa_000,
      "0.01"  = r_best_pa_001,
      "0.025" = r_best_pa_0025,
      "0.05"  = r_best_pa_005,
      "0.1"   = r_best_pa_01,
      "0.15"  = r_best_pa_015
  )
  
  raster_ones <- terra::clamp(rasterworking, lower=1, value=FALSE)
  
  dist_poly_raw <- as.polygons(raster_ones)
  dist_poly <- sf::st_as_sf(dist_poly_raw)
  
  #remove cell groups smaller than threshold
  dist_poly <- drop_crumbs(dist_poly, set_units(5, km^2))
  
  #fill holes
  dist_poly <- fill_holes(dist_poly, set_units(4, km^2))
  
  #smooth edges
  dist_poly <- smooth(dist_poly, method = "ksmooth", smoothness = 6)
  
  #plot for funsies
  #ggplot() +
  #    geom_sf(world, mapping = aes(), fill = "NA") +
  #    geom_sf(dist_poly, mapping = aes(), fill = "purple") +
  #    coord_sf(xlim = c(min(sp_df2$x) - 5000, max(sp_df2$x) + 5000),
  #             ylim = c(min(sp_df2$y) - 5000, max(sp_df2$y) + 5000)) +
  #    theme_classic()
  
  #write shapefile
  st_write(
      dist_poly,
      file.path(resultDir, paste0(binomial, "_range.shp")),
      driver = "ESRI Shapefile",
      delete_layer = TRUE
  )

  create_sdmFigure(
    spp = bw, r = r, r_pa = rasterworking, occ_df = sp_df2, world = world, resultDir = resultDir, bw = bw
  )
  
  print(paste(">>>>>>>>> Finished", binomial, "- outputs written to", resultDir, "<<<<<<<<<"))
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  gc()
  
}

# To run a specific species
#sdm_pipeline(1)

# Loop through pipeline for all species
for(i in seq_along(spp_list)){
  tryCatch(sdm_pipeline(x = i),
           error = function(e) print(paste(spp_list[i], "Error in Code Skipping for now!")))
}  

