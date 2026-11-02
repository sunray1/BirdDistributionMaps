library(magrittr)
library(dplyr)
library(terra)
library(ggplot2)

my_dir_path <- getwd()
wd <- list()
wd$R         <- file.path( my_dir_path, "R" )
wd$fun       <- file.path( wd$R, "functions" )
wd$data      <- file.path( my_dir_path, "data" )
wd$figures   <- file.path( my_dir_path, "figures" )
wd$shapefiles<- file.path( wd$data, "shapefiles")
wd$basemap   <- file.path( wd$shapefiles, "HI.shp")
wd$occs      <- file.path( wd$data, "BPBM_Aves_Oahu.csv")
wd$vars      <- file.path( wd$data, "modelVariablesStack_filtered.tif")
wd$species_params <- file.path(wd$data, "species_parameters.csv")
wd$out       <- file.path( my_dir_path, "outputs" )

proj4_aea <- "+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-157 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
proj4_wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

