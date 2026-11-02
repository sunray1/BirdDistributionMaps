library(terra, quietly = TRUE)

### Function to save SDM model results of interest
save_SDM_results <- function(ENMeval_output, AUCmin, resultDir, spp, occ_df){
  
  spp <- stringr::str_replace(string = spp, pattern = " ", replacement = "_")
  
  ### Find the best model
  bestmod <- ENMeval_output@results %>% dplyr::filter(delta.AICc == 0) %>% slice(1)
  if(bestmod$auc.train >= AUCmin & bestmod$auc.val.avg >= AUCmin){
    bestmod <- bestmod[1,]
  } else{
    bestmod <- ENMeval_output@results %>% 
      dplyr::filter(auc.train == max(auc.train) | auc.val.avg  == max(auc.val.avg)) 
    bestmod <- bestmod[1,]
  }
  write.csv(bestmod, file = file.path(resultDir, paste0( spp, "_bestModel.csv")), row.names = F)
  
  ### Build the best model
  maxent_args <- as.character(bestmod$tune.args)
  r_best <- ENMeval_output@predictions[[maxent_args]]
  raster::writeRaster(x = r_best, filename = file.path(resultDir, paste0(spp, "_SDM.tif")), overwrite = TRUE)
  
  varimp <- ENMeval_output@variable.importance[[maxent_args]]
  write.csv(x = varimp,               
            file = file.path(resultDir, paste0(spp, "_variableImportance.csv")),
            row.names = F)
}
# example
#save_SDM_results(ENMeval_output = ENMeval_output, 
#                 AUCmin = 0.7,
#                 resultDir = "ModelResults/",spp = "Astragalus bourgovii",
#                 occ_df = ab)