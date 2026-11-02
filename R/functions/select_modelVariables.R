#' Function to clip model variable layers to accessible area
#' Funciton is nested within a function that chooses what variables to keep 
#' by selecting the variables that lead to a VIF below a user specified value.
#' When choosing between two correlated variables, the model selects variable
#' that contributed the most to the maxent model when all variable were used 
#' in the model

# First we have the function to be nested.
library(usdm, quietly = TRUE)
remove_singleVariables <- function(pred_vars, maxent_mod, maxVIF){

  #calculate vifs of predictor variables
  pVIF <- usdm::vif(pred_vars)
  #calculate variable importance of maxent models
  m_results <- as.data.frame(as.table(maxent_mod@results)) %>% 
    dplyr::rename(variables = 1, rem = 2, permutation.importance = 3) %>% 
    dplyr::select(variables, permutation.importance)
  vIMP <- m_results %>% 
    dplyr::filter(str_detect(variables, '.permutation.importance')) %>% 
    dplyr::mutate(Variables = stringr::word(variables,  sep = fixed(".")))
  #join var importance with vifs
  jdf <- left_join(pVIF, vIMP, by = "Variables")
  # select the variables with the highest VIFs and select least important variable
  # based on permutation importance
  lowVar <- jdf %>% 
    dplyr::filter(VIF > maxVIF) %>% 
    dplyr::filter(VIF == sort(VIF, decreasing = TRUE)[1] |
             VIF == sort(VIF, decreasing = TRUE)[2]) %>% 
    dplyr::filter(permutation.importance == min(permutation.importance))
  pred_vars <- raster::stack(pred_vars)
  # make raster stack of variable without biggest VIF and lowest permutation imp
  message(as.character(lowVar$Variables))
  predictors_notInfl <- raster::dropLayer(pred_vars, as.character(lowVar$Variables))
  predictors_notInfl <- rast(predictors_notInfl)
  #p2_vifs <- usdm::vif(predictors_notInfl)
  
  return(predictors_notInfl)
}

# Function to automatically select SDM variables

select_sdmVariables <- function(pred_vars, maxent_mod, maxVIF){
  vv <- rast(pred_vars)
  
  while(max(usdm::vif(vv)$VIF) >= maxVIF){
    vv <- remove_singleVariables(pred_vars = vv,
                                 maxent_mod = maxent_mod,
                                 maxVIF = maxVIF)
  }
  
  return(vv)
}
