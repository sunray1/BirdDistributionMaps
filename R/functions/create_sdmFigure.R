################ This script will create a figure to visualize the models

### Load libraries
library(ggplot2, quietly = TRUE)
library(egg , quietly = TRUE)

create_sdmFigure <- function(world = world, spp, r, r_pa, occ_df, resultDir = resultDir, bw = bw){
  spp = gsub("_", " ", spp)
  tiledf <- as.data.frame(r, xy = T) %>% na.omit() %>% rename(ClogLog = 3)
  tiledf_pa <- as.data.frame(r_pa, xy = T) %>% na.omit() %>% rename(ClogLog = 3)
  
  if(nrow(occ_df > 5000)){
    
    ## Create a plot for the model
    p1 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(tiledf,
                mapping = aes(x = x, y = y, fill = ClogLog), alpha = 0.9) + 
      geom_point(occ_df, mapping = aes(x = x, y = y), color = 'black', shape = 21, alpha = 0.15) +
      coord_sf(xlim = c(min(tiledf$x) - 500, max(tiledf$x) + 500),
               ylim = c(min(tiledf$y) - 500, max(tiledf$y) + 500)) +
      scale_fill_viridis_c() +
      ggtitle(paste0(spp, " SDM", " (# coords = ", nrow(occ_df), ")")) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title = element_blank()
      )
    
    ## Create a plot for the presence-absence model
    p2 <- ggplot() +
      geom_sf(world, mapping = aes()) +
      geom_tile(tiledf_pa,
                mapping = aes(x = x, y = y, fill = as.character(ClogLog)), alpha = 0.9) + 
      geom_point(occ_df, mapping = aes(x = x, y = y), color = 'black', shape = 21, alpha = 0.15) +
      coord_sf(xlim = c(min(tiledf$x) - 500, max(tiledf$x) + 500),
               ylim = c(min(tiledf$y) - 500, max(tiledf$y) + 500)) +
      labs(fill = "Presence") +
      scale_fill_viridis_d() +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title = element_blank()
      ) 
    
    ## Combine these two plots into one figure
    e <- gridExtra::arrangeGrob(grobs = list(p1, p2), nrow = 2)
    
    ggsave(
      plot = e,
      filename = file.path(resultDir, paste0(bw,"_map.jpg")), 
      width = 8, height = 5
    )
    
  } else{
    
    ## Create a plot for the model
    p1 <- ggplot() +
      geom_sf(world, mapping = aes(), fill = NA) +
      geom_tile(tiledf,
                mapping = aes(x = x, y = y, fill = ClogLog)) + 
      geom_point(occ_df, mapping = aes(x = x, y = y), color = 'black', shape = 21, alpha = 0.5) +
      coord_sf(xlim = c(min(tiledf$x) - 500, max(tiledf$x) + 500),
               ylim = c(min(tiledf$y) - 500, max(tiledf$y) + 500)) +
      scale_fill_viridis_c() +
      ggtitle(paste0(spp, " SDM", " (# coords = ", nrow(occ_df), ")")) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title = element_blank()
      )
    
    ## Create a plot for the presence-absence model
    p2 <- ggplot() +
      geom_sf(world, mapping = aes()) +
      geom_tile(tiledf_pa,
                mapping = aes(x = x, y = y, fill = as.character(ClogLog))) + 
      geom_point(occ_df, mapping = aes(x = x, y = y), color = 'black', shape = 21, alpha = 0.5) +
      coord_sf(xlim = c(min(tiledf$x) - 500, max(tiledf$x) + 500),
               ylim = c(min(tiledf$y) - 500, max(tiledf$y) + 500)) +
      labs(fill = "Presence") +
      scale_fill_viridis_d() +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title = element_blank()
      ) 
    
    ## Combine these two plots into one figure
    e <- gridExtra::arrangeGrob(grobs = list(p1, p2), nrow = 2)
    
    ggsave(
      plot = e,
      filename = file.path(resultDir, paste0(bw,"_map.jpg")), 
      width = 8, height = 5
    )
    
  }
  
}
