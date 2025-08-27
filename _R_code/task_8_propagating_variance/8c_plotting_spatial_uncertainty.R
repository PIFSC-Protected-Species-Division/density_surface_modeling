library(here)
library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(picMaps)
library(mapview)
library(tidyterra)
library(ggspatial)
library(ggpubr)
library(nmfspalette)

local_wd <- here("_R_code","task_8_propagating_variance")

#' -----------------------------------------------------------------------------
#' Load varprop prediction output file names
#' -----------------------------------------------------------------------------

pred_files <- list.files(file.path(local_wd, "output","predictions"))
est_files <- pred_files[grep("est", pred_files)]
var_files <- pred_files[grep("var", pred_files)]

#' -----------------------------------------------------------------------------
#' Calculate yearly spatial CV rasters
#' -----------------------------------------------------------------------------

years <- NULL
cv_raster <- NULL
mean_raster <- NULL
var_raster <- NULL

for(i in 1:length(est_files)){
  est_r <- rast(file.path(local_wd, "output","predictions",est_files[i]))
  years[i] <- time(est_r)[1] %>% year()
  var_r <- rast(file.path(local_wd, "output","predictions",var_files[i]))
  var_raster[[i]] <- app(var_r, mean) + app(est_r,var)
  mean_raster[[i]] <- app(est_r,mean)
  cv_raster[[i]] <- app(var_raster[[i]], sqrt) / mean_raster[[i]]
}
cv_raster <- do.call(c, cv_raster)
names(cv_raster) <- years
mean_raster <- do.call(c, mean_raster)
names(mean_raster) <- years
var_raster <- do.call(c, var_raster)
names(var_raster) <- years

var_20_24 <- app(var_raster[[-1]], mean) + app(mean_raster[[-1]], var)
mean_20_24 <- app(mean_raster[[-1]], mean)
cv_20_24 <- app(var_20_24, sqrt)/mean_20_24
cv_raster <- c(cv_20_24, cv_raster)
names(cv_raster)[1] <- "Average 2020-2024"

writeRaster(cv_raster, file.path(local_wd, "output","yearly_cv.tif"), overwrite=T)

#' -----------------------------------------------------------------------------
#' Create yearly CV plots
#' -----------------------------------------------------------------------------
hi <- picMaps::hawaii_coast()
eez <- picMaps::hawaii_eez()
fkwm <- picMaps::pfkw_mgmt()

ppp <- vector("list",nlyr(cv_raster))

ppp[[1]] <- ggplot() + geom_spatraster(data=cv_raster[[1]]) +
  layer_spatial(hi, fill="white", color=NA) +
  layer_spatial(eez, color="white", fill=NA, linetype="dashed",lwd=0.5) +
  layer_spatial(fkwm, color="white", fill=NA,lwd=0.5) +
  scale_fill_viridis_b(option = "E", na.value = NA) + 
  # scale_fill_princess_c("maori",name="% Extrapolation") +
  # scale_fill_grass_c("reds",name="CV") +
  # scale_fill_nmfs(palette = "coral", discrete = FALSE, name="CV", na.value=NA) + 
  # scale_fill_gradientn(colors=c("white","bisque", "red","darkred"), name="CV", na.value = NA) +
  scale_x_continuous(breaks=c(175, 228)) + scale_y_continuous(breaks=c(0, 40)) + 
  theme_bw() + #theme(legend.position="none") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  coord_sf(expand = FALSE) + 
  ggtitle("Average 2020-2024")

lll <- get_legend(ppp[[1]])
ppp[[1]] <- ppp[[1]] + theme(legend.position="none")


for(i in 2:nlyr(cv_raster)){
  ppp[[i]] <- ggplot() + geom_spatraster(data=cv_raster[[i]]) +
    layer_spatial(hi, fill="white", color=NA) +
    layer_spatial(eez, color="white", fill=NA, linetype="dashed",lwd=0.5) +
    layer_spatial(fkwm, color="white", fill=NA,lwd=0.5) +
    scale_fill_viridis_b(option = "E", na.value = NA) + 
    # scale_fill_princess_c("maori",name="% Extrapolation") +
    # scale_fill_grass_c("reds",name="CV") +
    # scale_fill_nmfs(palette = "coral", discrete = FALSE, name="CV", na.value=NA) + 
    # scale_fill_gradientn(colors=c("white","bisque", "red","darkred"), name="CV", na.value = NA) +
    scale_x_continuous(breaks=c(175, 228)) + scale_y_continuous(breaks=c(0, 40)) + 
    theme_bw() + theme(legend.position="none") +
    theme(axis.ticks = element_blank()) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    coord_sf(expand = FALSE) + 
    ggtitle(names(cv_raster)[i])
}

pout <- ggarrange(ppp[[1]], ppp[[2]], ppp[[3]], ppp[[4]], ppp[[5]], ppp[[6]], ppp[[7]],lll)
ggsave(pout, file=file.path(local_wd, "output","yearly_cv_fig.png"), width=6.5, height=6.5)
