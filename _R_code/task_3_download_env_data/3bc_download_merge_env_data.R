
#' -----------------------------------------------------------------------------
#' Load packages 
#' -----------------------------------------------------------------------------
library(here)
library(terra)
library(sf)
library(lubridate)
library(picMaps)
library(reticulate) 
library(tidyverse)
library(foreach)
library(foreach)
library(future)
library(doFuture)
library(progressr)

#' -----------------------------------------------------------------------------
#' Specify output directory for nc files 
#' -----------------------------------------------------------------------------

localwd <- here("_R_code","task_3_download_env_data")

copernicus_data <- file.path(localwd, "copernicus_data")
if(!dir.exists(copernicus_data)) dir.create(copernicus_data)

reticulate::use_virtualenv("cm_py", required = TRUE)
cm <- reticulate::import("copernicusmarine")


#' -----------------------------------------------------------------------------
#' Load processed survey segment data 
#' -----------------------------------------------------------------------------

segout <- readRDS(file.path(here(), "_R_code","task_2_est_g0_esw","output","seg_sight_out_g0_ESW.rds"))
segdata <- segout$segdata
segdata <- st_as_sf(segdata, coords=c("mlon","mlat"), crs=4326) %>% st_shift_longitude()
segdata <- segdata %>% mutate(
  UTC = as.POSIXct(paste(year, month, day, sep = "-"), format = "%Y-%m-%d", tz = "UTC") + hms(mtime)
) %>% arrange(UTC)

#' -----------------------------------------------------------------------------
#' Loop over dates to download 
#' -----------------------------------------------------------------------------

years <- unique(segdata$year)

plan("multisession", workers = 8) 

with_progress({
  pb <- progressor(along = years) 
  env_download <- foreach(i = seq_along(years), 
                          .errorhandling = "pass", 
                          .options.future = list(seed = TRUE)) %dofuture% 
    {
      
      ### Set up Python environment in each worker
      cm <- reticulate::import("copernicusmarine")
      lg <- reticulate::import("logging")
      lg$getLogger("copernicusmarine")$setLevel('WARN')
      
      tmp_seg <- dplyr::filter(segdata, year==years[i]) 
      bound <- st_bbox(tmp_seg) + c(-1/6, -1/6, 1/6, 1/6)
      time_span <- date(tmp_seg$UTC) %>% range() %>% as.character()
      
      dataset_id <- ifelse(tmp_seg$year[1]<=2021,
                           "cmems_mod_glo_phy_my_0.083deg_P1D-m",
                           "cmems_mod_glo_phy_myint_0.083deg_P1D-m")
      
      nc_dl <- cm$subset(
        dataset_id = dataset_id,
        variables = as.list(c("mlotst","so","thetao","zos")),
        minimum_longitude = bound$xmin,
        maximum_longitude = bound$xmax,
        minimum_latitude = bound$ymin,
        maximum_latitude = bound$ymax,
        start_datetime = time_span[1], 
        end_datetime = time_span[2], 
        minimum_depth = 0.5,
        maximum_depth = 0.5,
        coordinates_selection_method = "outside",
        dataset_version = "202311",
        disable_progress_bar = TRUE,
        skip_existing = TRUE,
        output_directory = copernicus_data
      )
      
      pb()
      
      data.frame(years=years[i], filename=nc_dl$filename)
    }
  
}) 
plan("sequential")
env_download <- do.call(rbind, env_download)
write.csv(env_download, file.path(copernicus_data, "env_download.csv"), row.names = FALSE)


#' -----------------------------------------------------------------------------
#' Loop over downloaded files to attach to segment data
#' -----------------------------------------------------------------------------
#' env_data <- read.csv(file.path(copernicus_data, "env_download.csv"))
years <- env_download$years

segdata_env <- foreach(i = seq_along(years), 
                       .errorhandling = "pass", 
                       .options.future = list(seed = TRUE)) %do% 
  {
    filenm <- env_download$filename[i]
    env <- terra::rast(file.path(copernicus_data, filenm))
    vnms <- terra::varnames(env)
    
    tmp_seg <- dplyr::filter(segdata, year==years[i]) 
    
    for(j in seq_along(vnms)){
      tmp_env <- 1.0*env[vnms[j]]
      tmp_env <- terra::project(tmp_env, "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs")
      tmp_env <- focal(tmp_env, 3, mean, na.policy="only", na.rm=TRUE) # Make sure there are no cells with missing data
      
      ### Diagnostic plot if needed
      # library(ggplot2); library(ggspatial)
      # ggplot() +
      #   layer_spatial(tmp_env[[1]]) +
      #   layer_spatial(tmp_seg)
      
      tmp_ext <- terra::extract(tmp_env, vect(tmp_seg))[,-1]
      tmp_ext <- tmp_ext[outer(date(tmp_seg$UTC), date(time(tmp_env)), \(x,y) x-y==0)]
      if(any(is.na(tmp_ext))) print("FYI: NAs in sd extraction!") 
      tmp_seg <- cbind(tmp_seg, tmp_ext)
      colnames(tmp_seg)[colnames(tmp_seg)=="tmp_ext"] <- vnms[j]
      # Focal SD versions of habitat variables 
      tmp_env_sd <- focal(tmp_env, 3, sd, na.rm=TRUE)
      tmp_ext <- terra::extract(tmp_env_sd, vect(tmp_seg))[,-1]
      tmp_ext <- tmp_ext[outer(date(tmp_seg$UTC), date(time(tmp_env)), \(x,y) x-y==0)]
      if(any(is.na(tmp_ext))) print("FYI: NAs in sd extraction!") 
      tmp_seg <- cbind(tmp_seg, tmp_ext)
      colnames(tmp_seg)[colnames(tmp_seg)=="tmp_ext"] <- paste0(vnms[j],"_sd")
    }
    
    tmp_seg
  }

#' -----------------------------------------------------------------------------
#' Clean up data 
#' -----------------------------------------------------------------------------
## Convert list of env data by year to dataframe, add Long and Lat columns, drop geometry column
segdata_env <- do.call(rbind, segdata_env)
segdata_env$Longitude <- st_coordinates(segdata_env)[,1]
segdata_env$Latitude <- st_coordinates(segdata_env)[,2]
segdata_env <- st_drop_geometry(segdata_env)

## Arrange segments by time (UTC)
segdata_env <- arrange(segdata_env, UTC)

## Reorder and rename covariates
cols_to_modify <- c("thetao", "thetao_sd", "so", "so_sd", "mlotst", "mlotst_sd", "zos", "zos_sd")

segdata_env <- segdata_env %>% select(-all_of(cols_to_modify), all_of(cols_to_modify) ) %>% 
  rename(sst = thetao, sst_sd=thetao_sd, salinity=so, salinity_sd=so_sd, mld=mlotst, mld_sd=mlotst_sd, ssh=zos, ssh_sd=zos_sd)

#' -----------------------------------------------------------------------------
#' Save outputs 
#' -----------------------------------------------------------------------------
write.csv(segdata_env, file = file.path(localwd,"output","segdata_env.csv"), row.names = FALSE)
saveRDS(segdata_env, file= file.path(localwd,"output","segdata_env.rds") )

