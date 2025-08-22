library(mgcv)
library(here)
library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(picMaps)
library(mapview)

local_wd <- here("analysis_tasks","task_10_VarProp")
setwd(local_wd)

#' -----------------------------------------------------------------------------
#' Load model output
#' -----------------------------------------------------------------------------
varprop_list <- readRDS(here( "analysis_tasks","task_2_g0_ESW","output","varprop_list.rds"))
load(here( "analysis_tasks","task_4_models","output","sdm_model_results.RData"))
load(file.path("output","varprop_results.RData"))

#' -----------------------------------------------------------------------------
#' Locate environmental data for prediction
#' -----------------------------------------------------------------------------
env_dir <- here("analysis_tasks","task_4_models","2017_test")
nc_files <- list.files(path = env_dir, pattern = "\\.nc$", full.names = TRUE)

#' -----------------------------------------------------------------------------
#' Raster and storage tools
#' -----------------------------------------------------------------------------
reps <- nrow(bs_gs)

N_storage <- expand.grid(rep=1:reps, day=seq(15,31,3), month=7:12, year=c(2020:2024))
N_storage <- N_storage %>% mutate(
  date = paste(month,day,year,sep="/") %>% mdy(),
  N_cenpac = as.numeric(NA),
  N_HIeez = as.numeric(NA)
) %>% select(-month, -day, -year)

N_storage <- N_storage %>% group_nest(date)

env_tmp <- rast(nc_files[1])[[1]]

cp_area <- cellSize(env_tmp, unit="km")
eez <- picMaps::hawaii_eez()
eez_r <- rasterize(vect(eez), env_tmp, cover=T) 


#' ----------------------------------------------------------------------------
#' Loop over years and parameters to make predictions
#' ----------------------------------------------------------------------------

for(i in 1:nrow(N_storage)){
  # env_file <- paste0()
  env <- rast(nc_files[1])
  # env <- env[[ as.Date(time(env))==N_storage$date[i] ]]
  env <- 1.0*env[[time(env)==ymd("2017-07-30")]]
  names(env) <- c("mld","sal","sst","u","v","ssh")
  env$sst_sd <- focal(env$sst, 3, sd, na.rm=TRUE)
  env_df <- as.data.frame(env, xy=TRUE, cells=TRUE)
  env_df$Latitude <- env_df$y
  env_df[["Jmat"]] <- matrix(0, nrow(env_df), 5)
  # L_er <- predict(ergam.HIeez, newdata=env_df, type="lpmatrix")
  L_er_vp <-  predict(ergam_vp, newdata=env_df, type="lpmatrix")
  L_gs <- predict(gsgam.HIeez.log, newdata=env_df, type="lpmatrix")
  
  abund_r <- abund_r_HIeez <- rast(env, nlyr=reps)
  N_HIeez_orig <- rep(NA,reps)
  
  for(j in 1:nrow(N_storage$data[[i]])){
    pred_val <- as.vector(exp(L_er_vp %*% bs_er[j,])) * as.vector(exp(L_gs %*% bs_gs[j,]))
    # pred_val_HIeez <- as.vector(exp(L_er %*% bs_er_HIeez[j,])) * as.vector(exp(L_gs %*% bs_gs[j,]))
    abund_r[[j]][env_df$cell] <- pred_val
    abund_r[[j]] <- abund_r[[j]] * cp_area
    # abund_r_HIeez[[j]][env_df$cell] <- pred_val_HIeez
    abund_r_HIeez[[j]] <- abund_r_HIeez[[j]] * cp_area
    N_storage$data[[i]]$N_cenpac[j] <- as.numeric(global(abund_r[[j]], "sum", na.rm=TRUE))
    N_storage$data[[i]]$N_HIeez[j] <- as.numeric(global(abund_r[[j]], "sum", weights=eez_r, na.rm=TRUE))
    N_HIeez_orig[j] <- as.numeric(global(abund_r_HIeez[[j]], "sum", weights=eez_r, na.rm=TRUE))
  }
  
}



