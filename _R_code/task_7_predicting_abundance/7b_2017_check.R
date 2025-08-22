library(mgcv)
library(here)
library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(picMaps)
library(mapview)
library(doFuture)
library(foreach)

local_wd <- here("_R_code","task_7_predicting_abundance")
dir.create(file.path(local_wd, "output"))
dir.create(file.path(local_wd, "output","predictions"))

env_dir <- here("_R_code","task_5_download_prediction_data","output")

#' -----------------------------------------------------------------------------
#' Load model output
#' -----------------------------------------------------------------------------

load(here("_R_code","task_4_fit_dsm_models","output","dsm_fit_final.RData"))

#' -----------------------------------------------------------------------------
#' Refit GAMs without new data 2020-2023 and 2010, cruise 1004
#' -----------------------------------------------------------------------------

fkw_dsm_data <- fkw_dsm_data %>% filter(year<=2017, Cruise!=1004)
fkw_gs <- fkw_gs %>% filter(year<=2017, Cruise!=1004)

spline2use <- "ts"

er_fit <- gam(as.list(er_fit$call)$formula, 
              offset = log(effort),
              family = tw(),
              method="REML", 
              data = fkw_dsm_data, weights = ProbPel)
gs_fit <- gam(as.list(gs_fit$call)$formula,
              method="REML", 
              data = fkw_gs, weights = ProbPel)


#' -----------------------------------------------------------------------------
#' Raster and storage tools
#' -----------------------------------------------------------------------------

env_tmp <- rast(file.path(env_dir,"cenpac_2017_env.tif"))[[1]] 
area <- rast(ext(env_tmp), resolution=res(env_tmp)); crs(area) <- crs(env_tmp)
area <- cellSize(area, unit="km")
hi_r <-  picMaps::hawaii_coast() %>% rasterize(area, touches=TRUE, cover=TRUE, background=0)
hi_r <- 1-hi_r
hi_r <- ifel(hi_r==min(values(hi_r)), 0, hi_r)
cenpac_r <-  picMaps::cenpac() %>% rasterize(area, touches=TRUE, cover=TRUE)
cenpac_r <- cenpac_r * hi_r
eez_r <- picMaps::hawaii_eez() %>% rasterize(area, touches=TRUE, cover=TRUE)
fkw_mgmt_r <- picMaps::pfkw_mgmt() %>% rasterize(area, touches=TRUE, cover=TRUE)

#' -----------------------------------------------------------------------------
#' Determine clamping and truncating values 
#' -----------------------------------------------------------------------------

fit_var <- fkw_dsm_data %>%  select(sst:ssh_sd, Latitude)
clamp_df <- data.frame(var=colnames(fit_var), low=NA, high=NA, trunc_low=NA, trunc_high=NA)
for(i in 1:nrow(clamp_df)){
  v <- fkw_dsm_data[[clamp_df$var[i]]]
  clamp_df$low[i] <- quantile(v, 0.001)
  clamp_df$high[i] <- quantile(v, 0.999)
  clamp_df$trunc_low[i] <- 0.5*min(v)
  clamp_df$trunc_high[i] <- 1.5*max(v)
}


#' ----------------------------------------------------------------------------
#' Loop over years and parameters to make predictions
#' ----------------------------------------------------------------------------

# plan("multisession",workers=8)

N_data <- NULL

env <- rast(file.path(env_dir, "cenpac_2017_env.tif"))
times <- time(env) %>% unique() %>% sort()
for(j in seq_along(times)){
  sub_idx <- (terra::time(env)==times[j])
  env_j <- terra::subset(env, sub_idx)
  lat <- rast(env_j, nlyr=1); values(lat) <- terra::crds(lat)[,"y"]
  # clamp
  for(k in 1:nlyr(env_j)){
    clamp <- clamp_df[with(clamp_df, var==names(env_j)[k]),]
    env_j[[k]] <- ifel(env_j[[k]]<clamp$low, clamp$low, env_j[[k]])
    env_j[[k]] <- ifel(env_j[[k]]>clamp$high, clamp$high, env_j[[k]])
  }
  env_j$Latitude <- lat
  # Predict er and gs
  pred_er <- terra::predict(env_j, er_fit, type="response")
  pred_gs <- exp(terra::predict(env_j, gs_fit) + 0.5*gs_fit$sig2) * gs_bias_corr
  pred <- pred_er * pred_gs * area * cenpac_r
  pred <- ifel(is.na(pred), 0, pred)
  N_cenpac <- as.numeric(global(pred, "sum", na.rm=TRUE))
  N_eez <- as.numeric(global(pred, "sum", weights=eez_r, na.rm=TRUE))
  N_fkw <- as.numeric(global(pred, "sum", weights=fkw_mgmt_r, na.rm=TRUE))
  outfile <- paste0("N_",times[j],"_file.tif")
  N_data <- bind_rows(
    N_data,
    data.frame(year=2017, date=times[j], N_cenpac=N_cenpac, N_eez=N_eez, 
               N_fkw_mgmt=N_fkw, file=outfile
    )
  )
  writeRaster(pred, file.path(local_wd,"output","predictions","2017",outfile), overwrite=TRUE)
  cat(j," ")
}

save(N_data, clamp_df, file=file.path(local_wd, "output","2017","N_est_df.RData"))

N_summ <- N_data %>% group_by(year) %>% 
  summarize(
    N_cenpac_est = round(mean(N_cenpac)),
    N_cenpac_sd = round(sd(N_cenpac)),
    N_eez_est = round(mean(N_eez)),
    N_eez_sd = round(sd(N_eez)),
    N_fkw_mgmt_est = round(mean(N_fkw_mgmt)),
    N_fkw_mgmt_sd = round(sd(N_fkw_mgmt))
  )

write_csv(N_summ, file=file.path(local_wd,"output","2017","N_summ.csv"))
