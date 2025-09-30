library(mgcv)
library(here)
library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(picMaps)
library(mapview)
library(foreach)
library(doFuture)

local_wd <- here("_R_code","task_8_propagating_variance")

#' -----------------------------------------------------------------------------
#' Load model output
#' -----------------------------------------------------------------------------
varprop_list <- readRDS(here( "_R_code","task_2_est_g0_ESW","output","varprop_list.rds"))
load(here( "_R_code","task_4_fit_dsm_models","output","dsm_fit_final.RData"))
load(file.path(local_wd, "output","varprop_fit_sample.RData"))
load(here( "_R_code","task_7_predicting_abundance","output","clamp_df.RData"))

#' -----------------------------------------------------------------------------
#' Locate environmental data for prediction
#' -----------------------------------------------------------------------------
env_dir <- here("_R_code","task_5_download_prediction_data","output")
tif_files <- list.files(path = env_dir, pattern = "\\.tif$", full.names = TRUE)

# tif_files <- tif_files[c(1,5)]

#' -----------------------------------------------------------------------------
#' Raster and storage tools
#' -----------------------------------------------------------------------------
reps <- nrow(bs_gs)

env_tmp <- rast(tif_files[1])[[1]]
area <- rast(ext(env_tmp), resolution=res(env_tmp)); crs(area) <- crs(env_tmp)
area <- cellSize(area, unit="km")
hi_r <-  picMaps::hawaii_coast() %>% rasterize(area, touches=TRUE, cover=TRUE, background=0)
hi_r <- 1-hi_r
hi_r <- ifel(hi_r==min(values(hi_r)), 0, hi_r)
cenpac_r <-  picMaps::cenpac() %>% rasterize(area, touches=TRUE, cover=TRUE) * hi_r
names(cenpac_r) <- "cenpac_wt"
eez_r <- picMaps::hawaii_eez() %>% rasterize(area, touches=TRUE, cover=TRUE) * hi_r
names(eez_r) <- "hi_eez_wt"
fkw_mgmt_r <- picMaps::pfkw_mgmt() %>% rasterize(area, touches=TRUE, cover=TRUE) * hi_r
names(fkw_mgmt_r) <- "fkw_assess_wt"
pred_mult_r <- c(area, cenpac_r, eez_r, fkw_mgmt_r)
pred_mult_df <- as.data.frame(pred_mult_r, xy=TRUE, cells=TRUE)

L_gs_data <- predict(gs_fit, type="lpmatrix")
obs_gs <- weighted.mean(fkw_gs$ANI_033,w=fkw_gs$ProbPel)

#' -----------------------------------------------------------------------------
#' Model covariate names
#' -----------------------------------------------------------------------------
er_vars <- all.vars(as.list(er_fit$call)$formula[[3]])
er_vars <- er_vars[er_vars!="spline2use"]
gs_vars <- all.vars(as.list(gs_fit$call)$formula[[3]])
gs_vars <- gs_vars[gs_vars!="spline2use"]
var_nms <- c(er_vars, gs_vars)

clamp_df <- filter(clamp_df, var%in%var_nms, var!="Latitude")

#' -----------------------------------------------------------------------------

#' ----------------------------------------------------------------------------
#' Loop over years and parameters to make predictions
#' ----------------------------------------------------------------------------

plan("multisession",workers=length(tif_files))

N_storage <- foreach(i = 1:length(tif_files), .options.future=list(seed = TRUE)) %dofuture% { #i=year
  env_yr <- rast(tif_files[i])
  env_yr <- env_yr[[names(env_yr)%in%var_nms]]
  dates <- time(env_yr) %>% unique()
  year_i <- year(dates)[1]
  pred_est_r <- rast(ext(env_yr), resolution=res(env_yr), nlyr=length(dates)); crs(pred_est_r) <- crs(env_yr)
  pred_var_r <- pred_est_r
  N_stor <- expand.grid(par_rep=1:reps, date=dates) %>% mutate(
    N_cenpac = as.numeric(NA),
    N_eez = as.numeric(NA),
    N_assess = as.numeric(NA)
  ) %>% group_nest(date)
  for(j in 1:length(dates)){ # j=day
    env <- env_yr[[time(env_yr)==dates[j]]]
    env_df <- as.data.frame(env, xy=TRUE, cells=TRUE, na.rm=FALSE)
    pred_r <- rast(ext(env_yr), resolution=res(env_yr), nlyr=reps); crs(pred_r) <- crs(env_yr)
    for(l in 1:nrow(clamp_df)){
      env_df[clamp_df$var[l]] <- ifelse( env_df[[clamp_df$var[l]]]>clamp_df[l,"high"], clamp_df[l,"high"], env_df[[clamp_df$var[l]]])
      env_df[clamp_df$var[l]] <- ifelse( env_df[[clamp_df$var[l]]]<clamp_df[l,"low"], clamp_df[l,"low"], env_df[[clamp_df$var[l]]])
    }
    env_df$Latitude <- env_df$y
    env_df[["Jmat"]] <- matrix(0, nrow(env_df), 5)
    L_er_vp <-  predict(er_fit_vp, newdata=env_df, type="lpmatrix")
    L_gs <- predict(gs_fit, newdata=env_df, type="lpmatrix")
    for(k in 1:reps){ # k=par_rep
      pred_gs_bias <- exp(L_gs_data%*%bs_gs[k,])
      gs_bias_corr <- obs_gs/mean(pred_gs_bias)
      pred_k <-  as.vector(exp(L_er_vp %*% bs_er[k,]) * exp(L_gs %*% bs_gs[k,])) * gs_bias_corr
      pred_k <- pred_k * pred_mult_df$area
      pred_r[[k]][env_df$cell] <- pred_k
      N_stor$data[[j]]$N_cenpac[k] <-  sum(pred_k* pred_mult_df$cenpac_wt, na.rm=T)
      N_stor$data[[j]]$N_eez[k] <-  sum(pred_k* pred_mult_df$hi_eez_wt, na.rm=T)
      N_stor$data[[j]]$N_assess[k] <- sum(pred_k* pred_mult_df$fkw_assess_wt, na.rm=T)
    } # end k
    pred_est_r[[j]] <- app(pred_r, mean)
    pred_var_r[[j]] <- app(pred_r, var)
  } # end j
time(pred_est_r) <- dates
est_outfile <-  paste0("pred_est_",year_i,".tif")
var_outfile <- paste0("pred_var_",year_i,".tif")
writeRaster(pred_est_r, file.path(local_wd,"output","predictions", est_outfile), overwrite=TRUE)
writeRaster(pred_var_r, file.path(local_wd,"output","predictions", var_outfile), overwrite=TRUE)

tibble(year=year_i, data=list(N_stor), est_file=est_outfile, var_file=var_outfile)
} # end i
plan("sequential")

N_storage <- do.call(bind_rows, N_storage)

N_storage <- N_storage %>% unnest(cols=data) %>% unnest(cols=data)
write_csv(N_storage, file=file.path(local_wd,"output","N_varprop.csv"))

#' Calculate CV
N_cv <- N_storage %>% group_by(year) %>% 
  summarise(
    N_cenpac_m=mean(N_cenpac),  N_cenpac_sd=sd(N_cenpac),
    N_eez_m=mean(N_eez),  N_eez_sd=sd(N_eez),
    N_assess_m=mean(N_assess),  N_assess_sd=sd(N_assess),
  ) %>% mutate(
    N_cenpac_cv = N_cenpac_sd/N_cenpac_m,
    N_eez_cv = N_eez_sd/N_eez_m,
    N_assess_cv = N_assess_sd/N_assess_m
  ) %>% mutate(year = as.character(year))

N_cv_avg <- filter(N_storage, year!=2017) %>% 
  summarise(
    N_cenpac_m=mean(N_cenpac),  N_cenpac_sd=sd(N_cenpac),
    N_eez_m=mean(N_eez),  N_eez_sd=sd(N_eez),
    N_assess_m=mean(N_assess),  N_assess_sd=sd(N_assess),
  ) %>% mutate(
    N_cenpac_cv = N_cenpac_sd/N_cenpac_m,
    N_eez_cv = N_eez_sd/N_eez_m,
    N_assess_cv = N_assess_sd/N_assess_m
  ) %>% mutate(year="avg 2020-2024") %>% 
  relocate(year, .before=1)
N_cv <- bind_rows(N_cv_avg, N_cv)


N_cv <- N_cv %>% select(year, N_cenpac_cv,  N_eez_cv,  N_assess_cv) %>% 
  mutate(
    sigma_cenpac = sqrt(log(1+N_cenpac_cv^2)),
    sigma_eez = sqrt(log(1+N_eez_cv^2)),
    sigma_assess = sqrt(log(1+N_assess_cv^2))
  )
N_summ <- read_csv(here("_R_code","task_7_predicting_abundance","output","N_summ.csv"))
N_cv <- N_cv %>% full_join(N_summ)
N_cv <- N_cv %>% mutate(
  mu_cenpac = log(N_cenpac_est)-sigma_cenpac^2/2,
  mu_eez = log(N_eez_est)-sigma_eez^2/2,
  mu_assess = log(N_fkw_mgmt_est)-sigma_assess^2/2,
  CI_low_cenpac = exp(mu_cenpac - 1.96*sigma_cenpac),
  CI_hi_cenpac = exp(mu_cenpac + 1.96*sigma_cenpac),
  CI_low_eez = exp(mu_eez - 1.96*sigma_eez),
  CI_hi_eez = exp(mu_eez + 1.96*sigma_eez),
  CI_low_assess = exp(mu_assess - 1.96*sigma_assess),
  CI_hi_assess = exp(mu_assess + 1.96*sigma_assess)
)
N_cv <- N_cv %>% select(year, 
                        N_cenpac_est, CI_low_cenpac, CI_hi_cenpac, N_cenpac_cv,
                        N_eez_est, CI_low_eez, CI_hi_eez, N_eez_cv,
                        N_fkw_mgmt_est, CI_low_assess, CI_hi_assess, N_assess_cv)
write_csv(N_cv, file=file.path(local_wd,"output","N_cv_varprop.csv"))

