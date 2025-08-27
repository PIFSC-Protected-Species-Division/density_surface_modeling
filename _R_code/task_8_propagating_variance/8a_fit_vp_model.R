
#' -----------------------------------------------------------------------------
#' load packages
#' -----------------------------------------------------------------------------

library(mgcv)
library(dsm)
library(numDeriv)
library(here)
library(tidyverse)
library(Distance)

local_wd <- here("_R_code","task_8_propagating_variance")

#' -----------------------------------------------------------------------------
#' load support scripts
#' -----------------------------------------------------------------------------\

supp_files <- list.files(path = file.path(local_wd,"support_scripts"), pattern = "\\.R$", full.names = TRUE)
for(i in seq_along(supp_files)) source(supp_files[i])

#' -----------------------------------------------------------------------------
#' Load previous data and model objects
#' -----------------------------------------------------------------------------

varprop_list <- readRDS(here( "_R_code","task_2_est_g0_ESW","output","varprop_list.rds"))
load(here("_R_code","task_4_fit_dsm_models","output","dsm_fit_final.RData"))



#' -----------------------------------------------------------------------------
#' Fit model with Jacobian based RE
#' -----------------------------------------------------------------------------

er_fit_vp <- add_varprop(er_fit,
                        beaufort_data = select(fkw_dsm_data, beaufort), 
                        varprop_list, trace=TRUE)

#' -----------------------------------------------------------------------------
#' Obtain importance sample of parameters
#' -----------------------------------------------------------------------------

set.seed(8675309)

beta_er_samp <- gam.imp(er_fit_vp, 10000)
beta_er_samp$wts <- beta_er_samp$wts/sum(beta_er_samp$wts)
# 1/sum(beta_er_samp$wts^2) # Check effective sample size
idx <- sample(1:nrow(beta_er_samp$bs), 199, TRUE, prob = beta_er_samp$wts)
bs_er <- rbind(er_fit_vp$coefficients, beta_er_samp$bs[idx,])


beta_gs_samp <- gam.imp(gs_fit, 10000)
beta_gs_samp$wts <- beta_gs_samp$wts/sum(beta_gs_samp$wts)
# 1/sum(beta_gs_samp$wts^2) # Check effective sample size
idx <- sample(1:nrow(beta_gs_samp$bs), 199, TRUE, prob = beta_gs_samp$wts)
bs_gs <- rbind(gs_fit$coefficients, beta_gs_samp$bs[idx,])

save(er_fit_vp, bs_er, bs_gs, file=file.path(local_wd, "output","varprop_fit_sample.RData"))
