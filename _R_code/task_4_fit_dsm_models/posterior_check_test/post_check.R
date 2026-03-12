library(mgcv)
library(here)
library(tidyverse)
library(sf)
library(picMaps)
library(units)
library(WeightedROC)
library(gratia)
library(ROSE)

#' -----------------------------------------------------------------------------
#' Load data and AIC tables
#' -----------------------------------------------------------------------------

local_wd <- here("_R_code", "task_4_fit_dsm_models")
load(file.path(local_wd, "output","dsm_model_data.RData"))

load(here("_R_code","task_4_fit_dsm_models","output","dsm_fit_final.RData"))

supp_files <- list.files(path = here("_R_code","task_8_propagating_variance","support_scripts"), pattern = "\\.R$", full.names = TRUE)
for(i in seq_along(supp_files)) source(supp_files[i])

draw_tw <- function(mu, p, phi, ...){
  nc <- ifelse(is.null(ncol(mu)), 1, ncol(mu))
  nr <- ifelse(is.null(nrow(mu)), length(mu), nrow(mu))
  o <- sapply(1:nr,
              function(i, mu, p, phi){mgcv::rTweedie(mu[i,],p,phi)},
              mu=mu, p=p, phi=phi
  )
  return(as.vector(o))
}


beta_er_samp <- gam.imp(er_fit, 10000)
beta_er_samp$wts <- beta_er_samp$wts/sum(beta_er_samp$wts)
# 1/sum(beta_er_samp$wts^2) # Check effective sample size
idx <- sample(1:nrow(beta_er_samp$bs), 199, TRUE, prob = beta_er_samp$wts)
b_er <- rbind(er_fit$coefficients, beta_er_samp$bs[idx,])

beta_gs_samp <- gam.imp(gs_fit, 10000)
beta_gs_samp$wts <- beta_gs_samp$wts/sum(beta_gs_samp$wts)
# 1/sum(beta_gs_samp$wts^2) # Check effective sample size
idx <- sample(1:nrow(beta_gs_samp$bs), 199, TRUE, prob = beta_gs_samp$wts)
b_gs <- rbind(gs_fit$coefficients, beta_gs_samp$bs[idx,])



# Simulate data
p <- er_fit$family$getTheta(TRUE)
phi <- summary(er_fit)$dispersion
Ler <- predict(er_fit, newdata=fkw_dsm_data, type="lpmatrix")
Lgs <- predict(gs_fit, newdata=fkw_dsm_data, type="lpmatrix")

nz <- NULL
dens <- NULL
for(i in 1:nrow(b_er)){
  pred_er <- as.vector(exp(Ler %*% b_er[i,]+ log(fkw_dsm_data$effort)))
  pred_gs <- as.vector(exp(Lgs%*%b_gs[i,] + 0.5*gs_fit$sig2))
  modeldens <- pred_er * pred_gs * gs_bias_corr
  preddens <- sum(modeldens)/sum(fkw_dsm_data$effort)
  dens <- c(dens,preddens*A)
  tt <- draw_tw(mu, p, phi)
  nz <- c(nz, sum(tt>0))
}






