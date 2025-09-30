library(mgcv)
library(here)
library(tidyverse)
library(sf)
library(picMaps)
library(units)
library(WeightedROC)
library(gratia)

#' -----------------------------------------------------------------------------
#' Load data and AIC tables
#' -----------------------------------------------------------------------------

local_wd <- here("_R_code", "task_4_fit_dsm_models")
load(file.path(local_wd, "output","dsm_model_data.RData"))
aic_er <- read_csv(file.path(local_wd, "output","aic_er.csv"))
aic_gs <- read_csv(file.path(local_wd, "output","aic_gs.csv"))

#' -----------------------------------------------------------------------------
#' Fit models
#' -----------------------------------------------------------------------------

# Use top AIC model minus non-significant terms

er_fit <- gam(formula = nSI_033 ~ te(sst, Latitude, bs = spline2use, k = 5),
              offset = log(effort),
              family = tw(),
              method="REML", 
              data = fkw_dsm_data, weights = ProbPel)
gs_fit <- gam(formula = log(ANI_033) ~ s(ssh, bs = spline2use, k = 5) + s(mld, bs = spline2use, k = 5),
              method="REML", 
              data = fkw_gs, weights = ProbPel)

#' -----------------------------------------------------------------------------
#' Evaluate models
#' -----------------------------------------------------------------------------

### CenPac area
A <- picMaps::cenpac() %>% st_area() %>% set_units("km^2") %>% as.numeric()

### model predictions
pred_er <- exp(predict(er_fit) + log(fkw_dsm_data$effort))
pred_gs <- exp(predict(gs_fit) + 0.5*gs_fit$sig2)
gs_bias_corr <- weighted.mean(fkw_gs$ANI_033,w=fkw_gs$ProbPel)/mean(pred_gs)
pred_gs <- exp(predict(gs_fit, newdata=fkw_dsm_data) + 0.5*gs_fit$sig2)

### Design-based observed and predicted
obsdens <- (fkw_dsm_data$nSI_033 * fkw_dsm_data$ANI_033)
ltdens <- weighted.mean(obsdens,fkw_dsm_data$ProbPel)/weighted.mean(fkw_dsm_data$effort,fkw_dsm_data$ProbPel)
ltdens*A
# 23441.83

## model based
modeldens <- pred_er * pred_gs * gs_bias_corr
preddens <- sum(modeldens)/sum(fkw_dsm_data$effort)
preddens*A
# 22154.7

# Ratio
ltdens/preddens 
# 1.058097

### Compute AUC & TSS 

#  Set up presence and absence points as 0's and 1's
sightings <- as.numeric(fkw_dsm_data$nSI_033>=1)
pred_seg <- predict(er_fit, type="response") * pred_gs * gs_bias_corr 

rocr <- WeightedROC(pred_seg, sightings, weight = fkw_dsm_data$ProbPel)

## AUC
auc <- WeightedAUC(rocr)
auc
# 0.7119164

#  Plot ROC curve
ggplot()+geom_path(aes(FPR, TPR), data=rocr)+coord_equal()


## Compute TSS
d <- sqrt((rocr$TPR-1)^2 + (rocr$FPR-0)^2)
which(d==min(d))
opt_rocr <- rocr[which(d==min(d)),]

TSS <- opt_rocr$TPR + (1-opt_rocr$FPR) -1 
TSS
# 0.3626605

#' -----------------------------------------------------------------------------
#' ROSE Bootstrap
#' -----------------------------------------------------------------------------

rose_df <- fkw_dsm_data %>% select(
  nSI_033, sst, sst_sd, salinity, salinity_sd, mld, mld_sd, ssh, ssh_sd, Latitude, Longitude, effort
) %>% rename(cls = nSI_033)

auc_stor <- rep(NA, 50)
tss_stor <- rep(NA, 50)
for(i in 1:50){
  resamp <- ROSE(cls~., data=rose_df, hmult.majo=0.25, hmult.mino=0.25, seed=i)$data
  er_fit_rose <- gam(formula = cls ~ te(sst, Latitude, bs = spline2use, k = 5),
                     family = tw(),
                     method="REML", 
                     data = resamp)
  mu <- predict(er_fit_rose, newdata = fkw_dsm_data, type="response")
  phi <- er_fit_rose$sig2
  theta <- er_fit_rose$family$getTheta()
  p <- (1+2*exp(theta))/(1+exp(theta)); p <- max(1.01, p); p <- min(1.99, p)
  pnz <- 1 - exp(- (mu^(2-p)/(phi*(2-p))))
  rrr <-  WeightedROC(pnz, sightings, weight = fkw_dsm_data$ProbPel)
  auc_stor[i] <- WeightedAUC(rrr)
  tss_stor[i] <- max(rrr$TPR + (1-rrr$FPR) -1, na.rm=T)
  cat(i," ")
}
mean(auc_stor)
# 0.7626317
range(auc_stor)
# 0.7581797 0.7666302

mean(tss_stor)
# 0.4331774
range(tss_stor)
# 0.4144445 0.4537293


#' -----------------------------------------------------------------------------
#' Save results
#' -----------------------------------------------------------------------------

#' plots
sm_grid <- smooth_estimates(er_fit, select = "te(sst,Latitude)")
out <- exclude.too.far(sm_grid$sst, sm_grid$Latitude, fkw_dsm_data$sst, fkw_dsm_data$Latitude, 0.1)
sm_grid$.estimate[out] <- NA
ggplot(data = sm_grid) +
  geom_raster(aes(x = sst, y = Latitude, fill = .estimate)) +
  geom_contour(aes(x = sst, y = Latitude, z = .estimate), colour = "black") +
  geom_point(
    data = dat_res,
    aes(x = sst, y = Latitude),
    size = 0.25, 
    alpha = 0.25 
  ) +
  labs(title = "te(sst, Latitude)") + 
  scale_fill_distiller(palette = "RdBu", type = "div", name = "effect", na.value = NA) +
  coord_cartesian(expand = FALSE) +  theme_bw()
ggsave(file=file.path(local_wd,"output","er_fit.png"),width=6.5, height=4)

draw(gs_fit) & theme_bw()
ggsave(file=file.path(local_wd,"output","gs_fit.png"),width=6.5, height=4)

save(er_fit, gs_fit, gs_bias_corr, opt_rocr, fkw_dsm_data, fkw_gs, A, ltdens, spline2use,
     file=file.path(local_wd, "output", "dsm_fit_final.RData"))

     