
library(mgcv)
library(here)
library(tidyverse)
library(foreach)
library(future)
library(doFuture)
library(progressr)

#  For AUC calcs
library(ROCR)

#' -----------------------------------------------------------------------------
#' Load and process modeling data
#' -----------------------------------------------------------------------------

local_wd <- file.path(here(), "_R_code", "task_4_fit_dsm_models")
fkw_dsm_data <- readRDS(file.path(local_wd, "output", "fkw_dsm_data_1997to2024.rds"))

#' Extract species data for modeling
fkw_dsm_data <- fkw_dsm_data %>% mutate(
  effort = g0_033 * 2 * esw_033 * dist
)

#' Filter out years without any sightings
nz_sight <- fkw_dsm_data %>% group_by(year) %>% summarize(sightings = sum(nSI_033)) %>% 
  filter(sightings>0)

fkw_dsm_data <- fkw_dsm_data %>% filter(year%in%nz_sight$year, beaufort<7, effort>0)
fkw_gs <- fkw_dsm_data %>% filter(nSI_033 > 0)# Just segments with sightings

#' -----------------------------------------------------------------------------
#' Add augmented data for pelagic probabilty weighting
#' -----------------------------------------------------------------------------
nonpel <- fkw_dsm_data %>% filter(ProbPel<1)
nonpel <- nonpel %>% mutate(
  nSI_033 = 0,
  ProbPel = 1-ProbPel,
  augment=1
)

fkw_dsm_data <- bind_rows(fkw_dsm_data, nonpel) %>% 
  mutate(
    augment = ifelse(is.na(augment), 0, 1)
  )

#' -----------------------------------------------------------------------------
#' Extract sample size summaries
#' -----------------------------------------------------------------------------

# Need to modify with new dplyr tools!
# yearly.ss <- HIeez %>% group_by(year) %>%
#   summarize(
#     sightings = sum(.data[[nSI.col]])
#     )

#' -----------------------------------------------------------------------------
#' Create models to select from 
#' -----------------------------------------------------------------------------

hab_var <- c("sst","ssh","mld","salinity")
hab_var_sd <- c("sst_sd","ssh_sd","mld_sd","salinity_sd")

hab_inc <- expand.grid(rep(list(0:2), length(hab_var)))
hab_inc <- hab_inc[apply(hab_inc==2, 1, sum)<2,]
hab_sd_inc <- expand.grid(rep(list(0:1), length(hab_var_sd)))

form_hab <- apply(hab_inc, 1, function(row) {
  terms <- mapply(function(choice, var) {
    if (choice == 1) return(paste0("s(",var,",bs=spline2use,k=5)"))
    if (choice == 2) return(paste0("te(", var, ",Latitude,bs=spline2use,k=5)"))
    return(NULL)
  }, row, hab_var)
  terms <- unlist(terms)
  if (length(terms) == 0) return("nSI_033 ~ ")  # skip null model
  paste("nSI_033 ~", paste(terms, collapse = " + "))
}) %>% unlist()

form_hab_sd <- apply(hab_sd_inc, 1, function(row) {
  terms <- mapply(function(choice, var) {
    if (choice == 1) return(paste0("s(",var,",bs=spline2use,k=5)"))
    return(NULL)
  }, row, hab_var_sd)
  terms <- unlist(terms)
  # if (length(terms) == 0) return(NULL)  # skip null model
  paste(terms, collapse = " + ")
}) %>% unlist()

forms_er <- NULL
for(i in seq_along(form_hab)){
  if(i ==1) {
    tmp <- paste(form_hab[i], form_hab_sd[-1], sep="")
  } else {
    tmp <- c(form_hab[i], paste(form_hab[i], form_hab_sd[-1], sep=" + "))
  }
  forms_er <- c(forms_er, tmp)
}

# form_hab_sd <- paste(paste0("s(", hab_var_sd, ",bs=spline2use,k=5)"), collapse = " + ")
# forms_er <- paste(form_hab, form_hab_sd, sep=" + ")

#' -----------------------------------------------------------------------------
#' Loop over all model combinations
#' -----------------------------------------------------------------------------

#' Select spline type
#' 'ts'= thin plate splines ("shrinkage approach" applies additional smoothing penalty)
#' "tp" is default for s()
spline2use <- "ts" 
gam_select <- FALSE
# if(spline2use == "tp"){
#   gam_select <- TRUE
# } else gam_select <- FALSE

workers <- 8
plan("multisession", workers=workers)
# forms <- split(forms, cut(seq_along(forms), breaks = workers, labels = FALSE))

with_progress({
  pb <- progressor(along = forms_er) 
  ergam_list <- foreach(
    i = seq_along(forms_er), .options.future = list(seed = TRUE), 
    .errorhandling = "pass") %dofuture% {
      spline2use <- spline2use
      gam_fit <- gam(formula = as.formula(forms_er[i]), offset = log(effort),
                     family = tw(), #select=gam_select,
                     method="REML", 
                     data = fkw_dsm_data, weights = ProbPel)
      pb()
      list(summary = summary(gam_fit), aic=gam_fit$aic)
    }
})
plan("sequential")

mod_df <- data.frame(
  idx = 1:length(ergam_list),
  form = forms_er, 
  aic = sapply(ergam_list, \(x) x$aic), 
  dev_expl=sapply(ergam_list, \(x) x$summary$dev.expl*100)
) %>% mutate(
  waic = exp(-0.5*(aic-min(aic))),
  waic = waic/sum(waic)
) %>% arrange(desc(waic))

top_aic <- mod_df[1:10,]
mod_df <- mod_df %>% arrange(desc(dev_expl))
top_dev <- mod_df[1:10,]

write.csv(top_aic, file = file.path(local_wd,"top_aic.csv"))
write.csv(top_dev, file = file.path(local_wd,"top_dev.csv"))




old_ergam <- gam(nSI_033 ~ s(sst, bs="ts") + s(sst_sd, bs="ts") + te(ssh, Latitude, bs="ts", k=5)
             , offset = log(effort),
             family = tw(), #select=TRUE,
             method="REML", 
             data = fkw_dsm_data, weights = ProbPel)


#' Group Size GAM

jchxvb

form_hab2 <- paste(paste0("log(ANI_033) ~ s(", hab_var, ",bs=spline2use,k=5)"), collapse = " + ")


gsgam <- gam(formula = log(ANI_033) ~
               s(mld, bs=spline2use)
             + s(ssh, bs=spline2use)
             , family = gaussian,
             method="REML", 
             data =  fkw_gs, weights = ProbPel)


##############################################################
## MODEL SUMMARIES AND PLOTS  

summary(ergam.HIeez)
summary(gsgam.HIeez.log)

rqgam.check (ergam) # Function operates as gam.check but uses randomized quantile residuals
gam.check(gsgam.HIeez.log)

sqrt(mean(residuals.gam(ergam.HIeez,type="response")^2)) #Model error as root mean squared error (RMSE) 
sqrt(mean(residuals.gam(gsgam.HIeez.log,type="response")^2)) #Model error as root mean squared error (RMSE) 


# plot(ergam.HIeez,scale=0,residuals=TRUE)
# par(mfrow=c(2,2))
plot(ergam,scale=0, shade=TRUE, scheme=2, contour.col = 1, rug=TRUE)

# plot(gsgam.HIeez.log,scale=0,residuals=TRUE)
# par(mfrow=c(1,1))
plot(gsgam,scale=0,scheme=2, shade=TRUE)

save(HIeez, HIeez.no0, ergam.HIeez, gsgam.HIeez.log, spline2use, file=file.path("output","sdm_model_results.RData"))


#' -----------------------------------------------------------------------------
#' Attempt new predictions with terra::predict
#' -----------------------------------------------------------------------------

# eez <- picMaps::hawaii_eez()
# env <- rast(file.path(local_wd, "2017_test", "cenpac_pred_grid_2017.nc"))
# lat_r <- rast(env, nlyr=1)
# values(lat_r) <- crds(lat_r)[,"y"]
# pt <- read_xlsx(file.path(local_wd, "2017_test","Original_SDM_Tri_Daily_Pred_Dates.xlsx"), col_names = FALSE) %>% 
#   pull(2)
# pt <- pt[year(pt)==2017]
# 
# 
# cp_abund_r <- eez_abund_r <- rast(env, nlyrs=length(pt))
# time(cp_abund_r) <- pt
# time(eez_abund_r) <- pt
# cp_area <- cellSize(env, unit="km")
# eez_r <- rasterize(vect(eez), env, cover=T) 
# cp_abund <- eez_abund <- rep(NA, length(pt))
# 
# 
# 
# for(i in seq_along(pt)){
#   env_tmp <- env[[time(env)==pt[i]]]
#   names(env_tmp) <- c("mld","sal","sst","u","v","ssh")
#   env_tmp$sst_sd <- focal(env_tmp$sst, 3, sd, na.rm=TRUE)
#   env_tmp$ssh_sd <- focal(env_tmp$ssh, 3, sd, na.rm=TRUE)
#   env_tmp$sal_sd <- focal(env_tmp$sal, 3, sd, na.rm=TRUE)
#   env_tmp$Latitude <- lat_r
#   pred_er <- terra::predict(env_tmp, ergam.HIeez, type="response")
#   pred_gs <- terra::predict(env_tmp, gsgam.HIeez.log, type="response") %>% exp()
#   pred <- pred_er * pred_gs
#   cp_abund_r[[i]] <- cp_area * pred
#   cp_abund[i] <- as.numeric(global(pred, "sum", weights=cp_area, na.rm=TRUE))
#   eez_abund_r[[i]] <- eez_r * cp_abund_r[[i]]
#   eez_abund[i] <- as.numeric(global(pred, "sum", cp_area*eez_r, na.rm=TRUE)) #sum(values(eez_abund_r[[i]]), na.rm=T)
#   cat(i," ")
# }
# 
# cp_abund_r <- wrap(cp_abund_r)
# eez_abund_r <- wrap(eez_abund_r)
# 
# save(list=ls(), file="saved.RData")




