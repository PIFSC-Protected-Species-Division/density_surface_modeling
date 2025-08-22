
library(mgcv)
library(here)
library(tidyverse)
library(foreach)
library(future)
library(doFuture)
library(progressr)

#' -----------------------------------------------------------------------------
#' Load and process modeling data
#' -----------------------------------------------------------------------------

local_wd <- file.path(here(), "_R_code", "task_4_fit_dsm_models")
fkw_dsm_data <- readRDS(file.path(local_wd, "output", "fkw_dsm_data_1997to2024.rds"))

#' Add effort data
fkw_dsm_data <- fkw_dsm_data %>% mutate(
  effort = g0_033 * 2 * esw_033 * dist
)

#' Filter out years without any sightings
nz_sight <- fkw_dsm_data %>% group_by(year) %>% summarize(sightings = sum(nSI_033)) %>% 
  filter(sightings>0)

fkw_dsm_data <- fkw_dsm_data %>% filter(year%in%nz_sight$year, beaufort<7, effort>0)
fkw_gs <- fkw_dsm_data %>% filter(nSI_033 > 0)# Just segments with sightings

#' Filter out cruise 1004
# fkw_dsm_data <- fkw_dsm_data %>% filter(Cruise==1004)

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

#' -----------------------------------------------------------------------------
#' Loop over all model combinations
#' -----------------------------------------------------------------------------

#' Select spline type
#' 'ts'= thin plate splines ("shrinkage approach" applies additional smoothing penalty)
#' "tp" is default for s()
spline2use <- "ts" 

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
                     family = tw(),
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

#' -----------------------------------------------------------------------------
#' Group size model selection
#' -----------------------------------------------------------------------------

hab_inc_gs <- expand.grid(rep(list(0:1), length(hab_var)))
forms_gs <- apply(hab_inc_gs, 1, function(row) {
  terms <- mapply(function(choice, var) {
    if (choice == 1) return(paste0("s(",var,",bs=spline2use,k=5)"))
    return(NULL)
  }, row, hab_var)
  terms <- unlist(terms)
  # if (length(terms) == 0) return(NULL)  # skip null model
  paste(terms, collapse = " + ")
}) %>% unlist()
forms_gs <- paste("log(ANI_033)", forms_gs, sep=" ~ ")
forms_gs[1] <- "log(ANI_033) ~ 1"

gsgam_list <- foreach(
  i = seq_along(forms_gs), .options.future = list(seed = TRUE), 
  .errorhandling = "pass") %do% {
    spline2use <- spline2use
    gam_fit <- gam(formula = as.formula(forms_gs[i]),
                   method="REML", 
                   data = fkw_gs, weights = ProbPel)
    list(summary = summary(gam_fit), aic=gam_fit$aic)
  }

mod_df_gs <- data.frame(
  idx = 1:length(gsgam_list),
  form = forms_gs, 
  aic = sapply(gsgam_list, \(x) x$aic), 
  dev_expl=sapply(gsgam_list, \(x) x$summary$dev.expl*100)
) %>% mutate(
  waic = exp(-0.5*(aic-min(aic))),
  waic = waic/sum(waic)
) %>% arrange(desc(waic))

#' -----------------------------------------------------------------------------
#' Save model selection output
#' -----------------------------------------------------------------------------

write.csv(mod_df, file = file.path(local_wd, "output", "aic_er.csv"))
write.csv(mod_df_gs, file = file.path(local_wd, "output" "aic_gs.csv"))
save( fkw_dsm_data, fkw_gs, spline2use, file=file.path(local_wd, "output","dsm_model_data.RData"))





