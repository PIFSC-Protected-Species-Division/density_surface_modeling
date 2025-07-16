#' -----------------------------------------------------------------------------
#'  Prepare Data for SDM Models
#' -----------------------------------------------------------------------------

library(here)
library(tidyverse)


#' -----------------------------------------------------------------------------
#' Load Data
#' -----------------------------------------------------------------------------

localwd <- here("_R_code","task_4_fit_dsm_models")

segDataESWg0 <- readRDS( here("_R_code","task_2_est_g0_esw","output","seg_sight_out_g0_ESW.rds") )

#' Read table defining which sightings are from the pelagic population
fkw_popData <- readRDS( here("_R_code", "task_4_fit_dsm_models", "FKWsights_DSM_1997to2023.rds") ) 


#' Select sightings that are from pelagic 
spcode <- '033'
segDataESWg0$sightinfo <- segDataESWg0$sightinfo %>% filter(SpCode==spcode)%>%
  dplyr::mutate(CruiseSightNo = paste(Cruise,SightNo,sep="-"))%>%
  filter(CruiseSightNo %in% fkw_popData$CruiseSightNo)

#' Correct segdata based on pelagic sightings 
segDataESWg0$sightinfo <- segDataESWg0$sightinfo %>% left_join(.,fkw_popData[,c("CruiseSightNo","ProbPel","EffGrpSizeFull","EffGrpSizePelagic")], by="CruiseSightNo") 
segDataESWg0$sightinfo$nSI_033 <- 1
segDataESWg0$sightinfo$ANI_033 <- segDataESWg0$sightinfo$EffGrpSizeFull

sightinfo <- segDataESWg0$sightinfo %>% select(segnum, ProbPel, nSI_033, ANI_033)

segDataESWg0$segdata <- select(segDataESWg0$segdata, -nSI_033, -ANI_033)
segDataESWg0$segdata <- left_join(segDataESWg0$segdata, sightinfo)
segDataESWg0$segdata <- segDataESWg0$segdata %>% mutate(
  nSI_033 = ifelse(is.na(nSI_033), 0, nSI_033),
  ANI_033 = ifelse(is.na(ANI_033), 0, ANI_033),
  ProbPel = ifelse(is.na(ProbPel), 1, ProbPel)
)

write.csv(segDataESWg0$sightinfo, file.path(localwd, "output", "sightinfo_033.csv" ), row.names = FALSE)

#' -----------------------------------------------------------------------------
#'  Attach environmental data
#' -----------------------------------------------------------------------------

segEnvData <- readRDS(file.path(here(), "_R_code","task_3_download_env_data","output","segdata_env.rds")) %>% 
  select(segnum, Latitude, Longitude, UTC, sst:ssh_sd)
segdata <- left_join(segDataESWg0$segdata, segEnvData, by="segnum")

saveRDS(segdata, file.path(localwd, "output", "fkw_dsm_data_1997to2024.rds") )
write.csv(segdata, file.path(localwd, "output", "fkw_dsm_data_1997to2024.csv"), row.names = FALSE)

