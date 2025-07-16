
# Uses cruise data from 1986-2023 on false killer whales 
#  filtered for ETP, CNP (a broader CNP area than in the DSM) regions

#load packages
library(tidyverse)
library(Distance)
library(lubridate)
library(here)
library(swfscDAS)
library(MASS)
library(mvnfast)
library(sf)
library(picMaps)

localwd <- here("_R_code","task_2_est_g0_esw")
sp.code <- "033"

#read in input files
das_file <- file.path(localwd, "AllSurveys_g0_ESW_1986-2023.das")
y.proc <- swfscDAS::das_process(das_file)
strata_file_dir <- file.path(localwd, "Strata Files", "region_filter_esw_g0.csv")

y.eff <- swfscDAS::das_effort(
  y.proc, method = "equallength",
  seg.km = 10,
  dist.method = "greatcircle",
  num.cores = 1,
  strata_files = strata_file_dir)


y.eff.sight <- swfscDAS::das_effort_sight(y.eff, sp.codes = sp.code, sp.events = "S")

# Region filtering

y.eff.sight.strata <- swfscDAS::das_intersects_strata(y.eff.sight, list(InPoly = strata_file_dir))

seg<-y.eff.sight.strata$segdata%>%
  st_as_sf(coords = c("mlon","mlat"))%>%
  st_set_crs(4326)%>%
  st_shift_longitude()
region<-readRDS(file.path(localwd, "Strata Files", "region_filter_poly_esw_g0.rds"))%>%
  st_make_valid()%>%
  st_shift_longitude()


ggplot()+
  theme_classic()+
  geom_sf(data = seg, color = "blue")+
  geom_sf(data = region, color = "black", fill = NA)+
  geom_sf(data = seg[seg$InPoly==1,], color = "red")+
  geom_sf(data = picMaps::hawaii_coast(), fill = "grey")+
  geom_sf(data = picMaps::hawaii_eez(), fill = NA,color = "grey")


input_dat<-list()
input_dat$esw.dat <- y.eff.sight.strata$sightinfo %>% filter(SpCode == sp.code,
                                                         InPoly == 1,
                                                         OnEffort == TRUE)
input_dat$g0.dat <- y.eff.sight.strata$segdata %>% filter(InPoly == 1, EffType == "S")

saveRDS(input_dat, file.path(localwd, "output","esw_g0_input.rds"))
