## Create strata sf objects 
#load packages
library(tidyverse)
library(lubridate)
library(here)
library(swfscDAS)
library(MASS)
library(sf)
library(picMaps)

datadir<- here("analysis_tasks/task_2_g0_ESW")

CenPac <-read.csv(file.path(datadir, "Strata Files","CenPac2_broader.csv"))%>%
  st_as_sf(coords = c("Lon","Lat"))%>%
  st_bbox()%>%st_as_sfc()%>%
  st_set_crs(4326)


ETPCore <- read.csv(file.path(datadir, "Strata Files","ETPCore.csv"))%>%
  st_as_sf(coords = c("Lon","Lat"))%>%
  st_bbox()%>%st_as_sfc()%>%
  st_set_crs(4326)

ETPOuter <- read.csv(file.path(datadir, "Strata Files","ETPOuter.csv"))%>%
  st_as_sf(coords = c("Lon","Lat"))%>%
  st_bbox()%>%st_as_sfc()%>%
  st_set_crs(4326)


ggplot()+
  theme_classic()+
  geom_sf(data = CenPac%>%st_shift_longitude(), fill = NA, color = "green")+
  geom_sf(data = ETPCore%>%st_shift_longitude(), fill = NA, color = "blue")+
  geom_sf(data = ETPOuter%>%st_shift_longitude(), fill = NA, color = "red")+
  geom_sf(data = classify%>%st_shift_longitude(), size = 1, alpha = 0.5,color= "grey40")





strata_files<-list("CenPac" = CenPac%>%st_shift_longitude(), 
                   "ETPCore" = ETPCore%>%st_shift_longitude(), 
                   "ETPOuter" = ETPOuter%>%st_shift_longitude())
saveRDS(strata_files, (file.path(datadir, "Strata Files","strata_sf_files.rds")))

tmp<-st_union(strata_files[[1]],strata_files[[2]])
cookie<-st_union(tmp,strata_files[[3]])%>%
  st_shift_longitude()

saveRDS(cookie, (file.path(datadir, "Strata Files", "region_filter_poly_esw_g0.rds")))

cookie_coords<-st_coordinates(cookie)[,1:2]
colnames(cookie_coords)<-c("Lon","Lat")
cookie_coords[,"Lon"]<-cookie_coords[,"Lon"]-360

write.csv(cookie_coords, file = (file.path(datadir, "Strata Files", "region_filter_esw_g0.csv")), row.names = F)
