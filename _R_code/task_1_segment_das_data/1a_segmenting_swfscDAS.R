#Segmenting using swfscDAS

library(here)
library(swfscDAS)
library(tidyverse)

localwd <- here("_R_code","task_1_segment_das_data")
das_file <- file.path(localwd, "CenPac_Toolbox_1997-2023.das")
strata_file <- file.path(localwd, "CenPac2.csv")
# randPics <- readRDS(file.path(localwd, "output","seg_sight_out.rds"))[["randpicks"]] #if repeating same das file 
sp.code <- "033" # false killer whale sp.code


y.proc <- swfscDAS::das_process(das_file)

y.eff <- swfscDAS::das_effort(
  y.proc, method = "equallength", 
  seg.km = 10,
  dist.method = "greatcircle",
  num.cores = 1, 
  strata_files = strata_file
  # , randpicks.load = randPics #
)

y.eff.sight <- swfscDAS::das_effort_sight(y.eff, sp.codes = sp.code, sp.events = "S") # the sp.code here returns the number of included sightings and animals per segment. The sightings object is unaffected (has all sp.codes)

y.eff.sight.strata <- swfscDAS::das_intersects_strata(y.eff.sight, list(InPoly = strata_file))

#Filter for study area
y.eff.sight.strata$segdata <- y.eff.sight.strata$segdata %>% filter(InPoly == 1)
y.eff.sight.strata$sightinfo <- y.eff.sight.strata$sightinfo %>% filter(InPoly == 1)

#Filter for segment lengths > 0 that do not have an associated sighting
zero_dist_segs <- y.eff.sight.strata$segdata$segnum[y.eff.sight.strata$segdata$dist<0.1]
for(k in 1:length(zero_dist_segs)){
  segnum_k <- zero_dist_segs[k]
  y.eff.sight.strata$segdata$dist[y.eff.sight.strata$segdata$segnum==segnum_k] <- ifelse(y.eff.sight.strata$segdata$nSI_033[y.eff.sight.strata$segdata$segnum==segnum_k] > 0, 0.15,0)
}

y.eff.sight.strata$segdata <- y.eff.sight.strata$segdata %>% filter(dist > 0.1)

# fix NAs in aveBft

if(any(is.na(y.eff.sight.strata$segdata$avgBft))){
  fix_bft<-data.frame(NA_bft = which(is.na(y.eff.sight.strata$segdata$avgBft)), fixed = c(5,6,6,5,5,4,4)) #bft values provided from das file. Error in data entry resulted in NA values. 
  y.eff.sight.strata$segdata[fix_bft$NA_bft,"avgBft"] <- fix_bft$fixed
}
#filter sight info for S,G events

y.eff.sight.strata$sightinfo <- y.eff.sight.strata$sightinfo%>%filter(Event %in% c("S","G"))

if (!dir.exists(file.path(localwd,"output"))) {
  dir.create(file.path(localwd,"output"))
}

saveRDS(y.eff.sight.strata, file.path(localwd,"output", "seg_sight_out.rds"))


