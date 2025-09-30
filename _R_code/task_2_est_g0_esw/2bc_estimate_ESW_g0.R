#' -----------------------------------------------------------------------------
#' Estimate effective strip width (ESW) and trackline probability (g0) 
#' -----------------------------------------------------------------------------

#' load packages
library(tidyverse)
library(Distance)
library(lubridate)
library(here)
library(mvnfast)
library(picMaps)
library(future)
library(doFuture)
library(coda)

localwd <- here("_R_code","task_2_est_g0_esw")
sp.code <- "033"


#' Recall saved data from step 1
input_dat <- readRDS(file.path(localwd, "output","esw_g0_input.rds"))


#' -----------------------------------------------------------------------------
#' Step 1: Estimate ESW parameters 
#' -----------------------------------------------------------------------------
df <- data.frame(distance = input_dat$esw.dat$PerpDistKm, beaufort = input_dat$esw.dat$Bft)
df <- df[complete.cases(df),]

detfun <- ds(df, formula = ~beaufort, truncation = 5.5, key="hn") #for fkw
# plot(detfun, showpoints = FALSE)
summary(detfun)

# Predict ESW for modeling segments
modeling_data <- readRDS(file.path(here(), "_R_code", "task_1_segment_das_data", "output","seg_sight_out.rds"))
modeling_data$segdata <- modeling_data$segdata %>% rename(beaufort = avgBft)
modeling_data$segdata[paste0("esw","_",sp.code)] <- predict(detfun, modeling_data$segdata, esw=TRUE)[["fitted"]]


#' -----------------------------------------------------------------------------
#' Step 2: Estimate g0 for modeling segments using broad-scale data
#' -----------------------------------------------------------------------------
g0_dat <- input_dat$g0.dat
g0_dat$resp <- as.logical(g0_dat[,paste0("nSI_",sp.code)])
g0_dat <- g0_dat%>%dplyr::select(resp, mlat, mlon, avgBft, year,dist) %>% 
  rename(beaufort = avgBft)
g0_dat <- g0_dat[complete.cases(g0_dat), ]

g0_dat$esa <- predict(detfun, newdata = g0_dat, esw=TRUE)[["fitted"]] * 2 * g0_dat$dist
g0_dat <- filter(g0_dat, dist>0)

#' Combine beaufort=0,1 as the same g(0) = 1
g0_est <- scam::scam(formula = resp ~ s(I(pmax(beaufort-1,0)), bs='mpd', k=4)
                     + s(mlat, mlon, bs='tp') + s(year, k=4), 
                   offset = log(esa+0.0001),
                   family=poisson,
                   data=g0_dat,
                   not.exp=TRUE,
                   gamma=1.4)

#' predict relative g(0) for each Beaufort level
newdata_g0 <- data.frame(beaufort = 0:7,
                         mlat = mean(modeling_data$segdata$mlat, na.rm=TRUE),
                         mlon = mean(modeling_data$segdata$mlon, na.rm=TRUE),
                         year = round(mean(modeling_data$segdata$year, na.rm=TRUE)))

ps <- scam::predict.scam(g0_est, newdata = newdata_g0, type="response")
Rg0 <- ps/ps[1]

#' Linear interpolation for split segments with changing Beaufort level
BFlow <- floor(modeling_data$segdata$beaufort + 1)
BFhigh <- ceiling(modeling_data$segdata$beaufort + 1)
w <- BFhigh - (modeling_data$segdata$beaufort + 1)
modeling_data$segdata[,paste0("g0_w_",sp.code)] <- w
modeling_data$segdata[,paste0("g0","_",sp.code)] <- w*Rg0[BFlow] + (1-w)*Rg0[BFhigh]

saveRDS(modeling_data, file.path(localwd, "output","seg_sight_out_g0_ESW.rds"))

#' -----------------------------------------------------------------------------
#' Step 3: Simulate the previous process to obtain a covariance matrix for the 
#' detection parameters
#' -----------------------------------------------------------------------------

# Souce parameter simulation
source(file.path(localwd, "sim_esw_g0.R"))

# Design matrix for g(0) estimates 
Lp <- scam::predict.scam(g0_est, newdata = newdata_g0, type="lpmatrix")
# Subset for columns affecting g(0):
bft_idx <- grep("beaufort", colnames(Lp))
Lp <- Lp[,bft_idx]

# Draw detection parameter sample:
detfun_g0 <- detfun
v_theta <- solve(detfun_g0$ddf$hessian)
pars.esw <- detfun_g0$ddf$par

npar <- 10
samps <- 500 # total samples = npar * samps

plan("multisession", workers=npar)
set.seed(8675309)

g0_par_sim <- foreach(i = 1:npar,  .options.future = list(seed = TRUE), .errorhandling = "pass") %dofuture% {
  theta_star <- matrix(NA, samps, length(detfun$ddf$par))
  v_alpha_star <- matrix(NA, samps, length(bft_idx)^2)
  alpha_star <- matrix(NA, samps, length(bft_idx))
  Rg0_star <- matrix(NA, samps, 8)
  # Constructing variance-covariance matrix 
  for(j in 1:samps){
    g0_star_out <- sim_esw_g0(pars.esw, v_theta, detfun, g0_dat)
    g0_star <- g0_star_out$g0_star
    theta_star[j,] <-  g0_star_out$theta_star
    v_alpha_star[j,] <- as.vector(g0_star$Vp.t[bft_idx, bft_idx])
    alpha_star[j,] <- g0_star$coefficients.t[bft_idx]
    Rg0_star[j,] <- g0_star_out$Rg0
    
  }
  list(theta_star=theta_star, alpha_star=alpha_star, v_alpha_star=v_alpha_star, Rg0_star=Rg0_star)
}
plan("sequential")

theta_star <- lapply(g0_par_sim, \(x) x$theta_star) %>% do.call(rbind,.)
alpha_star <- lapply(g0_par_sim, \(x) x$alpha_star) %>% do.call(rbind,.)
v_alpha_star <- lapply(g0_par_sim, \(x) x$v_alpha_star) %>% do.call(rbind,.)
Rg0_star <- lapply(g0_par_sim, \(x) x$Rg0_star) %>% do.call(rbind,.)

## Compose full variance/covariance from law of total variance
v_alpha <- colMeans(v_alpha_star) %>% matrix(.,sqrt(length(.))) +
  var(alpha_star)
c_alpha_theta <- cov(alpha_star, theta_star)
V <- rbind(
  cbind(v_alpha, c_alpha_theta),
  cbind(t(c_alpha_theta), var(theta_star))
)

varprop_list <- list(
  newdata_g0 = newdata_g0,
  Lp_g0 = Lp,
  V_par = V,
  par_sim = cbind(alpha_star, theta_star),
  Rg0_sim = Rg0_star,
  detfun = detfun, 
  g0_dat = g0_dat,
  g0_est = g0_est,
  bft_idx = bft_idx
)

saveRDS(varprop_list, file = file.path(localwd, "output","varprop_list.rds"))

#' -----------------------------------------------------------------------------
#' Compute SDs and CIs for ESW and g(0) estimates
#' -----------------------------------------------------------------------------

g0_sim <- matrix(nrow=nrow(alpha_star), ncol=8)
for(i in 1:nrow(g0_sim)){
  g0_sim[i,] <- exp(Lp%*%alpha_star[i,])/exp(Lp%*%alpha_star[i,])[1]
}
## g(0) standard error
round(apply(g0_sim, 2, sd), 3)
## g(0) CI
round(HPDinterval(mcmc(g0_sim)),2)

esw_sim <- matrix(nrow=nrow(alpha_star), ncol=8)
for(i in 1:nrow(g0_sim)){
  detfun_sim <- detfun
  detfun_sim$ddf$par <- theta_star[i,]
  esw_sim[i,] <- predict(detfun_sim, newdata = newdata_g0, compute=TRUE, esw = TRUE)$fitted 
}
## ESW standard error
round(apply(esw_sim, 2, sd), 2)
## ESW CI
round(HPDinterval(mcmc(esw_sim)),2)


