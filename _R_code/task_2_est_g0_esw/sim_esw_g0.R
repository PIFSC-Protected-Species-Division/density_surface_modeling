sim_esw_g0 <- function(pars.esw, v_theta, detfun, g0_dat){
  i <- 1; good_ret <- 0
  while(good_ret==0 && i<6){
    detfun$ddf$par <- mvnfast::rmvn(1, pars.esw, v_theta)
    esw_star <- predict(detfun, newdata = g0_dat, compute=TRUE, esw = TRUE)$fitted 
    esa_star <- esw_star*2*g0_dat$dist
    g0_star <- try(scam::scam(formula = resp ~ s(I(pmax(beaufort-1,0)), bs='mpd', k=4)
                          + s(mlat, mlon, bs='tp') + s(year, k=4), 
                          offset = log(esa_star+0.0001),
                          family=poisson,
                          data=g0_dat,
                          not.exp=TRUE,
                          gamma=1.4), silent=TRUE)
    if(!inherits(g0_star, "try-error")){
      good_ret <- 1
      newdata_g0 <- data.frame(beaufort = 0:7,
                               mlat =mean(g0_dat$mlat, na.rm=TRUE),
                               mlon = mean(g0_dat$mlon, na.rm=TRUE),
                               year = round(mean(g0_dat$year, na.rm=TRUE)))
      
      ps <- scam::predict.scam(g0_star, newdata = newdata_g0, type="response")
      Rg0 <- ps/ps[1]
    } else{
      i <- i+1
    }
  }
  return(list(g0_star=g0_star, theta_star=detfun$ddf$par, Rg0=Rg0))
}
