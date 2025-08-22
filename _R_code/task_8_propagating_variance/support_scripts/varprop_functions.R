#' -----------------------------------------------------------------------------
#' Calculate effort offset components for various parameters
#' -----------------------------------------------------------------------------

make_g0_esw <- function(par=NULL, newdata, ds_model, Lp_g0){
  ndet <- length(ds_model$ddf$par)
  ng0 <- ncol(Lp_g0)
  if(length(par) != ndet + ng0) stop("'par' is not correct size in 'effort_offset()' function.")
  ds_model$ddf$par <- tail(par,ndet)
  esw <- predict(ds_model, newdata, compute = TRUE, esw=TRUE)[["fitted"]]
  g0 <- as.vector(exp(Lp_g0 %*% head(par,ng0))); g0 <- g0/g0[1]
  BFlow <- floor(newdata$beaufort + 1)
  BFhigh <- ceiling(newdata$beaufort + 1)
  w <- BFhigh - (newdata$beaufort + 1)
  g0_bft <- w*g0[BFlow] + (1-w)*g0[BFhigh]
  return(list(g0=g0_bft, esw=esw))
}



#' -----------------------------------------------------------------------------
#' Add detection parameter uncertainty to fitted GAM model
#' -----------------------------------------------------------------------------

add_varprop <- function(gam_model, beaufort_data, varprop_list, trace=FALSE){
  require(mgcv)
  require(dsm)
  require(numDeriv)
  ds_model <- varprop_list$detfun
  # Some checks
  if(gam_model$family$link != "log") stop("log link must be used!") 
  if(ds_model$ddf$ds$aux$ddfobj$scale$formula == "~1") stop("varprop doesn't work when there are no covariates in the detection function")
  if(gam_model$method != "REML") stop("REML must be used for smoothing parameter selection")
  
  linkfn <- gam_model$family$linkfun
  linkinvfn <- gam_model$family$linkinv
  
  # extract the call
  this_call <- as.list(gam_model$call)
  # remove the function
  this_call[1] <- NULL
  
  # log offset function
  par <- c(varprop_list$g0_est$coefficients.t[varprop_list$bft_idx], ds_model$ddf$par)
  
  log_effort_offset <- function(par, newdata, ds_model, Lp_g0){
    ndet <- length(ds_model$ddf$par)
    ng0 <- ncol(Lp_g0)
    if(length(par) != ndet + ng0) stop("'par' is not correct size in 'effort_offset()' function.")
    ds_model$ddf$par <- tail(par,ndet)
    esw <- predict(ds_model, newdata, compute = TRUE, esw=TRUE)[["fitted"]]
    g0 <- as.vector(exp(Lp_g0 %*% head(par,ng0))); g0 <- g0/g0[1]
    BFlow <- floor(newdata$beaufort + 1)
    BFhigh <- ceiling(newdata$beaufort + 1)
    w <- BFhigh - (newdata$beaufort + 1)
    g0_bft <- w*g0[BFlow] + (1-w)*g0[BFhigh]
    return(log(g0_bft* esw))
  }
  
  newdata <- mgcv::uniquecombs(beaufort_data)
  Jmat <- numDeriv::jacobian(log_effort_offset, par, newdata=newdata, ds_model=ds_model, Lp_g0= varprop_list$Lp_g0)
  Jmat <- Jmat[attr(newdata, "index"),]
  
  dat <- gam_model$model %>% rename(effort = `(offset)`) %>% 
    mutate(effort = exp(effort))
  if("(weights)" %in% colnames(dat)) colnames(dat)[colnames(dat) =="(weights)"] <- as.character(this_call$weights)
    
           
  dat[["Jmat"]] <- Jmat
  
  #' Refit model with Jacobian
  this_call$formula <- update.formula(this_call$formula,
                                      paste0(paste(as.character(this_call$formula)[c(2,1,3)],
                                                   collapse=" "), " + Jmat"))
  this_call$data <- dat
  this_call$paraPen <- c(this_call$paraPen, list(Jmat=list(S=solve(varprop_list$V_par))))
  this_call$fixed.priors <- "Jmat"
  this_call$scale.trace <- trace
  this_call$scale.estimated <- gam_model$scale.estimated
  gam.fixed.priors <- dsm:::gam.fixed.priors
  refit <- do.call("gam.fixed.priors", this_call)
  
}

#' -----------------------------------------------------------------------------
#' Add detection parameter uncertainty to fitted GAM model using simulation
#' -----------------------------------------------------------------------------

#' add_varprop_sim <- function(gam_model, beaufort_data, varprop_list, trace=FALSE){
#'   require(mgcv)
#'   require(dsm)
#'   ds_model <- varprop_list$detfun
#'   #' Some checks
#'   if(gam_model$family$link != "log") stop("log link must be used!") 
#'   if(ds_model$ddf$ds$aux$ddfobj$scale$formula == "~1") stop("varprop doesn't work when there are no covariates in the detection function")
#'   if(gam_model$method != "REML") stop("REML must be used for smoothing parameter selection")
#'   
#'   #' extract data
#'   dat <- gam_model$model %>% rename(effort = `(offset)`) %>% 
#'     mutate(effort = exp(effort))
#'   #' extract the call
#'   this_call <- as.list(gam_model$call); this_call[1] <- NULL
#'   
#'   #' Fitted parameters
#'   par <- c(varprop_list$g0_est$coefficients.t[varprop_list$bft_idx], ds_model$ddf$par)
#'   
#'   #' Loop over g0, esw parameters
#'   refit <- do.call("gam.fixed.priors", this_call)
#'   
#' }
