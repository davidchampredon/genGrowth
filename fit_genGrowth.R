source("genGrowth.R")

#' Least-square error function.
error.LS <- function(prm,dat,relative) {
  tvec <- dat$t
  r <- prm['r']
  p <- prm['p']
  target.inc <- dat$inc
  c0 <- target.inc[1]
  ggm <- genGrowth.inc(c0,r,p,tvec)
  res <- sqrt(sum((ggm - target.inc)^2))
  if(relative) res <- sqrt(sum((ggm/target.inc-1)^2 ))
  return(res)
}

#' Fit a generalized growth model.
#' @param dat Dataframe storing the incidence data. Must have column \code{t} for time and \code{inc} for incidence.
#' @param prm.init Named numerical vector. Initial values for the optimization algorithm. Must have elements named \code{r} and \code{p}.
#' @param fit.type String. Type of fit. \code{LS} for least-square, \code{LSconstraint} for least square with constraint on the parameters.
#' @param relative.err Boolean. Is the least square based on the relative error (between the data and the proposed incidence)?
#' @return A list with the value of the fitted parameters (r,p).
#' 
fit.genGrowth.inc <- function(dat, 
                              prm.init, 
                              fit.type="LS",
                              relative.err = FALSE) {
  # Fit parameters 'r' and 'p' 
  # for generalized growth model.
  
  if(fit.type=="LS"){
    fit <- optim(par = prm.init,
                 fn = error.LS, 
                 dat = dat,
                 relative = relative.err)
    res <- fit$par
  }
  if(fit.type=="LSconstraint"){
    fit <- nlminb(start = prm.init,
                  objective = error.LS,
                  lower = c(0,0), 
                  upper = c(99,0.999999), 
                  dat = dat,
                  relative = relative.err)
    res <- fit$par
  }
  
  return(res)
}


plot.data.fit <- function(dat,prm.fit,epichoice=NULL) {
  
  ggm.fit <- genGrowth.inc(c0 = dat$inc[1], 
                           r = prm.fit['r'], 
                           p = prm.fit['p'], 
                           tvec = dat$t)
  
  old.par <- par()
  par(mfrow=c(1,2))
  
  plot(dat$t, dat$inc, pch=16,
       main = paste(epichoice,"Fitted GGM\n p =",round(prm.fit['p'],3),"  r =",round(prm.fit['r'],3)),
       xlab = "time", ylab = "incidence", las = 1)
  lines(dat$t, dat$inc, typ="s")
  lines(dat$t, ggm.fit,col="red")
  grid()
  
  plot(dat$t, dat$inc, pch=16,
       main = paste("Fitted GGM (log scale)\n p =",round(prm.fit['p'],3),"  r =",round(prm.fit['r'],3)),
       xlab = "time", ylab = "incidence", log="y", las=1)
  lines(dat$t, dat$inc, typ="s")
  lines(dat$t, ggm.fit,col="red")
  
  par(old.par)
}
