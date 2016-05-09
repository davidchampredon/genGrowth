library(deSolve)

genGrowth.cuminc <- function(c0,r,p,tvec) {
	# Generalized growth formula 
	# for cumulative incidence:
	m <- 1/(1-p)
	a <- c0^(1/m)
	return( (r/m*tvec + a)^m )
}

genGrowth.inc <- function(c0,r,p,tvec) {
	# Generalized growth formula 
	# for period incidence:
	m <- 1/(1-p)
	a <- c0^(1/m)
	return( r*(r/m*tvec + a)^(m-1) )
}




ggm4.ode <- function(t,x,parms){
	with(as.list(c(parms,x)),
		 {
		 	dC <- r*exp(p*log(C))*(1-(C/K)^a)
		 	res <- c(dC)
		 	list(res)
		 })
}

genGrowth4.inc <- function(inc.init, r,p,K,a, tvec){
	
	parms <- c(r=r, p=p, K=K, a=a)
	times <- tvec
	cuminc <- lsoda(y = c(C=inc.init),
					times, 
					func = ggm4.ode, 
					parms)	
	inc <- c(cuminc[1,2],diff(cuminc[,2]))
	return(inc)
}

error.LS <- function(prm,dat,relative) {
	# Least-square error function.
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

error.ggm4.LS <- function(prm,dat,relative) {
	# Least-square error function.
	tvec <- dat$t
	r <- prm[['r']]
	p <- prm[['p']]
	K <- prm[['K']]
	a <- prm[['a']]
	target.inc <- dat$inc
	c0 <- target.inc[1]
	ggm <- genGrowth4.inc(inc.init = c0, r=r, p=p,K = K, a=a, tvec = tvec)
	res <- sqrt(sum((ggm - target.inc)^2))
	if(relative) res <- sqrt(sum((ggm/target.inc-1)^2 ))
	return(res)
}


fit.genGrowth.inc <- function(dat, prm.init, fit.type="LS",relative.err) {
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

fit.genGrowth4.inc <- function(dat, prm.init, fit.type="LS",relative.err) {
	# Fit parameters 'r' and 'p' 
	# for generalized growth model.
	
	if(fit.type=="LS"){
		fit <- optim(par = prm.init,
					 fn  = error.ggm4.LS, 
					 dat = dat,
					 relative = relative.err)
		res <- fit$par
	}
	if(fit.type=="LSconstraint"){
		fit <- nlminb(start = prm.init,
					  objective = error.ggm4.LS,
					  lower = c(0,0,0,0), 
					  upper = c(9.9,0.99,1E7,9.9), 
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
		 xlab = "time", ylab = "incidence")
	lines(dat$t, dat$inc, typ="s")
	lines(dat$t, ggm.fit,col="red")
	grid()
	
	plot(dat$t, dat$inc, pch=16,
		 main = paste("Fitted GGM (log scale)\n p =",round(prm.fit['p'],3),"  r =",round(prm.fit['r'],3)),
		 xlab = "time", ylab = "incidence", log="y")
	lines(dat$t, dat$inc, typ="s")
	lines(dat$t, ggm.fit,col="red")
	grid()
	par(old.par)
}


estimate.CI <- function(dat, CIwidth, n.MC, prm.init, relative.err) {
	
	inc <- dat$inc
	fit <- matrix(nrow = n.MC, ncol = 2)
	for (i in 1:n.MC) {
		# Sample an epidemic centred on real data:
		inc.sample <- rpois(n=length(inc),lambda = inc)
		dat.sample <- data.frame(t=dat$t, inc=inc.sample)
		# fit this new synthetic epidemic:
		fit[i,] <- fit.genGrowth.inc(dat = dat.sample,
									 prm.init = prm.init,
									 fit.type = "LSconstraint",
									 relative.err = relative.err)
	}
	a <- (1-CIwidth)/2
	r.md <- quantile(fit[,1],probs = 0.50)
	p.md <- quantile(fit[,2],probs = 0.50)
	r.ci <- quantile(fit[,1],probs = c(a,1-a))
	p.ci <- quantile(fit[,2],probs = c(a,1-a))
	
	return(list(r.md=r.md, p.md=p.md, r.ci=r.ci, p.ci=p.ci))
}