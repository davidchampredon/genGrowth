source("genGrowth.R")

error.LS <- function(prm,dat) {
	# Least-square error function.
	tvec <- dat$t
	r <- prm['r']
	p <- prm['p']
	target.inc <- dat$inc
	c0 <- target.inc[1]
	ggm <- genGrowth.inc(c0,r,p,tvec)
	return(sum((target.inc - ggm)^2 ))
}



fit.genGrowth.inc <- function(dat, prm.init, fit.type="LS") {
	# Fit parameters 'r' and 'p' 
	# for generalized growth model.
	
	if(fit.type=="LS"){
		fit <- optim(par = prm.init,
					 fn = error.LS, 
					 dat = dat)
		res <- fit$par
	}
	
	return(res)
}

plot.data.fit <- function(dat,prm.fit) {
	
	ggm.fit <- genGrowth.inc(c0 = dat$inc[1], 
							 r = prm.fit['r'], 
							 p = prm.fit['p'], 
							 tvec = dat$t)
	
	plot(dat$t, dat$inc, pch=16,
		 main = paste("Fitted GGM\n p =",round(prm.fit['p'],3),"  r =",round(prm.fit['r'],3)),
		 xlab = "time", ylab = "incidence")
	lines(dat$t, dat$inc, typ="s")
	lines(dat$t, ggm.fit,col="red")
	grid()
}