source("fit_genGrowth.R")

estimate.CI <- function(dat, CIwidth, n.MC, 
						prm.init, 
						relative.err = FALSE) {

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