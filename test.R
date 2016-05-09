source("ggm_lib.R")

inc <- c(1,3,2,7,12,19,25,22,36,77,55,66,45,38,12) #; plot(inc,log='y',typ='o')
dat <- data.frame(t=1:length(inc), inc=inc)

prm.init <- c(r=0.02 ,p=0.9,K=1E3,a=0.2)

fit <- fit.genGrowth4.inc(dat = dat, 
						  prm.init = prm.init,
						  fit.type = 'LSconstraint',
						  relative.err = F)


tvec <- 1:100
sim <- genGrowth4.inc(inc.init = 1, 
					  r = fit[['r']],p = fit[['p']], K = fit[['K']],a = fit[['a']], 
					  tvec = tvec)

par(mfrow=c(1,2))

myplot<- function(dolog){
	plot(tvec,sim,typ='l', log=dolog, lwd=3,ylim=range(sim,inc))
lines(inc,typ='o',col='blue')
}
myplot('')
myplot('y')
print(fit)



if(FALSE){
	fit <- estimate.CI(dat = data.frame(t=1:length(inc), inc=inc), 
					   CIwidth = 0.95, 
					   n.MC = 200, 
					   prm.init = c(r=0.001, p=0.1), 
					   relative.err = FALSE) 
	
	print(fit)
	
	tvec <- 1:(length(inc)+10)
	sim.md <- genGrowth.inc(c0=inc[1],
							r = fit[['r.md']],
							p = fit[['p.md']],
							tvec = tvec)
	
	sim.lo <- genGrowth.inc(c0=inc[1],
							r = fit[['r.ci']][1],
							p = fit[['p.ci']][1],
							tvec = tvec)
	
	sim.hi <- genGrowth.inc(c0=inc[1],
							r = fit[['r.ci']][2],
							p = fit[['p.ci']][2],
							tvec = tvec)
	
	plot(tvec, sim.md, typ='l',col='red',lwd=3, log='y',
		 ylim = range(sim.lo,sim.hi,inc))
	lines(tvec, sim.lo,col='red')
	lines(tvec, sim.hi,col='red')
	points(inc)
}