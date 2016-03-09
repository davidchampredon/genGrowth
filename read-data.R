read.inc <- function(filename, tstart,tend, doplot=FALSE, ...) {
	# Read incidence data.
	# 1st column must be dates or times
	# 2nd column must be incidence (not cumulative)
	
	dat0 <- read.csv(filename,...)
	names(dat0)<-c("date","inc")
	dat0$t <- 1:nrow(dat0)
	
	dat <- subset(dat0, t>tstart & t<tend)
	dat$t <- 1:nrow(dat)
	
	if(doplot){
		old.par <- par()
		par(mfrow=c(1,2))
		rng <- tstart:tend
		plot(dat0$t, dat0$inc, typ="s",
			 main = paste("Incidence Data\n",filename), 
			 xlab="time", ylab="incidence",
			 col="grey")
		lines(dat0$t[rng], dat0$inc[rng],typ="s",lwd=6)
		grid()
		plot(dat0$t, dat0$inc, typ="s",
			 main = "Incidence Data\n(Log scale)", 
			 xlab="time", ylab="incidence",
			 col="grey",log="y")
		lines(dat0$t[rng], dat0$inc[rng],typ="s",lwd=6)
		
		par(old.par)
	}
	return(dat)
}