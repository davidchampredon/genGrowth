read.inc <- function(filename, tstart,tend, doplot=TRUE) {
	# Read incidence data.
	# 1st column must be dates or times
	# 2nd column must be incidence (not cumulative)
	
	dat0 <- read.csv(filename, header = FALSE)
	names(dat0)<-c("date","inc")
	dat0$t <- 1:nrow(dat0)
	
	dat <- subset(dat0, t>tstart & t<tend)
	dat$t <- 1:nrow(dat)
	
	if(doplot){
		plot(dat0$t, dat0$inc, typ="s",
			 main = "Incidence Data", 
			 xlab="time", ylab="incidence",
			 col="grey")
		rng <- tstart:tend
		lines(dat0$t[rng], dat0$inc[rng],typ="s",lwd=6)
		grid()
	}
	return(dat)
}