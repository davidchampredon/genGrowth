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




# - - - - - - 
# tvec <- 0:30
# c0 <-2
# r <- 0.4
# p <- 0.9
# 
# cuminc <- genGrowth.cuminc(c0,r,p,tvec)
# inc <- genGrowth.inc(c0,r,p,tvec)
# 
# plot(tvec,ci,typ="s",log="y")
# lines(tvec,inc,typ="s",col="blue")
