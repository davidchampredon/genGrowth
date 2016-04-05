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


if(FALSE){
tt <- 1:30
c0 <- 1
r <- 0.1
p <- 0.7
yy <- genGrowth.inc(c0,r,p,tt)
yy2 <- genGrowth.inc(c0,r*1.9,p*1.0,tt)
plot(tt,yy,typ='l',log='y')
lines(tt,yy2,col='red')
}