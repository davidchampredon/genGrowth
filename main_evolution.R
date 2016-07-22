
source('estimCI.R')
source('fit_evolution.R')
source('read-data.R')
source("../Datsid/read_db.R")
source("../Datsid/utils.R")


get.db.data <- function(dbname, epichoice){
	
	# Read data:
	dat.db <- get.epi.ts(dbname,
						 synthetic = 0,
						 country = NULL,
						 disease = NULL) 
	
	# make sure it's incidence:
	dat.db <- subset(dat.db, eventtype == 'incidence')
	dat.db$epikey <- paste(sep=".",
						   dat.db$disease_name,
						   dat.db$country,
						   dat.db$adminDiv1)
	dat.db <- subset(dat.db, epikey %in% trimws(epichoice))
	
	# reformat in simple data frame:	
	dd <- dat.db$reportdate
	
	dat <- data.frame(t = date.to.duration(dd),
					  inc = dat.db$count)
	return(dat)
}

dbname <- "../Datsid/a.db"
episelect <- read.csv('episelect.csv')

ee <- 5   # 1, 4, *5*

epichoice <- as.character(episelect$key[ee])
sub.start <- episelect$start[ee]
sub.end <- episelect$end[ee]


dat <- get.db.data(dbname, epichoice)
dat.sub <- subset(dat, t>=sub.start & t <= sub.end)

prm.init <- c(r=1, p=0.8)
CIwidth <- 0.95
n.MC <- 100

dat.fit <- dat.sub
dat.fit$t <- dat.fit$t - dat.fit$t[1] +1

res <- list()

tmin = 6

cnt = 1
for (i in seq(tmin,nrow(dat.fit), by = 1)) {
	tmp <- dat.fit[1:i, ]
	res[[cnt]] <- estimate.CI(dat = tmp, 
							  CIwidth, 
							  n.MC, 
							  prm.init, 
							  relative.err = FALSE) 
	cnt <- cnt + 1
}

n <- length(res)
r <- numeric(n)
p <- numeric(n)
r.lo <- numeric(n)
p.lo <- numeric(n)
r.hi <- numeric(n)
p.hi <- numeric(n)
for (i in 1:n) {
	r[i] <- res[[i]]$r.md
	p[i] <- res[[i]]$p.md
	r.lo[i] <- res[[i]]$r.ci[1]
	r.hi[i] <- res[[i]]$r.ci[2]
	p.lo[i] <- res[[i]]$p.ci[1]
	p.hi[i] <- res[[i]]$p.ci[2]
}

par(mfrow=c(1,2), cex.axis = 2, cex.lab = 2,las = 1)

# data 
plot(dat$t, dat$inc, log='y', typ='s', 
	 xlab = 'days',
	 ylab = 'incidence',
	 main = epichoice)
points(dat.sub$t,dat.sub$inc, pch=16)

# fit 
plot(x=r, y=p,
	 typ='o',
	 cex=1, pch=16, lty =2, 
	 main = epichoice,
	 ylim = c(0,1),
	 xlim = range(0,r.hi))
arrows(x0=r, x1=r, y0=p.lo, y1=p.hi, col='darkgrey',code = 3,angle = 90, length = 0.05, lwd = 0.5)
arrows(x0=r.lo, x1=r.hi, y0=p, y1=p, col='darkgrey',code = 3,angle = 90, length = 0.05, lwd = 0.5)
text(x=r[1], y=p[1], labels = "start", cex = 2,pos = 3)
text(x=r[n], y=p[n], labels = "end", cex = 2,pos = 3)
