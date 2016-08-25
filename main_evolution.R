
source('estimCI.R')
source('fit_evolution.R')
source('read-data.R')
source("../Datsid/read_db.R")
source("../Datsid/utils.R")

save.to.file <- T

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

# ==== PLOTS ==== 

if(save.to.file) png("fit_evol.png",width = 2000,height = 900)

par(mfrow=c(1,2), 
	mar = c(5, 6, 4, 2) + 0.1,
	cex.axis = 2, 
	cex.lab = 4,
	cex.main = 4,
	las = 1)

# ==== data  =====
plot(dat$t, dat$inc, log='y', typ='s', 
	 xlab = 'days',
	 ylab = '', #incidence',
	 main = "Incidence",# epichoice
	 lwd = 6
	 )
grid(equilogs = F)
points(dat.sub$t,dat.sub$inc, pch=15,cex=2.5)

mxinc = 1.15*max(dat.sub$inc)
shiftdown = 0.70
arrows(x0=dat.sub$t[1], x1=dat.sub$t[tmin],
	   y0 = shiftdown*mxinc, y1=shiftdown* mxinc, length = 0.1, angle = 60,code = 3)
arrows(x0=dat.sub$t[1], x1=dat.sub$t[length(dat.sub$t)],
		 y0 = mxinc, y1=mxinc, length = 0.1, angle = 60,code = 3)

text(x=0.5*(dat.sub$t[1]+dat.sub$t[tmin]), y=shiftdown*mxinc, 
	 labels = "first fit", 
	 cex = 1.7,pos = 1)

text(x=0.5*(dat.sub$t[1]+dat.sub$t[length(dat.sub$t)]), y=mxinc, 
	 labels = "last fit", 
	 cex = 1.7,pos = 1)

# ==== fit ====

lwd.ci <- 0.7

plot(x=r, y=p,
	 typ='o',
	 cex=1.7, pch=1, lty =2, 
	 main = "Parameters Fit", #epichoice,
	 ylim = c(0,1),
	 xlim = range(0,r.hi))
arrows(x0=r, x1=r, y0=p.lo, y1=p.hi, col='darkgrey',code = 3,angle = 90, length = 0.05, lwd = lwd.ci)
arrows(x0=r.lo, x1=r.hi, y0=p, y1=p, col='darkgrey',code = 3,angle = 90, length = 0.05, lwd = lwd.ci)
text(x=r[1], y=p[1], labels = "first fit", cex = 2,pos = 4)
text(x=r[n], y=p[n], labels = "last fit", cex = 2,pos = 4)

# text(x=r, y=p, labels = c(1:length(r)), cex = 0.8,pos = 4)

points(x=r, y=p, cex=1.2, pch=16)

if(save.to.file) dev.off()