source("read-data.R")
source("fit_genGrowth.R")
source("fit_evolution.R")
source("estimCI.R")

source("../Datsid/read_db.R")
source("../Datsid/utils.R")

set.seed(12)

main.analysis <- function(filename,tstart,tend, CIwidth, n.MC=1000) {
	# Read data:
	dat <- read.inc(filename, tstart, tend, 
					doplot = TRUE, header=TRUE)
	
	# Fitting (no confidence intervals)
	prm.init= c(r=0.6, p=0.8)
	prm.fit <- fit.genGrowth.inc(dat,prm.init,fit.type = "LSconstraint")
	plot.data.fit(dat,prm.fit)
	
	# Estimates (with confidence intervals)
	est.prm <- estimate.CI(dat = dat, 
						   CIwidth = CIwidth, 
						   n.MC = 10, 
						   prm.init = prm.init)
	print(filename)
	print(est.prm)
	
	fit.evolution.CI(filename = filename, 
					 tstart = tstart, 
					 tend = tend,
					 prm.init = prm.init, 
					 CIwidth = CIwidth, 
					 n.MC = n.MC, 
					 tstep=1, 
					 min.length = 5)
}


main.analysis.df <- function(dat,tstart,tend, CIwidth, n.MC=1000) {
	
	# Fitting (no confidence intervals)
	prm.init= c(r=0.6, p=0.8)
	prm.fit <- fit.genGrowth.inc(dat,prm.init,fit.type = "LSconstraint")
	plot.data.fit(dat,prm.fit)
	
	# Estimates (with confidence intervals)
	est.prm <- estimate.CI(dat = dat, 
						   CIwidth = CIwidth, 
						   n.MC = 10, 
						   prm.init = prm.init)
	return(est.prm)
	
	# fit.evolution.CI(filename = filename, 
	# 				 tstart = tstart, 
	# 				 tend = tend,
	# 				 prm.init = prm.init, 
	# 				 CIwidth = CIwidth, 
	# 				 n.MC = n.MC, 
	# 				 tstep=1, 
	# 				 min.length = 5)
}


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

main.analysis.db <- function(dbname, epichoice,
							 tstart, tend, 
							 CIwidth, n.MC=1000) {
	# Read data:
	dat <- get.db.data(dbname,epichoice)
	dat <- subset(dat, t<=tend & tstart<=t)
	
	# Fitting (no confidence intervals)
	prm.init <- c(r=0.6, p=0.8)
	prm.fit <- fit.genGrowth.inc(dat, prm.init,
								 fit.type = "LSconstraint")
	plot.data.fit(dat,prm.fit,epichoice)
	
	# Estimates (with confidence intervals)
	est.prm <- estimate.CI(dat = dat, 
						   CIwidth = CIwidth, 
						   n.MC = 10, 
						   prm.init = prm.init)
	
	if(FALSE){
		fit.evolution.CI(filename = filename,
						 tstart = tstart,
						 tend = tend,
						 prm.init = prm.init,
						 CIwidth = CIwidth,
						 n.MC = n.MC,
						 tstep = 1,
						 min.length = 5)
	}
	return(est.prm)
}

dbname <- "../Datsid/a.db"
epilist <- get.list.existing(dbname)
epilist

epis <- read.csv('episelect.csv', header = TRUE)
epikey <- as.character(epis$key)
tstart <- epis$start
tend <- epis$end

res <- list()

pdf(paste0("plot_ggm.pdf"), width = 18, height = 10)
for(i in 1:length(epikey)){
	res[[i]] <- main.analysis.db(dbname = dbname,
								 epichoice = epikey[i],
								 tstart = tstart[i],
								 tend = tend[i], 
								 CIwidth = 0.95, 
								 n.MC = 1000)
}
dev.off()

plot.p.r.fit <- function(res,epikey){

	for(i in 1:length(res)){
		if(i==1) {
			plot(x=res[[i]]$r.md, y=res[[i]]$p.md,
				 xlab = "r",
				 ylab = "p",
				 las = 1,
				 xlim = c(0.001,2), 
				 ylim = c(0.001,1),
				 log='')
			grid()
		}
		if(i>1) {
			points(x=res[[i]]$r.md, y=res[[i]]$p.md)
		}
		text(x=res[[i]]$r.md, y=res[[i]]$p.md,
			 labels = epikey[i],
			 pos = 3, cex = 0.5,offset = 1)
		segments(x0 = res[[i]]$r.ci[1], x1=res[[i]]$r.ci[2],
				 y0 = res[[i]]$p.md, y1=res[[i]]$p.md)
		segments(x0 = res[[i]]$r.md, x1=res[[i]]$r.md,
				 y0 = res[[i]]$p.ci[1], y1=res[[i]]$p.ci[2])
	}	
}

plot.p.r.fit(res,epikey)

epichoice <- "influenza.FRANCE."
z <- get.db.data(dbname, epichoice)

par(mfrow=c(1,1))
slice <- slice.recurrent.timeseries(dat = z,w = 15)

res.slice <- list()
pdf(paste0("plot_slice_",epichoice,".pdf"), width = 18, height = 10)
for(i in 1:length(slice)){
	res.slice[[i]] <- main.analysis.df(dat = slice[[i]],
									   tstart = 1 ,
									   tend = length(res.slice[[i]]$inc), 
									   CIwidth = 0.95, 
									   n.MC=1000)
}
dev.off()

plot.p.r.fit(res.slice,epikey = NULL)

# par(mfrow=c(5,6))
# for(i in 1:length(aa)) plot(aa[[i]]$t,aa[[i]]$inc, typ='o',pch=1,cex=0.6,log='y')


