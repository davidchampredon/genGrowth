library(ggplot2)
source("read-data.R")
source("fit_genGrowth.R")
source("fit_evolution.R")
source("estimCI.R")
source("../Datsid/read_db.R")
source("../Datsid/utils.R")

t1 <- as.numeric(Sys.time())
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
	return(c(epi=epichoice,est.prm))
}


plot.p.r.fit <- function(res,epilabel,title=""){
	
	if(is.null(epilabel)) epilabel <- 1:length(res)
	
	tmp <- vector()
	for(i in 1:length(res)) tmp[i] <- res[[i]]$r.ci[2]
	rmax <- max(tmp)
	
	for(i in 1:length(res)){
		if(i==1) {
			plot(x=res[[i]]$r.md, y=res[[i]]$p.md,
				 xlab = "r",
				 ylab = "p",
				 las = 1,
				 main = title,
				 xlim = c(0.001,rmax), 
				 ylim = c(0.001,1),
				 log='')
			grid()
		}
		if(i>1) {
			points(x=res[[i]]$r.md, y=res[[i]]$p.md)
		}
		text(x=res[[i]]$r.md, y=res[[i]]$p.md,
			 labels = epilabel[i],
			 pos = 3, cex = 0.5,offset = 1)
		segments(x0 = res[[i]]$r.ci[1], x1=res[[i]]$r.ci[2],
				 y0 = res[[i]]$p.md, y1=res[[i]]$p.md)
		segments(x0 = res[[i]]$r.md, x1=res[[i]]$r.md,
				 y0 = res[[i]]$p.ci[1], y1=res[[i]]$p.ci[2])
	}	
}


fit.ggm.recurrent.epi <- function(dbname, epichoice, w) {
	### Fit 'p,r' on recurrent epidemic time series.
	### Identify every growth phases and fit individualy
	### on this 'slice'
	
	pdf(paste0("plot_slice_",epichoice,"_details.pdf"), 
		width = 18, height = 10)
	
	dat <- get.db.data(dbname, epichoice)
	par(mfrow=c(1,2))
	slice <- slice.recurrent.timeseries(dat = dat,w)
	
	res.slice <- list()
	for(i in 1:length(slice)){
		res.slice[[i]] <- main.analysis.df(dat = slice[[i]],
										   tstart = 1 ,
										   tend = length(slice[[i]]$inc), 
										   CIwidth = 0.95, 
										   n.MC = 1000)
		res.slice[[i]] <- c(epi=epichoice,slice=i,res.slice[[i]])
	}
	dev.off()
	
	pdf(paste0("plot_slice_",epichoice,"_final.pdf"), 
		width = 18, height = 10)
	par(mfrow=c(1,1))
	plot.p.r.fit(res.slice, epilabel = NULL, title = epichoice)
	dev.off()
	return(res.slice)
}


list.to.df <- function(L){
	return(data.frame(r.md = L$r.md,
					  r.ci.lo = L$r.ci[1],
					  r.ci.hi = L$r.ci[2],
					  p.md = L$p.md,
					  p.ci.lo = L$p.ci[1],
					  p.ci.hi = L$p.ci[2],
					  epi = L$epi))
}

# -  -  - - - - -   - - -  - - - - - -  - - ---    --   -  - -- 
# -  -  - - - - -   - - -  - - - - - -  - - ---    --   -  - -- 

dbname <- "../Datsid/a.db"
# epilist <- get.list.existing(dbname)
# epilist

epis <- read.csv('episelect.csv', header = TRUE)
epis
epikey <- as.character(epis$key)
tstart <- epis$start
tend <- epis$end


# aa <- get.epi.ts(db.path = dbname,country = 'LIBERIA',
# 				 disease = 'ebola',synthetic = 0)
# a <- subset(aa,eventtype=='incidence')
# a$t <- date.to.duration(a$reportdate)
# par(mfrow=c(1,1))
# z <- subset(a,100<=t & t<=235)
# plot(z$t, z$count,typ='o')

res <- list()
pdf(paste0("plot_singleEpi_details.pdf"), width = 18, height = 10)
for(i in 1:length(epikey)){
	res[[i]] <- main.analysis.db(dbname = dbname,
								 epichoice = epikey[i],
								 tstart = tstart[i],
								 tend = tend[i], 
								 CIwidth = 0.95, 
								 n.MC = 1000)
}
dev.off()
pdf(paste0("plot_singleEpi_final.pdf"), width = 18, height = 10)
plot.p.r.fit(res,epikey)
dev.off()

cityUK <- paste0("measles.UK.",
				 c("Birmingham","Bristol","Liverpool","London",
				   "Manchester","Newcastle","Sheffield"))

res.recur.fr <- fit.ggm.recurrent.epi(dbname = dbname,epichoice = "influenza.FRANCE.",w = 15)
res.recur.uk <- sapply(X = cityUK, FUN = fit.ggm.recurrent.epi, dbname=dbname, w=50)

# Make a big dataframe with all results:
for(i in 1:length(res)){
	if(i==1) df <- list.to.df(res[[i]])
	if(i>1) df <- rbind(df, list.to.df(res[[i]]) )
}
for(i in 1:length(res.recur.fr)){
	df <- rbind(df, list.to.df(res.recur.fr[[i]]) )
}
for(i in 1:length(res.recur.uk)){
	for(j in 1:length(res.recur.uk[[i]]))
		df <- rbind(df, list.to.df(res.recur.uk[[i]][[j]]) )
}

z <- gregexpr(pattern = ".",text = df$epi,fixed=TRUE)
stopchar <- matrix(unlist(z),ncol = 2)[,1]

df$disease <- substr(as.character(df$epi),
					 start=1,
					 stop=5)

pdf("plot_ALL_final.pdf",width=18,height = 10)
alpha <- 0.8
g <- ggplot(df)
g <- g + geom_pointrange(aes(x=r.md,y=p.md,
							 ymin=p.ci.lo,ymax=p.ci.hi,
							 colour = disease, shape=disease), alpha=alpha)
g <- g + geom_segment(aes(x=r.ci.lo,xend=r.ci.hi,y=p.md,yend=p.md,
						  colour = disease), alpha=alpha)
g <- g + scale_x_log10()
plot(g)

g <- g + facet_wrap(~epi)
plot(g)
dev.off()

# - - - - - - 
t2 <- as.numeric(Sys.time())
message(paste0("\n\n --- Completed in ",round((t2-t1)/60,1)," min ---\n"))
