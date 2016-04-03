source("read-data.R")
source("fit_genGrowth.R")
source("fit_evolution.R")
source("estimCI.R")

source("../Datsid/read_db.R")

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


main.analysis.db <- function(dbname, disease, country,
							 tstart,tend, 
							 CIwidth, n.MC=1000) {
	# Read data:
	# dat.old <- read.inc(filename, tstart, tend, 
	# 				doplot = TRUE, header=TRUE)
	
	dat.db <- get.epi.ts(dbname,
						 synthetic = 0,
						 country = country,
						 disease = disease) 
	
	dat <- data.frame(date = dat.db$reportdate, 
					  t = )
	
	# Fitting (no confidence intervals)
	prm.init <- c(r=0.6, p=0.8)
	prm.fit <- fit.genGrowth.inc(dat, prm.init,
								 fit.type = "LSconstraint")
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



dbname <- "../Datsid/bcktest.db"
epilist <- get.list.existing(dbname)
epilist

disease <- c("ebola","MERS-CoV")
country <- c("RDC","SOUTHKOREA")

tstart <- c(1,30,1)
tend <- tstart + c(30,15,22)

for(i in 1:length(fname)){
	pdf(paste0("plot_",i,".pdf"),width = 18, height = 10)
	main.analysis(filename = fname[i],
				  tstart = 	tstart[i],tend = tend[i], 
				  CIwidth=0.95, n.MC=1000)
	
	main.analysis.db(filename = fname[i],
				  tstart = 	tstart[i],tend = tend[i], 
				  CIwidth=0.95, n.MC=1000)
	
	dev.off()
}
