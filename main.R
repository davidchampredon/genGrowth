source("read-data.R")
source("fit_genGrowth.R")

filename <- "./data/2009_Flu_Mexico.csv"
tstart <- 30
tend <- 44

dat <- read.inc(filename, tstart,tend)

# Fitting

prm.init= c(r=0.6, p=0.8)

prm.fit <- fit.genGrowth.inc(dat,prm.init)

plot.data.fit(dat,prm.fit)
