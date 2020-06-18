source('estimCI.R')
source('fit_genGrowth.R')

set.seed(1234)


# --- Create simulated data:

n.data = 20
dat = data.frame(t = 1:n.data,
                 inc = rpois(n = n.data, 
                             lambda = exp(0.2 * c(1:n.data)) ))
plot(dat, typ='b', 
     main = 'Simulated data',
     xlab='time', ylab='incidence', las=1) 


# --- Fit a general growth model

# -- Fitting without confidence intervals:

prm.init = c(r=0.6, p=0.8)  # initial values for the optimization

prm.fit = fit.genGrowth.inc(dat,
                            prm.init,
                            fit.type = "LSconstraint",
                            relative.err = FALSE)

plot.data.fit(dat,prm.fit)
print(prm.fit)

# -- Fitting with confidence intervals:

est.prm = estimate.CI(dat      = dat, 
                      CIwidth  = 0.99, 
                      n.MC     = 500, 
                      prm.init = prm.init)

print(est.prm)








