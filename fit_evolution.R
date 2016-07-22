source("fit_genGrowth.R")



fit.evolution <- function(filename, tstart, tend, prm.init, 
						  relative.err,
						  tstep=1, min.length = 5) {
	# Plots the values of fitted parameters
	# as the number of data increases
	
	tend2 <- tstart+min.length
	i <- 1
	prm.fit <- list()
	while (tend2 <= tend) {
		dat <- read.inc(filename, tstart,tend2)
		prm.fit[[i]] <- fit.genGrowth.inc(dat,prm.init,fit.type = "LSconstraint")
		tend2 <- tend2 + tstep
		i <- i + 1
	}
	
	M <- matrix(unlist(prm.fit), ncol = 2, byrow = TRUE)
	par(mfrow=c(1,1))
	plot(x=M[,1], y=M[,2],
		 xlim = c(0,max(M[,1])),
		 ylim = c(0,1),
		 xlab="r",ylab="p",typ="o",
		 col=rgb(0,0,0,0.5), pch=16)
	text(M[1,1],M[1,2],labels = "start",pos=3,cex=0.5)
	points(M[nrow(M),1],M[nrow(M),2],pch=0,col="red",cex=2,lwd=5)
	grid()
}


fit.evolution.CI <- function(filename, 
							 tstart, tend, 
							 prm.init, 
							 CIwidth, 
							 n.MC, 
							 tstep=1, min.length = 5) {
	# Plots the values of fitted parameters
	# as the number of data increases
	
	tend2 <- tstart+min.length
	i <- 1
	prm.fit <- list() 
	while (tend2 <= tend) {
		dat <- read.inc(filename, tstart,tend2)
		prm.fit[[i]] <- estimate.CI(dat = dat,
									CIwidth = 0.95, 
									n.MC = 10, 
									prm.init = prm.init)
		tend2 <- tend2 + tstep
		i <- i + 1
	}
	
	M <- unlist(prm.fit[[1]])
	for(i in 2:length(prm.fit)){
		M <- rbind(M,unlist(prm.fit[[i]]))
	}
	colnames(M) <- c("r","p","r.lo","r.hi","p.lo","p.hi")
	
	par(mfrow=c(1,1))
	alphavec <- 0.3 #pmax(0.2,1:nrow(M)/(nrow(M)))
	title <- paste0("Estimates (r,p) with growing data set\n from ",
					min.length," to ",tend-tstart," points (",filename,")")
	
	plot(x=M[,'r'], y=M[,'p'],
		 main = title,
		 xlim = c(0,max(M[,'r.hi'])),
		 ylim = c(0,1),
		 xlab="r",ylab="p",typ="o",lty=2,
		 col=rgb(0,0,0,0.3), cex=2,
		 pch=16, las=1)
# 	rect(xleft = M[,'r.lo'],ybottom = M[,'p.lo'],
# 		 xright = M[,'r.hi'],ytop = M[,'p.hi'],
# 		 density = 30,
# 		 col = rgb(0,0,0,0.1),border = rgb(0,0,0,0))
	lwdseg <- 2
	segments(x0=M[,'r'],x1=M[,'r'],y0=M[,'p.lo'],y1=M[,'p.hi'],col=rgb(0,0,0,alphavec),lwd = lwdseg)
	segments(x0=M[,'r.lo'],x1=M[,'r.hi'],y0=M[,'p'],y1=M[,'p'],col=rgb(0,0,0,alphavec),lwd = lwdseg)
	
 	text(M[1,1],M[1,2],labels = paste0("start(",min.length,")"),pos=3,cex=0.7)
#  	rect(xleft = M[nrow(M),'r.lo'],ybottom = M[nrow(M),'p.lo'],
#  		 xright = M[nrow(M),'r.hi'],ytop = M[nrow(M),'p.hi'],
#  		 density = 0,
#  		 lwd = 3,
#  		 col = rgb(1,0,0))
 	segments(x0=M[nrow(M),'r'],x1=M[nrow(M),'r'],y0=M[nrow(M),'p.lo'],y1=M[nrow(M),'p.hi'],col="red",lwd = 4)
 	segments(x0=M[nrow(M),'r.lo'],x1=M[nrow(M),'r.hi'],y0=M[nrow(M),'p'],y1=M[nrow(M),'p'],col="red",lwd = 4)
 	points(x=M[nrow(M),'r'],y=M[nrow(M),'p'],col=rgb(1,0,0,0.5),cex = 2, pch=15)
	# grid(lty=1)
}