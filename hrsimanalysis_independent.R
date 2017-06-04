rm(list=ls(all=T))
#setwd('~/Dropbox/Projects/Microclustering/Code')
setwd('~/Dropbox (Personal)/Projects/Microclustering/Code')
require('ggplot2')
require('reshape2')
require('LCMCR')

samptype <- 'mult'
sigmanames <- c('1_10n','1_5n','1_4n','1_3n','1_2n','2_3n','1_n','2_n')
nsig <- length(sigmanames)
PRIGHT <- matrix(0,nsig,1)
COV <- matrix(0,nsig,1)
MSE.N <- matrix(0,nsig,1)
MSE.IS <- matrix(0,nsig,1)
ESTS <- matrix(0, nsig, 250)

for (sigid in 1:nsig) {
  load(paste('Outputs/hrsim_independent_',samptype,'_',sigmanames[sigid],'.RData',sep=''))
  print(mean((Ns-Nhats)^2))
  PRIGHT[sigid] <- per.right
  COV[sigid] <- cov.Ntot
  MSE.N[sigid] <- mse.Ntot
  MSE.IS[sigid] <- mse.Nhat
  ESTS[sigid,] <- mNtot
}

dat <- data.frame(cbind(PRIGHT,COV,MSE.N,MSE.IS))
names(dat) <- c('per_right','coverage','mse_pop','mse_isect')
dat$c <- c(1/10,1/5,1/4,1/3,1/2,2/3,1,2)

dat <- melt(dat,id='c')

png(paste('Figures/hrsim_results_independent',samptype,'.png',sep=''),width=900,height=600)
ggplot(dat,aes(x=c,y=value)) + geom_point(size=4) + facet_wrap(~variable,scales="free") + theme(text=element_text(size=24))
dev.off()





