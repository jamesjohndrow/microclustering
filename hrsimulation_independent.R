rm(list=ls(all=T))
setwd('~/Dropbox (Personal)//Projects/Microclustering/Code')
require('dga')
data(graphs5)
require('LCMCR')
source('hrfuncs.R')
require('lineprof')
require('shiny')
require(Rcapture)
# require('compiler')
# enableJIT(3)


n <- 5000 # unique *observed* entities if samptype!= mult, else unique total entities
samptype <- 'mult'

p <- 5
a <- 1
a2 <- 3
b <- 1
k <- 1

nsim <- 250

sigmas <- c(2/n,1/n,2/(3*n),1/(2*n),1/(3*n),1/(4*n),1/(5*n),1/(10*n))
sigmastrs <- c('2_n','1_n','2_3n','1_2n','1_3n','1_4n','1_5n','1_10n')
nsig <- length(sigmas)

for (ct in 1:nsig) {
  sigma <- sigmas[ct]
  sigmastr <- sigmastrs[ct]
  #  sigma <- 2/(n)
  #  sigmastr <- '2_n'
  
  Nhats <- matrix(0,2^p-1,nsim)
  Ns <- matrix(0,2^p-1,nsim)
  N0s <- matrix(0,nsim,1)
  NR <- matrix(0,nsim,1)
  TR <- matrix(0,nsim,1)
  NTOT <- matrix(0,1,nsim)
  NCI <- matrix(0, nsim, 2)
  ## sample for Dunson-Xing model
  tmp <- dxparam(p,a,b,k)
  B <- tmp$B
  
  
  for (s in 1:nsim) {
    if ((s%%10)==0) {
      print(paste(sigmastr,s))
    }
    rsim <- rlink(p,a,b,k,n,sigma,samptype, a2=a2) # samples individuals
    rsim$Nhat[is.na(rsim$Nhat)] <- 0
    
    Nhats[,s] <- rsim$Nhat #linked/estimated capture histories (no missing cell)
    Ns[,s] <- rsim$N # true capture histories
    N0s[s] <- rsim$N0 # number not captured on any list
    print(paste("N0: ", rsim$N0))
    NR[s] <- rsim$Nright #number correctly assigned to entity 
    TR[s] <- rsim$totrec # total number of records
    
## do we still need to check that every list had at least one capture?    
    dat <- rsim$icap
    
    ##remove NAs
    dat <- dat[!is.na(dat[,2]),-1]

    rcap <- closedpCI.t(dat, m = "Mt")
    NTOT[s] <-rcap$results[1] 
    NCI[s,] <-rcap$CI[2:3]
   }
  
  
  mNtot <- mean(NTOT)
  Ntrue <- apply(Ns,2,sum) + N0s
  mse.Ntot <- mean((NTOT-t(Ntrue))^2)
  cov.Ntot <- mean((NCI[,1] < Ntrue) & (NCI[,2] > Ntrue))
  
  Nhats[is.na(Nhats)] <- 0
  mse.Nhat <- mean((Nhats-Ns)^2,na.rm=T)
  per.right <- sum(NR)/sum(TR)
  
  save.image(paste('Outputs/hrsim_independent_',samptype,'_',sigmastr,'.RData',sep=''))
}



