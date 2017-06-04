require('climtrends')

dxparam <- function(p,a,b,k, a2=a) { #change this to allow multiple values of a so not beta(a,a)
  nu = rgamma(k,b,1)
  nu = nu/sum(nu)
  lam = array(0,dim=c(2^p,2,k))
  P = matrix(0,2^p,1)
  B = matrix(0,2^p,p)
  for (j in 1:p) {
    B[,j] = rep(sort(rep(c(0,1),2^(j-1))),2^(p-j))
  }
  
  for (j in 1:k) {
    lam = rbeta(p,a,a2)
    pr = t(matrix(rep(c(lam),2^p),p,2^p))
    pr = apply(pr*B+(1-pr)*(1-B),1,prod)
    P = P + nu[j]*pr
  }
  return(list(P=P,B=B))
}

rlink <- function(p,a,b,k,n,sigma,samptype, a2=a) {
  Pi <- dxparam(p,a,b,k, a2=a2)
  P <- Pi$P
  B <- Pi$B
  B <- B[2:(2^p),]
  
  if (samptype=='mult') {
    N <- rmultinom(1,n,P)
    N0 <- N[1]
    N <- N[2:(2^p)]
    nob <- sum(N)
  } else {
    N <- rmultinom(1,n,P[2:(2^p)]/sum(P[2:(2^p)]))
    N0 <- rnbinom(1,n,1-P[1])
    nob <- n
  }
  csN <- cumsum(N)
  csN <- c(0,csN)
  
  pid <- sample.int(n,nob,replace=F)
  theta <- pid/n
  
  nrowR <- 0 
  for (j in 1:(2^p-1)) {
    nrowR <- nrowR + sum(B[j,])*N[j]
  }
  
  R <- matrix(0,nrowR,2+p)
  Nright <- 0
  totrec <- 0
  rloc <- 1
  for (j in 1:(2^p-1)) {
    if (N[j]>0) {
      nrec <- sum(B[j,])
      totrec <- totrec + nrec*N[j]
      B0 <- matrix(0,N[j]*nrec,p)
      whichcap <- which(B[j,]>0)
      for (l in 1:nrec) {
        B0[((l-1)*N[j]+1):(l*N[j]),whichcap[l]] <- 1
      }
      
      y <- matrix(rnorm(nrec*N[j],theta[(csN[j]+1):csN[j+1]],sigma),N[j],nrec)
      idmax <- round(n*y)
      idmax[idmax>n] <- idmax[idmax>n]%%n
      idmax[idmax==0] <- n
      idmax[idmax<0] <- (n+idmax[idmax<0])%%n
      isright <- (idmax==matrix(rep(pid[(csN[j]+1):csN[j+1]],nrec),N[j],nrec))
      Nright <- Nright+sum(c(isright*1))
      R[rloc:(rloc+N[j]*nrec-1),1] <- c(idmax)
      #R[rloc:(rloc+N[j]*nrec-1),2:(p+1)] <- t(matrix(rep(B[j,],N[j]*nrec),p,N[j]*nrec))
      R[rloc:(rloc+N[j]*nrec-1),2:(p+1)] <- B0
      R[rloc:(rloc+N[j]*nrec-1),p+2] <- j
      rloc <- rloc+N[j]*nrec
    }
  }
  
  
  M <- data.frame(seq(n))
  names(M) <- 'id'
  for (j in 1:p) {
    ctj <- tapply(R[,j+1],R[,1],max)
    ctj <- data.frame(as.numeric(rownames(ctj)),ctj)
    names(ctj) <- c('id',paste('iscap',j,sep='')) 
    M <- merge(M,ctj,by='id',all.x=T)
  }
  icap <- M
  
  Bid <- cbind(B,seq(2^p-1))
  Bid <- data.frame(Bid)
  names(Bid) <- c(paste('iscap',seq(p),sep=''),'bid')
  
  M <- merge(M,Bid,by=c(paste('iscap',seq(p),sep='')),all.x=T)
  Nhat <- table(M$bid)
  
  Nhat <- data.frame(Nhat)
  names(Nhat) <- c('bid','nhat')
  Nhat$nhat <- as.numeric(Nhat$nhat)
  Nhat2 <- data.frame(seq(2^p-1))
  names(Nhat2) <- 'bid'
  Nhat2 <- merge(Nhat2,Nhat,by='bid',all.x=T)
  Nhat <- as.matrix(Nhat2$nhat)
    
  return(list(Nhat=Nhat,N=N,N0=N0,Nright=Nright,totrec=totrec,icap=icap))
}



rlink_samp <- function(p,a,b,k,n,sigma,nsamp,ncheck,samptype,a2=a) {
  Pi <- dxparam(p,a,b,k,a2=a2)
  P <- Pi$P
  B <- Pi$B
  B <- B[2:(2^p),]
  
  if (samptype=='mult') {
    N <- rmultinom(1,n,P)
    N0 <- N[1]
    N <- N[2:(2^p)]
    nob <- sum(N)
  } else {
    N <- rmultinom(1,n,P[2:(2^p)]/sum(P[2:(2^p)]))
    N0 <- rnbinom(1,n,1-P[1])
    nob <- n
  }
  csN <- cumsum(N)
  csN <- c(0,csN)
  
  pid <- sample.int(n,nob,replace=F)
  theta <- pid/n
  #theta.mat <- matrix(rep(theta,p),n,p)
  theta.ord <- seq(n)/n
  
  checkmid <- (ncheck-1)/2+1
  theta.check <- matrix(0,nob,ncheck)
  for (j in 1:ncheck) {
    idxt <- (pid-checkmid+j)%%n
    idxt[idxt==0] <- n
    theta.check[,j] <- theta.ord[idxt]
  }
  
  #theta.ord.mat <- matrix(rep(theta.ord,each=n),n,n)
  
  NHAT <- matrix(0,2^p-1,nsamp)
  NRIGHT <- matrix(0,nsamp,1)
  ICAP <- array(NA,dim=c(n,p+1,nsamp))
  
  for (r in 1:nsamp) {
    print(r)
    nrowR <- 0
    for (j in 1:(2^p-1)) {
      nrowR <- nrowR + sum(B[j,])*N[j]
    }
    #Z <- matrix(rnorm(n*p,0,sigma),n,p)

    R <- matrix(0,nrowR,2+p)
    Nright <- 0
    totrec <- 0
    rloc <- 1
    for (j in 1:(2^p-1)) {
      if (N[j]>0) {
        nrec <- sum(B[j,])
        totrec <- totrec + nrec*N[j]
        B0 <- matrix(0,N[j]*nrec,p)
        whichcap <- which(B[j,]>0)
        for (l in 1:nrec) {
          B0[((l-1)*N[j]+1):(l*N[j]),whichcap[l]] <- 1
        }

        y <- matrix(rnorm(nrec*N[j],theta[(csN[j]+1):csN[j+1]],sigma),N[j],nrec)
        #y <- matrix(Z[(csN[j]+1):csN[j+1],1:nrec] + theta.mat[(csN[j]+1):csN[j+1],1:nrec],N[j],nrec)
        #idmax <- round(n*y)
        idmax <- matrix(0,N[j],nrec)
        for (l in 1:nrec) {
          #rtheta <- t(matrix(rep(theta.ord,N[j]),n,N[j]))
          #rtheta <- theta.ord.mat[1:N[j],]
          #ry <- matrix(rep(y[,l],n),N[j],n)
          ry <- matrix(rep(y[,l],ncheck),N[j],ncheck)
          #ll <- -(rtheta-ry)^2/(2*sigma^2)
          ll <- -(theta.check[(csN[j]+1):csN[j+1],]-y[,l])^2/(2*sigma^2)
          if (N[j] > 1) {
            ll.max <- apply(ll,1,max)
            ll.max <- matrix(rep(ll.max,ncheck),N[j],ncheck)
            ll <- ll-ll.max
            wts <- exp(ll)
            wts.sum <- apply(wts,1,sum)
            wts.sum <- matrix(rep(wts.sum,ncheck),N[j],ncheck)
            wts <- wts/wts.sum
            wts.cs <- t(apply(wts,1,cumsum))
            u <- runif(N[j])
            u <- matrix(rep(u,ncheck),N[j],ncheck)
            isgtr <- u>wts.cs
            sgtr <- apply(isgtr,1,sum)+1
            #idmax[,l] <- sgtr
            persid <- (pid[(csN[j]+1):csN[j+1]]+(sgtr-checkmid))%%n
            persid[persid==0] <- n
            idmax[,l] <- persid
          } else {
            ll.max <- max(ll)
            ll.max <- matrix(rep(ll.max,ncheck),N[j],ncheck)
            ll <- ll-ll.max
            wts <- exp(ll)
            wts.sum <- sum(wts)
            wts.sum <- matrix(rep(wts.sum,ncheck),N[j],ncheck)
            wts <- wts/wts.sum
            wts.cs <- cumsum(wts)
            u <- runif(N[j])
            u <- matrix(rep(u,ncheck),N[j],ncheck)
            isgtr <- u>wts.cs
            sgtr <- sum(isgtr)+1
            #idmax[,l] <- sgtr
            persid <- (pid[(csN[j]+1):csN[j+1]]+(sgtr-checkmid))%%n
            persid[persid==0] <- n
            idmax[,l] <- persid
          }
        }



        idmax[idmax>n] <- idmax[idmax>n]%%n
        idmax[idmax==0] <- n
        idmax[idmax<0] <- (n+idmax[idmax<0])%%n
        isright <- (idmax==matrix(rep(pid[(csN[j]+1):csN[j+1]],nrec),N[j],nrec))
        Nright <- Nright+sum(c(isright*1))
        R[rloc:(rloc+N[j]*nrec-1),1] <- c(idmax)
        #R[rloc:(rloc+N[j]*nrec-1),2:(p+1)] <- t(matrix(rep(B[j,],N[j]*nrec),p,N[j]*nrec))
        R[rloc:(rloc+N[j]*nrec-1),2:(p+1)] <- B0
        R[rloc:(rloc+N[j]*nrec-1),p+2] <- j
        rloc <- rloc+N[j]*nrec
      }
    }

    M <- data.frame(seq(n))
    names(M) <- 'id'
    for (j in 1:p) {
      ctj <- tapply(R[,j+1],R[,1],max)
      ctj <- data.frame(as.numeric(rownames(ctj)),ctj)
      names(ctj) <- c('id',paste('iscap',j,sep=''))
      M <- merge(M,ctj,by='id',all.x=T)
    }
    icap <- M

    Bid <- cbind(B,seq(2^p-1))
    Bid <- data.frame(Bid)
    names(Bid) <- c(paste('iscap',seq(p),sep=''),'bid')

    M <- merge(M,Bid,by=c(paste('iscap',seq(p),sep='')),all.x=T)
    Nhat <- table(M$bid)

    Nhat <- data.frame(Nhat)
    names(Nhat) <- c('bid','nhat')
    Nhat$nhat <- as.numeric(Nhat$nhat)
    Nhat2 <- data.frame(seq(2^p-1))
    names(Nhat2) <- 'bid'
    Nhat2 <- merge(Nhat2,Nhat,by='bid',all.x=T)
    Nhat <- as.matrix(Nhat2$nhat)
    

    NHAT[,r] <- Nhat
    NRIGHT[r] <- Nright
    ICAP[,,r] <- as.matrix(icap)
  }
  
  return(list(NHAT=NHAT,N=N,N0=N0,NRIGHT=NRIGHT,totrec=totrec,ICAP=ICAP))
}



samp_index <- function(theta,n,N,yl,sigma) {
  rtheta <- t(matrix(rep(theta,N[j]),n,N[j]))
  ry <- matrix(rep(yl,n),N[j],n)
  ll <- -(rtheta-ry)^2/(2*sigma^2)
  ll.max <- apply(ll,1,max)
  ll.max <- matrix(rep(ll.max,n),N[j],n)
  ll <- ll-ll.max
  wts <- exp(ll)
  wts.sum <- apply(wts,1,sum)
  wts.sum <- matrix(rep(wts.sum,n),N[j],n)
  wts <- wts/wts.sum
  wts.cs <- t(apply(wts,1,cumsum))
  u <- runif(N[j])
  u <- matrix(rep(u,n),N[j],n)
  isgtr <- u>wts.cs
  sgtr <- apply(isgtr,1,sum)+1
  idmax[,l] <- sgtr
}






