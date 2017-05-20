
pitaun <- function(nmax, theta){
  nrows <- length(theta)  
  mu <- cos(theta)
  # Initialize recurrence pi_0=0, pi_1=1
  # pinm1 contains pi_{n-1} n=1..nNmax+1
  pinm1 <- matrix(0, nrows, nmax + 1)
  pinm1[ , 2] <- 1
  # Get pi_2 to pi_nNmax by recurrence (see SPlaC guide)
  # pi_n is pinm1(:,n+1)
  for (n in seq(2, nmax)) {
    pinm1[ ,n+1] <- (2*n-1)/(n-1) * mu * pinm1[ ,n]- n/(n-1)*pinm1[,n-1]
  }
  pin <- pinm1[,-1L]
  # return tau_n matrix
  nvec <- seq_len(nmax)
  nmat <- matrix(nvec+1, ncol=nmax, nrow=nrows, byrow = TRUE) 
  taun <- outer(mu,nvec) * pin - nmat * pinm1[,nvec]
  
  # return pi_n matrix (except n=0)
  list(pin = pin, taun = taun)
}

bessel_Z <- function(rho, n_max, type){
  
  list(Z0, Z1, Z2)
}

mie_nf <- function(wavelength, epsilon, abcd, r, theta){
  
  n_max <- ncol(abcd[[1]])
  pt <- pitaun(n_max, theta)
  
  Ereg <- Etheta(wavelength, epsilon, an, bn, r, theta, pt, "regular")
  Eirr <- Etheta(wavelength, epsilon, cn, dn, r, theta, pt, "irregular")
  
  list(Ecr = Ereg[["Ecr"]] + Eirr[["Ecr"]],  
             Ereg[["Ect"]] + Eirr[["Ect"]],
             Ereg[["Esf"]] + Eirr[["Esf"]])
}

Etheta <- function(wavelength, epsilon, cn, dn, r, theta, pt, type){
  
  rho <- 2*pi* sqrt(epsilon) / wavelength * r
  Z <- bessel_Z(rho, n_max, type)
  
  n_max <- ncol(cn)
  nn <- seq_len(n_max)
  
  cffnr <- sqrt((2*nn+1)/(4*pi))
  mun <- cffnr / (nn*(nn+1))
  
  dn1Z1 <- dn * Z[["Z1"]]
  icn1Z0 <- 1i*cn * Z[["Z0"]]
  dn1Z2 <- dn *  Z[["Z2"]]
  
  for(ii in seq_along(wavelength)){
    
    vecNdep1 <- dn1Z1[ii,] * cffnr
    vecNdep0 <- icn1Z0[ii,] * mun
    vecNdep2 <- dn1Z2[ii,] * mun
    Ersum[ii,] <- pt[["pin"]] %*% vecNdep1
    
    # % for Et and Ef
    tmp1 <- pt[["pin"]] %*% vecNdep0
    tmp2 <- pt[["taun"]] %*% vecNdep2
    Etsum[ii,] <- tmp1 %*% tmp2
   
    tmp1 <- pt[["taun"]] %*% vecNdep0
    tmp2 <- pt[["pin"]] %*% vecNdep2
    Efsum[ii,] <- tmp1 %*% tmp2
  }
  
  
  list(Ecr = -2*sin(theta) * Ersum,
       Ect = -2*Etsum,
       Esf = 2 * Efsum)
}



# 
# stEAllPhiReg=PweEgenThetaAllPhi(lambda,epsilon,stAbcdn1.an1,stAbcdn1.bn1,r0,theta,'j',stPinTaun);
# stEAllPhiIrr=PweEgenThetaAllPhi(lambda,epsilon,stAbcdn1.cn1,stAbcdn1.dn1,r0,theta,'h1',stPinTaun);
#
# % disp 'PweEsurf: Compiling results and averages...'
# % Ecr, Ect, and Esf are [L x T]
# stEsurf.Ecr=GenCheckSum2Mat(stEAllPhiReg.Ecr,stEAllPhiIrr.Ecr,'Ecr','PweEsurf');
# stEsurf.Ect=GenCheckSum2Mat(stEAllPhiReg.Ect,stEAllPhiIrr.Ect,'Ect','PweEsurf');
# stEsurf.Esf=GenCheckSum2Mat(stEAllPhiReg.Esf,stEAllPhiIrr.Esf,'Esf','PweEsurf');
# 
# % computes surface-averages
# stEsurf=PweEFaverages(stEsurf);
