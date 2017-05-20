library(Bessel)


ricatti_bessel <- function(rho, nmax){
  
  rho <- as.complex(rho)
  sq <- sqrt((pi/2)*rho)
  bj <- sq * BesselJ(z = rho, nu = 0.5, nSeq = nmax+1)
  bh <- sq * BesselH(z = rho, nu = 0.5, nSeq = nmax+1, m = 1)
  
  psi <- bj[, -1L]
  xi <- bh[, -1L]
  norho <- outer(1/rho, seq_len(nmax))
  
  dpsi <- bj[, -(nmax+1)] - norho * psi
  dxi <- bh[, -(nmax+1)] - norho * xi
  
  list(psi=psi, xi=xi, dpsi=dpsi, dxi=dxi)
}

rho <- 2*pi/seq(500, 800, by=100) * 50 * sqrt(-10+1i)
nmax <- 10
ricatti_bessel(rho, nmax)

susceptibility <- function(s, x, nmax){
  
  smat <- matrix(s, ncol=nmax, nrow=length(x), byrow=FALSE)
  z <- s*x
  rbz <- ricatti_bessel(z, nmax)
  rbx <- ricatti_bessel(x, nmax)
  PP1 <- rbz[["psi"]] * rbx[["dpsi"]];
  PP2 <- rbx[["psi"]] * rbz[["dpsi"]];
  PP3 <- rbz[["psi"]] * rbx[["dxi"]];
  PP4 <- rbx[["xi"]] * rbz[["dpsi"]];
  
  G_numerator   <-  - PP1 + smat * PP2
  D_numerator   <-    PP2 - smat * PP1
  B_denominator <-  - PP4 + smat * PP3
  A_denominator <-    PP3 - smat * PP4
  
  list(G = G_numerator / A_denominator,
       D = D_numerator / B_denominator,
       A = 1i * smat   / A_denominator,
       B = 1i * smat   / B_denominator)
}

efficiencies <- function(x, GD, mode=c("EM", "Magnetic", "Electric"), order = NULL){
  
  mode <- match.arg(mode)
  
  nmax <- NCOL(GD$G)
  nvec <- seq.int(nmax)
  nvec2 <- 2 * nvec + 1
  if(all(is.numeric(order)) && all(is.finite(order))) {
    nvec2[-order] <- 0
  }
  
  G2 <- Mod(GD$G)^2
  D2 <- Mod(GD$D)^2
  
  GR <- Re(GD$G)
  DR <- Re(GD$D)
  
  if(mode == "EM"){
    scatmat <- G2 + D2
    ext_coeff <- GR + DR
  } else {
    if(mode == "Electric"){
      scatmat <- D2
      ext_coeff <- DR
    } else {
      scatmat <- G2 
      ext_coeff <- GR
    }
  }
  
  Qsca <- 2 / x^2 * scatmat %*% nvec2
  Qext <- - 2 / x^2 * ext_coeff %*% nvec2
  Qabs <- Qext - Qsca
  
  cbind(Qext = Qext, Qsca = Qsca, Qabs = Qabs)
}

susceptibilities <- function(ls, lx, n_max){
  
  n_w <- max(sapply(lx, length)) # wavelengths
  n_lay <- length(lx) # layers
  
  # browser()
  lz <- Map("*", ls, lx)
  
  lrbz <- lapply(lz, ricatti_bessel, nmax = n_max)
  lrbx <- lapply(lx, ricatti_bessel, nmax = n_max)
  
  # allocate results
  init <- matrix(0+0i, n_w, n_max)
  GD <- replicate(n_lay, list(Gamma = init, Delta = init, A = init, B = init), 
                  simplify=FALSE)
  
  for (kk in seq_len(n_lay)){
    
    if(kk==1){
      
      # first layer, coeffs at 0
      Gamkm1 = init
      Delkm1 = init
      
    } else {
      
      Gamkm1 = GD[[kk-1]][['Gamma']]
      Delkm1 = GD[[kk-1]][['Delta']]
      
    }
    
    
    smat <- matrix(ls[[kk]], ncol=n_max, nrow=n_w)
    
    # auxiliary functions for Gamma and A
    
    PP1 <- lrbz[[kk]][["psi"]] + Gamkm1 * lrbz[[kk]][["xi"]]
    PP2 <- lrbz[[kk]][["dpsi"]] + Gamkm1 * lrbz[[kk]][["dxi"]]
    
    Nnk <- lrbx[[kk]][["dpsi"]] * PP1 - smat * lrbx[[kk]][["psi"]] * PP2
    Dnk <- lrbx[[kk]][["xi"]] * PP2 - lrbx[[kk]][["dxi"]] * PP1 
    
    GD[[kk]][['Gamma']] <- Nnk / Dnk
    GD[[kk]][['A']] = -1i*smat / Dnk
    
    # same logic for Delta and B
    
    PP1 <- lrbz[[kk]][["psi"]] + Delkm1 * lrbz[[kk]][["xi"]]
    PP2 <- lrbz[[kk]][["dpsi"]] + Delkm1 * lrbz[[kk]][["dxi"]]
    
    Nnk <- lrbx[[kk]][["psi"]] * PP2 - smat * lrbx[[kk]][["dpsi"]] * PP1
    Dnk <- smat * lrbx[[kk]][["dxi"]] * PP1 - lrbx[[kk]][["xi"]] * PP2 
    
    GD[[kk]][['Delta']] <- Nnk / Dnk
    GD[[kk]][['B']] = 1i*smat / Dnk
    
  }
  
  GD
}

incident_PWE <- function(n_max){
  
  nn <- seq_len(n_max)
  bn1 <-  1i^(nn+1) * sqrt(pi*(2*nn+1))
  list(bn1 = bn1, an1 =bn1)
}

mie_ml <- function(wavelength, epsilon, radii, 
                nmax = 10,
                efficiency = FALSE, mode=c("EM", "Magnetic", "Electric"),
                order = Inf){
  
  mode <- match.arg(mode)
  n_lay <- length(radii)
  # print(n_lay)
  medium <- epsilon[[n_lay + 1]]
  
  # eps_k / eps_k+1
  ls <- Map(function(e1,e2) sqrt(e1)  / sqrt(e2), epsilon[-(n_lay+1)], epsilon[-1])
  # 2pi/lambda * sqrt(eps+) * a
  lx <- Map(function(r, e) 2 * pi / wavelength * sqrt(e) * r, radii, epsilon[-1])

  GD <- susceptibilities(ls, lx, nmax)
  
  Einc <- incident_PWE(n_max)
  
  # % calculate incident PW coefficient an1 and bn1 [1 x nNmax]
  # stIncEabn1=PweIncEabn1(nNmax);
  # an1mat=repmat(stIncEabn1.an1,length(lambda),1); % [L x nNmax]
  # bn1mat=repmat(stIncEabn1.bn1,length(lambda),1); % [L x nNmax]
  
  # % calculate Mie coefficients using a downward recurrence (pp. 623,624)
  # % each cell of coeff corresponds to a given kk=0..nK
  
  
  # allocate results
  init <- rep(0+0i, n_w)
  abcd <- replicate(n_lay + 1, 
                    list(gamma = init, delta = init, alpha = init, beta = init), 
                  simplify=FALSE)
  
  # Initialize recurrence 
  abcd[[n_lay+1]][["alpha"]] <- Einc[["an1"]]
  abcd[[n_lay+1]][["beta"]] <- Einc[["bn1"]]
  
  for (kk in seq(from=n_lay, to=1, by=-1)){
    
    abcd[[kk+1]][["gamma"]] <- GD[[kk]][["Gamma"]] * abcd[[kk+1]][["alpha"]]
    abcd[[kk+1]][["delta"]] <- GD[[kk]][["Delta"]] * abcd[[kk+1]][["beta"]]
    
    abcd[[kk]][["alpha"]] <- GD[[kk]][["A"]] * abcd[[kk+1]][["alpha"]]
    abcd[[kk]][["beta"]] <- GD[[kk]][["B"]] * abcd[[kk+1]][["beta"]]
    
  }
  
  # % get theta-dep functions to avoid repeated computations
  # theta=linspace(0,pi,nNbTheta); % row [1 x T]
  # stPinTaun=PwePinTaun(stM.nNmax,transpose(theta)); % fields are [T x nNmax]
  # 
  # % Spherical averages and theta dependence on the largest sphere surface (outside)
  # stEsurf=PweSurfProperties(stM,stM.a,nNbTheta,stPinTaun);
  
  
  
  
  # use last one for efficiencies
  Q <- efficiencies(lx[[n_lay]], GD[[n_lay]], mode=mode, order=order)
  if(!efficiency) Q <- Q * (pi*radii[[n_lay]]^2)
  results <- data.frame(wavelength, Q)
  names(results) <- c("wavelength", "extinction", "scattering", "absorption")
  invisible(results)
}

library(dielectric)
gold <- epsAg(seq(300, 800))
a <- 30
b <- 34
c <- 35

bare <- mie(gold$wavelength, gold$epsilon, radius=a, medium=1.33, efficiency=FALSE)

leps <- list(gold$epsilon, 1.33^2, 1.5^2, 1.33^2)
# leps <- list(1.5^2, gold$epsilon,  1.33^2)
la <- list(a,b,c)
coated <- mie_ml(gold$wavelength, leps, radii=la, efficiency=FALSE)

matplot(bare$wavelength, bare[, -1], type="l", lty=1,
        xlab=expression(lambda/mu*m), ylab=expression(sigma/mu*m^2))
matlines(coated$wavelength, coated[, -1],  type="l",  lty=2)

legend("topright", c(names(bare)[-1], "coated"), col=1:3, lty=c(1,1,1,2))

