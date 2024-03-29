
##' Far-field cross-sections
##'
##' Homogeneous sphere illuminated by a plane wave
##' @title mie
##' @param wavelength real vector
##' @param epsilon complex vector
##' @param radius scalar
##' @param medium scalar, refractive index of surrounding medium
##' @param n_max truncation order
##' @param efficiency logical, scale by geometrical cross-sections
##' @param mode type of mode
##' @param order order of multipoles
##' @return data.frame
##' @author Baptiste Auguie
##' @family user
##' @export
##' @examples 
##' gold <- epsAu(seq(400, 800))
##' cross_sections <- with(gold, mie(wavelength, epsilon, radius=50, medium=1.33, efficiency=FALSE))
##' matplot(cross_sections$wavelength, cross_sections[, -1], type="l", lty=1,
##'         xlab=expression(lambda/mu*m), ylab=expression(sigma/mu*m^2))
##' legend("topright", names(cross_sections)[-1], col=1:3, lty=1)
mie <- function(wavelength, epsilon, radius, medium = 1.0,
                n_max=ceiling(2 + max(x) + 4 * max(x)^(1/3)),
                efficiency = FALSE, mode=c("EM", "Magnetic", "Electric"),
                order = Inf){
  
  mode <- match.arg(mode)
  
  s <- sqrt(epsilon) / medium
  x <- 2 * pi / wavelength * medium * radius
  
  ## lazy evaluation rules.. default n_max evaluated now
  coeffs <- susceptibility(s, x, n_max)
  Q <- efficiencies(x, coeffs, mode=mode, order=order)
  if(!efficiency) Q <- Q * (pi*radius^2)
  results <- data.frame(wavelength, Q)
  names(results) <- c("wavelength", "extinction", "scattering", "absorption")
  invisible(results)
}


##' Average Mloc
##'
##' Enhancement factor averaged over the surface 
##' @title average_Mloc
##' @param wavelength real vector
##' @param epsilon complex vector
##' @param radius scalar
##' @param medium scalar, refractive index of surrounding medium
##' @param n_max truncation order
##' @return data.frame with wavelength and Mloc
##' @author Baptiste Auguie
##' @export
average_Mloc <- function(wavelength, epsilon, radius, medium = 1.0, n_max=10){
  
  s <- sqrt(epsilon) / medium
  x <- 2 * pi / wavelength * medium * radius
  
  rbx <- ricatti_bessel(x, n_max)
  GD <- susceptibility(s, x, n_max)
  
  tmp1 <- rbx[["psi"]] + GD[["G"]] * rbx[["xi"]]
  tmp2 <- rbx[["dpsi"]] + GD[["D"]] * rbx[["dxi"]]
  summat1 <- abs(tmp1)^2 + abs(tmp2)^2
  
  nlen <- seq_len(n_max)
  cc1 <- t(2*nlen + 1)
  cc2 <- cc1*nlen*(nlen+1)
  
  tmp1 <- rbx[["psi"]] + GD[["D"]] * rbx[["xi"]] 
  summat2 <- abs(tmp1)^2
  
  # From Eq. H.79
  data.frame(wavelength=wavelength, Mloc = 1/(2*x^2) * tcrossprod(summat1, cc1) + 1/(2*x^4) * tcrossprod(summat2, cc2))
}


##' Efficiencies
##'
##' Calculates the far-field efficiencies for plane-wave illumination
##' @title efficiencies
##' @param x real vector, size parameter
##' @param GD list with Gamma, Delta, A, B
##' @param mode type of mode
##' @param order order of multipoles
##' @return matrix of Qext, Qsca, Qabs
##' @author Baptiste Auguie
##' @export
efficiencies <- function(x, GD, mode=c("EM", "Magnetic", "Electric"), order = NULL){
  
  mode <- match.arg(mode)
  
  n_max <- NCOL(GD$G)
  nvec <- seq_len(n_max)
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


##' Riccati-Bessel function psi and its derivative
##'
##' Obtained from BesselJ, converted to spherical Bessel, and scaled
##' @title psi
##' @param rho complex vector, argument
##' @param n_max integer, maximum order
##' @return a list with psi_n and psi'_n
##' @author Baptiste Auguie
##' @export
psi <- function(rho, n_max){

  nvec <- seq.int(n_max)
  nmat <- matrix(nvec, ncol=n_max, nrow=length(rho), byrow=TRUE)
  rhomat <- matrix(rho, ncol=n_max, nrow=length(rho), byrow=FALSE)
  psi <- sqrt(rho * pi/2) * BesselJ(rho, 1/2, expon.scaled = FALSE, nSeq = n_max+1)
  psip <- psi[ , nvec] - nmat * psi[ , nvec + 1] / rhomat
  
  list(psi = psi[ , -1], psip = psip)
  
}


##' Riccati-Bessel function xi and its derivative
##'
##' Obtained from BesselH (Hankel function), converted to spherical Hankel, and scaled
##' @title xi
##' @param rho complex vector, argument
##' @param n_max integer, maximum order
##' @return a list with psi_n and psi'_n
##' @author Baptiste Auguie
##' @export
xi <- function(rho, n_max){
  
  nvec <- seq.int(n_max)
  nmat <- matrix(nvec, ncol=n_max, nrow=length(rho), byrow=TRUE)
  rhomat <- matrix(rho, ncol=n_max, nrow=length(rho), byrow=FALSE)
  xi <- sqrt(rho * pi/2) * BesselH(1, rho+0i, 1/2, expon.scaled = FALSE, nSeq = n_max+1)
  xip <- xi[ , nvec] - nmat * xi[ , nvec + 1] / rhomat
  
  list(xi = xi[ , -1], xip = xip)

  
}


ricatti_bessel <- function(rho, n_max){
  
  rho <- as.complex(rho)
  sq <- sqrt((pi/2)*rho)
  bj <- sq * BesselJ(z = rho, nu = 0.5, nSeq = n_max+1)
  bh <- sq * BesselH(z = rho, nu = 0.5, nSeq = n_max+1, m = 1)
  
  psi <- bj[, -1L]
  xi <- bh[, -1L]
  norho <- outer(1/rho, seq_len(n_max))
  
  dpsi <- bj[, -(n_max+1)] - norho * psi
  dxi <- bh[, -(n_max+1)] - norho * xi
  
  list(psi=psi, xi=xi, dpsi=dpsi, dxi=dxi)
}



##' Generalised susceptibility for the Mie theory
##'
##' Corresponds to the usual coefficients a_n, b_n, c_n, d_n
##' @title susceptibility
##' @param n_max integer, maximum order
##' @param s complex vector, relative refractive index
##' @param x real vector, size parameter
##' @return list with Gamma, Delta, A, B
##' @author Baptiste Auguie
##' @export
susceptibility <- function(s, x, n_max){
  
  smat <- matrix(s, ncol=n_max, nrow=length(x), byrow=FALSE)
  z <- s*x
  rbz <- ricatti_bessel(z, n_max)
  rbx <- ricatti_bessel(x, n_max)
  
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
