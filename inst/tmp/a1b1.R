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

rb1 <- function(x){
  
  j0 <- sin(x)/x
  j1 <- sin(x)/x^2 - cos(x)/x
  y0 <- -cos(x)/x
  y1 <- -cos(x)/x^2 - sin(x)/x
  h0 <- j0 + 1i*y0
  h1 <- j1 + 1i*y1
  
  psi0 <- x*j0
  psi1 <- x*j1
  xi0 <- x*h0
  xi1 <- x*h1
  
  dpsi1 = psi0 - 1/x * psi1
  dxi1 = xi0 - 1/x * xi1
  
  list(psi=psi1, xi=xi1, dpsi=dpsi1, dxi=dxi1)
}

rho <- 2*pi/seq(500, 800, by=100) * 50 * sqrt(-10+1i)
nmax <- 10
ricatti_bessel(rho, 1)
rb1(rho)
ricatti_bessel(rho, 1)$dpsi
  
  
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


# 
# x <- 2*pi/seq(500, 800, by=100)
# s <- 1.5
# susceptibility1(s, x, 10)
# susceptibility(10, s, x)

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

##' Far-field cross-sections
##'
##' Homogeneous sphere illuminated by a plane wave
##' @title mie
##' @param wavelength real vector
##' @param epsilon complex vector
##' @param radius scalar
##' @param medium scalar, refractive index of surrounding medium
##' @param nmax truncation order
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
##'
##'gold <- epsAu(seq(200, 1500))
##'library(ggplot2)
##'
##'params <- expand.grid(order = c(1, 2, Inf), mode = c("EM", "Magnetic", "Electric"), stringsAsFactors=FALSE)
##'
##'all <- plyr::mdply(params, mie, wavelength=gold$wavelength, 
##'              epsilon=gold$epsilon, radius=80, medium=1.5,
##'             .progress="text")
##'
##'m <- tidyr::pivot_longer(all, meas = c("extinction", "scattering", "absorption"))
##'
##'ggplot(m) +
##'  facet_grid(mode~variable, scales="free") +
##'  geom_path(aes(wavelength, value, colour = factor(order))) +
##'  scale_linetype_manual(values = c(2, 3, 1)) +
##'  labs(x = expression(wavelength / nm),
##'       y = expression(sigma / nm^2),
##'       colour = "Mode",
##'       linetype = "Order")
##' 
mie <- function(wavelength, epsilon, radius, medium = 1.0,
                nmax=ceiling(2 + max(x) + 4 * max(x)^(1/3)),
                efficiency = FALSE, mode=c("EM", "Magnetic", "Electric"),
                order = Inf){
  
  mode <- match.arg(mode)
  
  s <- sqrt(epsilon) / medium
  x <- 2 * pi / wavelength * medium * radius
  
  ## lazy evaluation rules.. default nmax evaluated now
  coeffs <- susceptibility(s, x, nmax)
  Q <- efficiencies(x, coeffs, mode=mode, order=order)
  if(!efficiency) Q <- Q * (pi*radius^2)
  results <- data.frame(wavelength, Q)
  names(results) <- c("wavelength", "extinction", "scattering", "absorption")
  invisible(results)
}

library(dielectric)
gold <- epsAg(seq(300, 800))
a <- 30
cross_sections <- with(gold, mie(wavelength, epsilon, radius=a, medium=1.33, efficiency=TRUE, order=1))
matplot(cross_sections$wavelength, cross_sections[, -1], type="l", lty=1,
        xlab=expression(lambda/mu*m), ylab=expression(sigma/mu*m^2))
# M <- with(gold, average_Mloc(wavelength, epsilon, radius=a, medium=1.33, nmax=100))
# lines(M$wavelength, M$Mloc/max(M$Mloc)*max(cross_sections[,-1]), lty=2)
legend("topright", c(names(cross_sections)[-1], "<Mloc>"), col=1:3, lty=c(1,1,1,2))

 
gold <- epsAu(seq(200, 1500))
library(ggplot2)

params <- expand.grid(order = c(1, 2, 3, Inf), mode = c("EM", "Magnetic", "Electric"), stringsAsFactors=FALSE)

all <- purrr::pmap_df(params, mie, wavelength=gold$wavelength, 
                   epsilon=gold$epsilon, radius=80, medium=1.5,
                   .id='id')
params$id <- as.character(1:nrow(params))
all <- left_join(params, all, by='id')

m <- tidyr::pivot_longer(all, cols = c("extinction", "scattering", "absorption"))

ggplot(m) +
  facet_grid(mode~name, scales="free") +
  geom_path(aes(wavelength, value, colour = factor(order))) +
  scale_linetype_manual(values = c(2, 3, 1)) +
  labs(x = expression(wavelength / nm),
       y = expression(sigma / nm^2),
       colour = "Mode",
       linetype = "Order")

