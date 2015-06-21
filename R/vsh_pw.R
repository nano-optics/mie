vsh_pw <- function(wavelength, epsilon, cn1, dn1, r0, theta, type=c("j","h"), pitau){
  
  NT <- length(theta)
  NL <- length(wavelength)
  dummy <- matrix(0, nrow=NL, ncol=NT)
  E <- list(Er = dummy, Et = dummy, Ep = dummy)
  
  # E = E_{cr} cos(phi) e_r + E_{ct} cos(phi) e_theta + E_{sf} sin(phi) e_phi
  return(E)
}