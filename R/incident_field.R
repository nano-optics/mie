
incident_field <- function(wavelength,nMedium,r0,theta){
  nL <- length(wavelength)
  nT <- length(theta)
  k <- 2*pi*nMedium / wavelength
  kr <- r0 * outer(k, cos(theta))
  # exp(ikM z) is [L x T], obtained by matrix product of [L x 1] by [1 x T]
  phasefact <- exp(1i * kr)
  matsin <- matrix(sin(theta), ncol=nT, nrow=nL, byrow = TRUE) 
  matcos <- matrix(cos(theta), ncol=nT, nrow=nL, byrow = TRUE) 
  phasefact * matsin 
  # Results are all [L x T] matrices
  # They result from Eq. H.16 for e_x and H.76 for E_inc
  list(Ecr= phasefact * matsin ,
       Ect = phasefact * matcos,
       Esf = - phasefact)
}