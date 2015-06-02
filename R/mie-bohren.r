## internal wrapper
bhcoat <- function(x, y, rrcore, rrshell, lmax=10){

  call <- .Fortran("bhcoat", 
                   as.double(x),
                   as.double(y),
                   as.complex(rrcore),
                   as.complex(rrshell),
                   as.integer(lmax),
                   qext=double(1),
                   qsca=double(1),
                   qback=double(1))
  
  c(qext = call$qext, qsca = call$qsca, qabs = call$qext - call$qsca)
}


##' Far-field cross-sections
##'
##' Coated sphere illuminated by a plane wave
##' @title mie_bh
##' @param wavelength real vector
##' @param epsilon.core complex vector
##' @param epsilon.coating complex vector
##' @param radius scalar
##' @param thickness scalar
##' @param medium scalar, refractive index of surrounding medium
##' @param lmax truncation order
##' @param efficiency logical, scale by geometrical cross-sections
##' @return data.frame
##' @author Baptiste Auguie
##' @family user
##' @export
##' @examples
##' gold <- epsAu(seq(400, 800))
##' coated <- with(gold, mie_bh(wavelength, epsilon, radius=50, medium=1.33, 
##' efficiency=FALSE))
##' bare <- with(gold, mie(wavelength, epsilon, radius=50, medium=1.33, 
##' efficiency=FALSE))
##' matplot(coated$wavelength, coated[, -1], type="l", lty=1,
##'         xlab=expression(lambda/nm), ylab=expression(sigma/nm^2))
##' matlines(bare$wavelength, bare[, -1], type="l", lty=2)
##' legend("topright", c(names(coated)[-1], "bare"), col=1:3, lty=c(1,1,1,2))
mie_bh <- function(wavelength, epsilon.core, epsilon.coating = medium^2,
                   radius, thickness = 0, medium = 1,
                   lmax=ceiling(2 + max(y) + 4 * max(y)^(1/3)),
                   efficiency = TRUE){

  k0 <- 2*pi / wavelength
  k <- k0*medium
  x <- k*radius
  y <- k*(radius+thickness)
  
  rrcore <- sqrt(epsilon.core)/medium
  rrshell <- sqrt(epsilon.coating)/medium
  
  tmp <- do.call(rbind, mapply(bhcoat, x=x, y=y,
                               rrcore=rrcore, rrshell=rrshell, lmax=lmax, SIMPLIFY=FALSE))
  if(!efficiency) tmp <- tmp * (pi*(radius+thickness)^2)
  res <- data.frame(cbind(wavelength, tmp))
  names(res) <- c("wavelength", "extinction", "scattering", "absorption")
  
  res
}
