\name{mie}
\alias{mie}
\title{mie}
\usage{
  mie(wavelength, epsilon, radius, medium = 1,
    nmax = ceiling(2 + max(x) + 4 * max(x)^(1/3)),
    efficiency = TRUE)
}
\arguments{
  \item{wavelength}{real vector}

  \item{epsilon}{complex vector}

  \item{radius}{scalar}

  \item{medium}{scalar, refractive index of surrounding
  medium}

  \item{nmax}{truncation order}

  \item{efficiency}{logical, scale by geometrical
  cross-sections}
}
\value{
  data.frame
}
\description{
  Far-field cross-sections
}
\details{
  Homogeneous sphere illuminated by a plane wave
}
\examples{
gold <- epsAu(seq(400, 800))
cross_sections <- with(gold, mie(wavelength, epsilon, radius=0.05, medium=1.33, efficiency=FALSE))
matplot(cross_sections$wavelength, cross_sections[, -1], type="l", lty=1,
        xlab=expression(lambda/mu*m), ylab=expression(sigma/mu*m^2))
legend("topright", names(cross_sections)[-1], col=1:3, lty=1)
}
\author{
  Baptiste Auguie
}
\seealso{
  Other user: \code{\link{mie_bh}},
  \code{\link{mie-package}}
}

