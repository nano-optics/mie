\name{epsAg}
\alias{epsAg}
\title{epsAg}
\usage{
  epsAg(lambda, epsilon.infty = 4, lambda.p = 282,
    mu.p = 17000)
}
\arguments{
  \item{lambda}{wavelength in nm}

  \item{epsilon.infty}{background dielectric constant}

  \item{lambda.p}{plasma wavelength}

  \item{mu.p}{damping constant}
}
\value{
  data.frame
}
\description{
  permittivity silver
}
\details{
  analytical dielectric function of Silver (Drude model)
}
\examples{
wvl <- seq(300, 900)
silver <- planar:::epsAg(wvl)

matplot(silver$wavelength, cbind(Re(silver$epsilon), Im(silver$epsilon)), t="l", lty=1)
if(require(dielectric)) { data(AgPalik)
matpoints(AgPalik$wavelength, cbind(Re(AgPalik$epsilon), Im(AgPalik$epsilon))) }
}
\author{
  baptiste Auguie
}
\references{
  Principles of surface-enhanced Raman spectroscopy and
  related plasmonic effects Eric C. Le Ru and Pablo G.
  Etchegoin, published by Elsevier, Amsterdam (2009).
}

