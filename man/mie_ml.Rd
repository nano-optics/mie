% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mie_layer.R
\name{mie_ml}
\alias{mie_ml}
\title{mie_ml}
\usage{
mie_ml(
  wavelength,
  epsilon,
  radii,
  n_max = 10,
  efficiency = FALSE,
  mode = c("EM", "Magnetic", "Electric"),
  order = Inf
)
}
\arguments{
\item{wavelength}{real vector}

\item{epsilon}{list of dielectric functions, from inner to outer}

\item{radii}{concentric radii of each interface, from smaller to larger}

\item{efficiency}{logical, scale by geometrical cross-sections}

\item{mode}{type of mode}

\item{order}{order of multipoles}

\item{medium}{scalar, refractive index of surrounding medium}

\item{nmax}{truncation order}
}
\value{
data.frame
}
\description{
Far-field cross-sections
}
\details{
Multilayered sphere illuminated by a plane wave
}
\seealso{
Other user: 
\code{\link{mie-package}},
\code{\link{mie_approximation}()},
\code{\link{mie_bh}()},
\code{\link{mie}()}
}
\author{
Baptiste Auguie
}
\concept{user}