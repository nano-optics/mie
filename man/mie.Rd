\name{mie}
\alias{mie}
\title{mie}
\usage{
  mie(wavelength, epsilon, radius, medium = 1,
    nmax = ceiling(2 + max(x) + 4 * max(x)^(1/3)),
    efficiency = FALSE,
    mode = c("EM", "Magnetic", "Electric"), order = Inf)
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

  \item{mode}{type of mode}

  \item{order}{order of multipoles}
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

gold <- epsAu(seq(200, 1500))
library(ggplot2)

params <- expand.grid(order = c(1, 2, Inf), mode = c("EM", "Magnetic", "Electric"), stringsAsFactors=FALSE)

all <- mdply(params, mie, wavelength=gold$wavelength, epsilon=gold$epsilon, radius=0.08, medium=1.5,
            .progress="text")

m <- melt(all, meas = c("extinction", "scattering", "absorption"))

ggplot(m) +
 facet_grid(mode~variable, scales="free") +
 geom_path(aes(wavelength, value, colour = mode, linetype = factor(order),
               group=interaction(mode, order, variable))) +
 cda::theme_minimal() +
 scale_linetype_manual(values = c(2, 3, 1)) +
 labs(x = expression(wavelength / mu*m),
      y = expression(sigma / mu*m^2),
      colour = "Mode",
      linetype = "Order")
}
\author{
  Baptiste Auguie
}
\seealso{
  Other user: \code{\link{mie-package}},
  \code{\link{mie_bh}}
}

