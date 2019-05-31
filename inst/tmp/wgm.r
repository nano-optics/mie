## ----setup, echo=FALSE,results='hide'------------------------------------
library(mie)
require(reshape2)
require(ggplot2)

## ----comparison, echo=TRUE------------------------------------

wavelength <- seq(300, 1000)
gold <- epsAu(wavelength)

model <- function(a=30, d=10){
test <- mie_bh(wavelength, epsilon.core=gold$epsilon, epsilon.coating = 1.5^2,
               thickness = d, radius = a, medium = 1.33, efficiency = TRUE)

cross_sections <- data.frame(wavelength=gold[, 1], scattering = test[, 3], 
                              absorption=test[, 4], extinction=test[, 2])

m <- melt(cross_sections, id=1, measure=c("scattering", "absorption", "extinction"))
}
params <- expand.grid(a=seq(5, 40, by=5), d=seq(100, 150, by=5))
all <- mdply(params, model, .progress='text')
p <- 
ggplot(all, aes(wavelength, value, linetype=variable))+
  facet_grid(.~a) +
  geom_line(aes(colour=d,group=interaction(variable,d))) +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_viridis_c()

p

# s <- susceptibility(2, 1, 3)
# s$G
# -0.00814193+0.08986456i -4.023e-06+2.005721e-03i -9.1e-10+3.009717e-05i
# -0.1241427+0.3297443i -0.00030793+0.01754531i -1.964e-07+4.432085e-04i
