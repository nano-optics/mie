library(mie)
require(reshape2)
require(ggplot2)

gold <- epsAu(seq(300, 1000))

test <- with(gold, mie_bh(wavelength, epsilon, radius = 40,
                     medium = 1.33, efficiency = FALSE))


cross_sections <- with(gold, mie(wavelength, epsilon, radius=40, 
                                 medium=1.33, efficiency=FALSE))
cross_sections2 <- data.frame(wavelength=gold[, 1], scattering = test[, 3], 
                              absorption=test[, 4])

m <- melt(cross_sections, id=1, measure=c("scattering", "absorption"))
m2 <- melt(cross_sections2, id=1, measure=c("scattering", "absorption"))

p <- 
ggplot(m, aes(wavelength, value,  group=variable))+
  geom_line(aes(colour=variable),alpha=0.5, size=1.5) +
  geom_line(aes(colour=variable), size=1.5,  data=m2, linetype=2) +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="") +
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set2")

p
