library(mie)


gold <- epsAu(seq(300, 1000))

test <- with(gold, mie_spectrum(wavelength, epsilon, radius = 0.04,
                     medium = 1.333, efficiency = FALSE, order = "all"))


cross_sections <- with(gold, mie(wavelength, epsilon, radius=0.04, medium=1.333, efficiency=FALSE))
cross_sections2 <- data.frame(wavelength=gold[, 1], scattering = test[, 3], absorption=test[, 2])

m <- melt(cross_sections, id=1, measure=c("scattering", "absorption"))
m2 <- melt(cross_sections2, id=1, measure=c("scattering", "absorption"))

p <- 
ggplot(m, aes(wavelength*1e3, value,  group=variable))+
  geom_line(aes(colour=variable), size=1.5) +
  geom_line(aes(colour=variable), size=1.5, data=m2, linetype=2) +
  theme_minimal() +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="") +
  theme(legend.position="none")+
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set2")

p
