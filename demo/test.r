library(mie)
## ?mie
 gold <- epsAu(seq(300, 1000))
 cross_sections <- with(gold, mie(wavelength, epsilon, radius=1, medium=1, efficiency=TRUE))
 cross_sections2 <- rbind( 
                         data.frame(wavelength=300e-3, absorption=0, scattering=0, extinction=0),cross_sections,
                         data.frame(wavelength=1000e-3, absorption=0, scattering=0, extinction=0))
 ## cross_sections
 ## matplot(cross_sections$wavelength, cross_sections[, -1], type="l", lty=1,
 ##         xlab=expression(lambda/mu*m), ylab=expression(sigma/mu*m^2))
 ## legend("topright", names(cross_sections)[-1], col=1:3, lty=1)

m <- melt(cross_sections, id=1, measure=c("scattering", "absorption"))
m2 <- melt(cross_sections2, id=1, measure=c("scattering", "absorption"))
p <- 
ggplot(m, aes(wavelength*1e3, value,  group=variable))+
  geom_polygon(aes(fill=variable), data=m2) +
  geom_line(aes(colour=variable), size=1.5) +
  theme_minimal() +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="") +
  theme(legend.position="none")+
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set2")

p

ggsave("mielarge.pdf", p, width=4, height=4)


library(planar)

parametersAu <- list(epsilon=list(1.0^2, gold$epsilon, 1.0^2),
                   lambda=gold$wavelength*1e3, thickness=c(0, 2000, 0),
                   theta=0, polarisation='p')


dAu <- cbind(wavelength=gold$wavelength, data.frame(do.call(recursive.fresnel2, parametersAu))[, c("R", "T")])
dAu$A <- 1 - dAu$R - dAu$T
dAu2 <- rbind( 
                         data.frame(wavelength=300e-3, R=0, T=0, A=0), dAu,
                         data.frame(wavelength=1000e-3, R=0, T=0, A=0))


m3 <- melt(dAu, id="wavelength", meas=c("R", "A"))
m4 <- melt(dAu2, id="wavelength", meas=c("R", "A"))

p2 <- 
ggplot(m3, aes(wavelength*1e3, value,  group=variable))+
  geom_polygon(aes(fill=variable), data=m4) +
  geom_line(aes(colour=variable), size=1.5) +
  theme_minimal() +
  labs(x="wavelength /nm", y="", colour="") +
  theme(legend.position="none")+
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set2")

p2


ggsave("planar.pdf", p2, width=4, height=4)

