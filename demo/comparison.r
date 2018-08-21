## ----setup, echo=FALSE,results='hide'------------------------------------
library(mie)
require(reshape2)
require(ggplot2)

## ----comparison, echo=TRUE------------------------------------

gold <- epsAg(seq(200, 600))

test <- with(gold, mie_bh(wavelength, epsilon, radius = 25,
                     medium = 1.0, efficiency = FALSE))


cross_sections <- with(gold, mie(wavelength, epsilon, radius=25, 
                                 medium=1.0, efficiency=FALSE, nmax = 4))
cross_sections2 <- data.frame(wavelength=gold[, 1], scattering = test[, 3], 
                              absorption=test[, 4])

m <- melt(cross_sections, id=1, measure=c("scattering", "absorption"))
m2 <- melt(cross_sections2, id=1, measure=c("scattering", "absorption"))

p <- 
ggplot(m, aes(wavelength, value,  group=variable))+
  geom_line(aes(colour=variable),alpha=0.5, size=1.0) +
  geom_line(aes(colour=variable), size=1.0,  data=m2, linetype=2) +
  labs(x="wavelength /nm", y=expression(sigma), colour="") +
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set2")

p


write.table(file='mie_25nm_Ag.txt', cross_sections, row.names = F)
system("open .")
