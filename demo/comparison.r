## ----setup, echo=FALSE,results='hide'------------------------------------
library(mie)
require(tidyr)
require(ggplot2)

## ----comparison, echo=TRUE------------------------------------

gold <- epsAg(seq(300, 800))
radius <- 50
medium <- 1.33

test <- with(gold, mie_bh(wavelength, epsilon, radius = radius,
                     medium = medium, efficiency = FALSE))


cross_sections <- with(gold, mie(wavelength, epsilon, radius=radius, 
                                 medium=medium, efficiency=FALSE, n_max = 4))
cross_sections2 <- data.frame(wavelength=gold[, 1], scattering = test[, 3], 
                              absorption=test[, 4])

m <- pivot_longer(cross_sections, cols=c("scattering", "absorption"))
m2 <- pivot_longer(cross_sections2, cols=c("scattering", "absorption"))

p <- 
ggplot(m, aes(wavelength, value,  group=name))+
  geom_line(aes(colour=name,linetype='ELR'),alpha=0.5, size=1.0) +
  geom_line(aes(colour=name,linetype='BH'), size=1.0,  data=m2) +
  # scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set2") +
  scale_linetype_manual(values=c(2,1)) +
  labs(x="wavelength /nm", y=expression(sigma/nm^2), colour="",
       subtitle="Comparison of Bohren and Huffman's Fortran code\n and Matlab code by Etchegoin and Le Ru (ported to R)") 
  

p


