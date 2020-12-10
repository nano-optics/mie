## ----setup, echo=FALSE,results='hide'------------------------------------
library(mie)
require(reshape2)
require(ggplot2)
library(dplyr)
library(plyr)

## ----comparison, echo=TRUE------------------------------------

wavelength <- seq(200, 1000)
gold <- epsAu(wavelength)

model <- function(a=30, d=10){
test <- mie_bh(wavelength, epsilon.core=gold$epsilon, epsilon.coating = 1.33^2,
               thickness = 0, radius = a, medium = 1.33, efficiency = FALSE)

km <- 2*pi/wavelength
V <- 4/2*pi*a^3
cross_sections <- data.frame(wavelength=gold[, 1], scattering = test[, 3], 
                              absorption=test[, 4], extinction=test[, 2])

# cross_sections$ratio1 <- (cross_sections$scattering / cross_sections$absorption ) / V / km^3

cross_sections$ratio1 <- (cross_sections$extinction / cross_sections$scattering - 2/3) * V * km^3

cross_sections$ratio2 <- (cross_sections$extinction / cross_sections$absorption - 1) / V / km^3
m <- melt(cross_sections, id=1, 
          measure=c("scattering", "absorption", "extinction",'ratio1','ratio2'))
}
params <- expand.grid(a=seq(10, 100, by=10))

params <- expand.grid(a=seq(5, 30, by=5))

all <- mdply(params, model, .progress='text')


p <- 
ggplot(subset(all, variable!='ratio'), aes(wavelength, value, linetype=variable))+
  facet_wrap(~a,scales='free') +
  geom_line() +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_viridis_c()


perm <- data.frame(wavelength = gold$wavelength, ratio = -pi^2*Im(1.33^3 / (gold$epsilon - 1.33^2)))
p1 <- ggplot(subset(all, variable == 'ratio1' & wavelength < 1000), 
       aes(wavelength, value, colour=factor(a),group=a))+
  geom_line() +
  labs(x="wavelength /nm", y=expression((sigma[scat]/sigma[abs] )%*%lambda^3/R^3), colour="R /nm") +
  scale_fill_brewer(palette="Pastel2") +
  # scale_color_viridis_c() +
  # scale_y_log10()+
  # guides(colour='none')
  theme() +
  geom_line(data=perm, aes(wavelength,ratio), inherit.aes = FALSE)

p1
# ggsave(p1, file='ratio.pdf',width=6,height=6)


p2 <- ggplot(subset(all, variable == 'ratio2' & wavelength < 1000), 
             aes(wavelength, value, colour=factor(a),group=a))+
  geom_line() +
  labs(x="wavelength /nm", y=expression((sigma[ext]/sigma[abs] - 1)%*%lambda^3/R^3), colour="R /nm") +
  scale_fill_brewer(palette="Pastel2") +
  # scale_color_viridis_c() +
  # scale_y_log10()+
  theme()

# g <- egg::ggarrange(p1,p2, ncol=2)
# ggsave(g, file='ratios.pdf',width=12,height=6)

# system('open .')# s <- susceptibility(2, 1, 3)
# s$G
# -0.00814193+0.08986456i -4.023e-06+2.005721e-03i -9.1e-10+3.009717e-05i
# -0.1241427+0.3297443i -0.00030793+0.01754531i -1.964e-07+4.432085e-04i
