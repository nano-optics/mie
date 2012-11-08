library(mie)

library(ggplot2)
gold <- epsAu(seq(400, 900))
silver <- epsAg(seq(300, 800))



exact <- function(..., material){
  material <- get(material)
  mie(...,  wavelength=material$wavelength, epsilon=material$epsilon)

}
params <- expand.grid(order = c(1, 2, Inf), radius = c(0.03, 0.07), material = c("gold", "silver"), 
                      mode = c("EM", "Magnetic", "Electric"), stringsAsFactors=FALSE)

all <- mdply(params, foo, medium=1.33,
            .progress="text")


approximate <- function(radius, material="gold", medium=1.336){

  material <- get(material)
  
  s <- sqrt(epsilon) / medium
  k <- 2*pi*medium / material$wavelength
  x <- k*radius
  s2 <- s^2
  D1_0 <- 2i/3 * (s2 - 1)/(s2 + 2) * x^3
  DE1_RC <-  D1_0 / (1 - 3/5*x^2*(s2 - 2)/(s2 + 2) - D1_0 - 3/350*x^4 * (s^4 - 24*s2 + 16) / (s2 + 2))
  ## Delta <- D1_0
  Delta <- DE1_RC
  ## Delta <- (2i/3*x^3*(s2 - 1)) / (s2+2 - 3/5*x^2*(s2-2) - 2i/3*x^3*(s2-1) - 3/350*x^4*(s^4 - 24*s2+16))
  dip <- 
  transform(data.frame(wavelength=material$wavelength,
                       extinction = -2/x^2 * 3* Re(Delta),
                       scattering = 2/x^2 * 3* Mod(Delta)^2),
            absorption = extinction - scattering,
            order = 1)

  Delta2 <- (1i*x^5/30*(s2-1)) / (s2+3/2+5/14*x^2-5/2646*x^4*(s^4+30*s2-45)-1i*x^5/30*(s2-1))
  quad <- 
  transform(data.frame(wavelength=material$wavelength,
                       extinction = -2/x^2 * 3* Re(Delta2),
                       scattering = 2/x^2 * 3* Mod(Delta2)^2),
            absorption = extinction - scattering,
            order = 2)
  all <- transform(dip,
                   extinction = extinction + quad$extinction, 
                   scattering = scattering + quad$scattering,
                   absorption = absorption + quad$absorption, order=Inf)
                   
  rbind(dip, quad, all)
}


rad <- 0.1
test <- approximate(rad, gold$wavelength, gold$epsilon)


cross_sections <- with(gold, mie(wavelength, epsilon, radius=rad, medium=1.33, efficiency=TRUE,
                                 mode = c("Electric"), order=1))

m <- melt(cross_sections, id=c("wavelength", "order"), measure=c("scattering", "absorption"))
m2 <- melt(test, id=c("wavelength", "order"), measure=c("scattering", "absorption"))

p <- 
ggplot(m, aes(wavelength*1e3, value))+
  facet_grid(variable~order) +
  geom_line(aes(colour=variable), size=1.5) +
  geom_line(aes(colour=variable), size=1.5, data=m2, linetype=2) +
  theme_minimal() +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="") +
  theme(legend.position="none")+
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set2")

p
