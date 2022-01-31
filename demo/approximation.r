## ----setup, echo=FALSE,results='hide'------------------------------------

library(mie)
require(reshape2)
require(plyr)
library(ggplot2)
theme_set(theme_bw() + theme(strip.background = element_blank()))
gold <- epsAu(seq(400, 900))
silver <- epsAg(seq(300, 800))

exact <- function(..., material){
  material <- get(material)
  mie(...,  wavelength=material$wavelength, epsilon=material$epsilon, efficiency=TRUE)
}

approximate <- function(..., material){
  material <- get(material)
  mie_approximation(...,  wavelength=material$wavelength, epsilon=material$epsilon, efficiency=TRUE)
}


## ----comparison, echo=TRUE------------------------------------

params <- expand.grid(order = c(1, 2, Inf), radius = c(60), 
                      material = c("gold", "silver"), 
                      mode = c( "EM"), stringsAsFactors=FALSE)

all <- mdply(params, exact, medium=1.33,
             .progress="text")

all2 <- mdply(params, approximate, medium=1.33,
             .progress="text")


m <- melt(list(exact = all, approximate=all2), 
          id=c("wavelength", "order", "material", "radius"),
          measure=c("scattering", "absorption"))
m <- transform(m, order=factor(order, levels=c(Inf, 1, 2)))
mAu <- subset(m, material == "gold")
mAg <- subset(m, material == "silver")

p<- 
ggplot(m, aes(wavelength, value, colour=order, group=interaction(order, L1)))+
  facet_grid(variable~material, scales="free") +
  geom_line(aes(linetype=L1)) +
  ## geom_line(aes(colour=variable), data=m2) +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="", 
       subtitle="Benchmark of approximation formulas\n for a 60nm metal sphere in water") +
  labs(colour="order", linetype="method")+
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set1") +
  scale_linetype_manual(values=c( "dashed", "solid")) 

p
