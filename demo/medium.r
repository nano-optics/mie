## ----setup, echo=FALSE,results='hide'------------------------------------
library(mie)
require(reshape2)
require(ggplot2)

## ----comparison, echo=TRUE------------------------------------

wavelength <- seq(400, 500, length=100)
mat <- data.frame(wavelength = wavelength, epsilon = 2.1^2 +0i)

model <- function(a=20, medium=1.33){
  
  d <- mie(wavelength, mat$epsilon, radius=a, medium=medium, efficiency=FALSE, nmax = 4)
  
  cross_sections <- data.frame(wavelength=wavelength, scattering = d[, 3], 
                               absorption=d[, 4], extinction=d[, 2])
  
  m <- melt(cross_sections, id=1, measure=c("extinction"))
}

params <- expand.grid(medium=seq(1, 2, by=0.1))
library(dplyr)
all <- purrr::pmap_df(params, model, .id = 'id')
params$id <- as.character(1:nrow(params))
d <- left_join(params, all, by='id')
str(d)
p <- 
  ggplot(d, aes(wavelength, value, colour=factor(medium))) +
  geom_line() + 
  scale_color_hue(h=c(60,360)) +
  labs(x="wavelength /nm", y=expression(sigma[ext]/nm^2), colour=expression(n[medium])) 

p

