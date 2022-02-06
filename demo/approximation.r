## ----setup, echo=FALSE,results='hide'------------------------------------

library(mie)
require(tidyr)
require(dplyr)
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

all <- purrr::pmap_df(params, exact, medium=1.33,
             .id = 'id')
all$model <- 'exact'
all2 <- purrr::pmap_df(params, approximate, medium=1.33,
              .id = 'id')
all2$model <- 'approx'
params$id <- as.character(1:nrow(params))
all <- left_join(params, all, by='id')
all2 <- left_join(params, all2, by='id')
str(all)
str(all2)

  m <- pivot_longer(all2, 
          cols=c("scattering", "absorption"))
str(all2)
str(m)
# m <- transform(m, order=factor(order, levels=c(Inf, 1, 2)))

p<- 
ggplot(m, aes(wavelength, value, colour=factor(order.x), group=interaction(order.x, model)))+
  facet_grid(name~material, scales="free") +
  geom_line(aes(linetype=model)) +
  ## geom_line(aes(colour=variable), data=m2) +
  labs(x="wavelength /nm", y=expression(sigma/(pi*a^2)), colour="", 
       subtitle="Benchmark of approximation formulas\n for a 60nm metal sphere in water") +
  labs(colour="order", linetype="method")+
  scale_fill_brewer(palette="Pastel2")+
  scale_colour_brewer(palette="Set1") +
  scale_linetype_manual(values=c( "dashed", "solid")) 

p
