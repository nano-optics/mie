## ----setup, echo=FALSE,results='hide'------------------------------------
library(mie)
require(tidyr)
require(ggplot2)

## ----comparison, echo=TRUE------------------------------------

wavelength <- seq(500, 600, length=2e3)
mat <- data.frame(wavelength = wavelength, epsilon = 1.5^2 +0i)

model <- function(a=30, d=0, medium=1.33){
  test <- mie_bh(wavelength, epsilon.core=mat$epsilon, epsilon.coating = 1.5^2,
                 thickness = d, radius = a, medium = medium, efficiency = FALSE)
  
  cross_sections <- data.frame(wavelength=mat[, 1], scattering = test[, 3], 
                               absorption=test[, 4], extinction=test[, 2])
  
  # m <- pivot_longer(cross_sections, id=1, measure=c("scattering", "absorption", "extinction"))
  m <- pivot_longer(cross_sections, cols=c( "extinction"))
}

params <- expand.grid(a=seq(500, 1500, by=500))
# params <- expand.grid(a=seq(5000, 8000, by=500))
all <- purrr::pmap_df(params, model, medium=1.0, .id = 'id')
params$id <- as.character(1:nrow(params))
all <- left_join(params, all, by='id')

p <- 
  ggplot(all, aes(wavelength, value, linetype=name))+
  facet_wrap(~a, scales="free") +
  geom_line(aes(colour=factor(a),group=interaction(name))) +
  labs(x="wavelength /nm", y=expression(sigma[ext]/nm^2), colour="radius /nm") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Set1")

p

# ggsave("glassinair.pdf",width=10,height=3)
# s <- susceptibility(2, 1, 3)
# s$G
# -0.00814193+0.08986456i -4.023e-06+2.005721e-03i -9.1e-10+3.009717e-05i
# -0.1241427+0.3297443i -0.00030793+0.01754531i -1.964e-07+4.432085e-04i
