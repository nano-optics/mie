library(mie)

gold <- epsAu(seq(300, 800))
silver <- epsAg(seq(300, 800))
library(ggplot2)


foo <- function(..., material){
  material <- get(material)
  mie(...,  wavelength=material$wavelength, epsilon=material$epsilon)

}
params <- expand.grid(order = c(1, 2, Inf), radius = c(0.03, 0.07), material = c("gold", "silver"), 
                      mode = c("EM", "Magnetic", "Electric"), stringsAsFactors=FALSE)

all <- mdply(params, foo, medium=1.0,
            .progress="text")

m <- melt(all, meas = c("extinction", "scattering", "absorption"))
m <- subset(m, mode %in% c("Electric", "Magnetic") & variable != "extinction")

ggplot(m) +
 facet_wrap(~variable + material+ radius, scales="free", ncol=4, as.table=FALSE) +
 geom_line(aes(wavelength, value, colour = factor(order), linetype=mode)) +
 ## scale_linetype_manual(values = c(2, 3, 1)) +
 labs(x = expression(wavelength / mu*m),
      y = expression(sigma / mu*m^2),
      colour = "Order",
      linetype = "Mode") +
  theme_bw()

## ggsave("mie.pdf")
