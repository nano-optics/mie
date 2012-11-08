library(mie)
library(plyr)
library(reshape2)
library(grid)
library(ggplot2)
gold <- epsAu(seq(400, 900))
silver <- epsAg(seq(300, 800))



exact <- function(..., material){
  material <- get(material)
  mie(...,  wavelength=material$wavelength, epsilon=material$epsilon, efficiency=TRUE)

}

approximate <- function(radius, material="gold", medium=1.336, order=Inf, ...){

  material <- get(material)
  
  s <- sqrt(material$epsilon) / medium
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
  if(order == 1)
    return(dip)

  Delta2 <- (1i*x^5/30*(s2-1)) / (s2+3/2+5/14*x^2-5/2646*x^4*(s^4+30*s2-45)-1i*x^5/30*(s2-1))
  quad <- 
  transform(data.frame(wavelength=material$wavelength,
                       extinction = -2/x^2 * 5 * Re(Delta2),
                       scattering = 2/x^2 * 5* Mod(Delta2)^2),
            absorption = extinction - scattering,
            order = 2)
  if(order == 2)
    return(quad)
  Delta3 <- (4i/4725*x^7*(s2-1)) / (s2 +4/3+7*x^2/135*(s2+4) - 7*x^4/10692*(s^4+8*s^2-32) - 4i*x^2/4725*(s2-1))
  oct <- 
  transform(data.frame(wavelength=material$wavelength,
                       extinction = -2/x^2 * 7 * Re(Delta3),
                       scattering = 2/x^2 * 7 * Mod(Delta3)^2),
            absorption = extinction - scattering,
            order = 3)
  
  if(order == 3)
    return(oct)
  
  all <- transform(dip,
                   extinction = extinction + quad$extinction + oct$extinction, 
                   scattering = scattering + quad$scattering + oct$scattering,
                   absorption = absorption + quad$absorption + oct$absorption, order=6)
                
  if(order == Inf)
    return(all)
}


params <- expand.grid(order = c(1, 2, 3, Inf), radius = c(0.05, 0.07), material = c("gold", "silver"), 
                      mode = c( "Electric"), stringsAsFactors=FALSE)

all <- mdply(params, exact, medium=1.33,
             .progress="text")

all2 <- mdply(params, approximate, medium=1.33,
             .progress="text")


m <- melt(list("Exact (Mie)" = all, Approximation=all2), id=c("wavelength", "order", "material", "radius"),
          measure=c("scattering", "absorption"))
m <- transform(m, order=factor(order, levels=c(1, 2, 3, 6, Inf)))
mAu <- subset(m, material == "gold")
mAg <- subset(m, material == "silver")

dummy <- expand.grid(wavelength = 0.5, order=factor(1, levels=c(1, 2, 3, Inf)), L1="Approximate", 
                     variable=factor(c("absorption", "scattering")),radius = c(0.05, 0.07), 
                     stringsAsFactors=FALSE)
dummy$value <- c( 2.5,8, 1.7,   7)

maxi <- ddply(mAg, .(radius, variable), summarize, 
              wavelength=wavelength[which.max(value)], value=round(max(value), 1))
p<- 
ggplot(mAg, aes(wavelength*1e3, value))+
  facet_wrap(~radius+variable, as.table=TRUE, scales="free") +
  geom_line(aes(linetype=L1, colour=order, group=interaction(order, L1))) +
  geom_blank(data=dummy) +
  geom_segment(aes(xend=wavelength*1e3+20, yend=value), data=maxi) +  
  geom_text(aes(x=wavelength*1e3+25, label=value), data=maxi, hjust=0, size=3) +  
  ## geom_line(aes(colour=variable), data=m2) +
  theme_minimal() +
  scale_y_continuous(expand=c(0.005, 0)) +
  scale_x_continuous(expand=c(0, 0)) +
  labs(x="Wavelength [nm]", y="", colour="") +
  labs(colour="", linetype="")+
  scale_colour_manual(values=c( RColorBrewer::brewer.pal(4, "Set1"), "black"),
                      guide="none") +
  scale_linetype_manual(values=c( "dashed", "solid")) +
  ## scale_size_manual(values=c(0, 1)) +
  theme_minimal() +
  theme(line =               element_line(colour = "black", size = 0.2, linetype = 1,
                            lineend = "butt"),
        panel.border=element_rect(fill=NA), legend.position=c(0.8, 0.32),
        panel.margin =       unit(2, "mm"),
        legend.justification = "center",
        strip.text.x =       element_blank(),
        axis.ticks.length =  unit(-0.12, "cm"),
        axis.ticks.margin =  unit(0.25, "cm"),
         axis.text.y =        element_blank(),
        ## panel.grid.major =  element_blank(),
        panel.grid.minor =   element_blank(),
         axis.title.x =       element_text(face="bold", vjust=0))

p

ggsave("fig5tmp.pdf", p, width=4, height=4.5)
