
library(mie)
require(reshape2)
require(ggplot2)

n_sp = 2.5
n_m = 1.5

scale = 10

a = scale*1
wavelength = scale*seq(2,20,length=1e3)

cross_sections =  mie(wavelength, n_sp^2, radius=a, 
                       medium=n_m, efficiency=TRUE)

m = melt(cross_sections, id=1, measure=c("extinction"))
m$x = n_m*2*pi/m$wavelength * a

ggplot(m, aes(x, value))+
  geom_line() +
  labs(x="x", y=expression(sigma/(pi*a^2)), colour="") 
