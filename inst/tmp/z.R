# function stZnAll=GenZnAll(nn_max, rho,sBessel)
# computes the three Zn(rho) auxiliary functions for the radial
# dependence of VSHs for n=1 to nn_max
# can be used for both regular VSHs (based on j(rho)) or
# irregular VSHs (based on h1(rho)).

bessel_z <- function(rho, n_max, type = "h"){
  
  rho <- as.complex(rho)
  
  if(type == "j") {
    f <-  sq * BesselJ(z = rho, nu = 0.5, nSeq = n_max+1) 
    } else if(type == "h") {
    f <-  sq * BesselH(z = rho, nu = 0.5, nSeq = n_max+1, m = 1)
    }
  
  sq <- sqrt((pi/2)/rho)
  # matrix of spherical Bessel
  Zn <- outer(sq, f)
  Z0 <- Zn[,2:(nn_max+1)]
  Z1 <- outer(Z0, 1/rho)
  Z2 <- Zn[,n] - outer(n,Z1)
  
  list(Z0=Z0, Z1=Z1, Z2=Z2)
}




