md <- function(.p, m=2){

		pl.list <- as.polylist(.p)
			if(m == 1) return(pl.list)
			for(n in seq(1, m-1)){
				pl.list[[n+1]] <- deriv(pl.list[[n]])
			}
	
	(pl.list)
}
#  
# md(Pl[[1]], 1)
# 
# md(Pl[[3]], 2)
# md(Pl[[4]], 4)


Nm <- function(l, phi){
t.m <- seq(0, l)
matrix( sqrt((2*l+1) / (4*pi) * factorial(l-t.m) / factorial(l+t.m)), ncol = l+1, nrow = length(phi), byrow=T)
}


Expm <- function(l, phi){
 t.m <- seq(0, l)
 sapply(t.m, function(.m) exp(1i* .m*phi))
}
# l <- 3
# ee <- Expm(l, phi)
	
sphericalHarmonics <- function(l, theta = seq(0, pi, length=10), phi=seq(0, 2*pi, length=10)){
	require(orthopolynom)
	require(lattice)
	require(grid)

	# 0 < theta < pi : polar angle
	# 0 < phi < 2 pi : azimuth angle
	# l <- 3
	 # cannot index from 0

Pl <- as.polylist(legendre.polynomials(l))

Plm <- lapply(seq_along(Pl), function(ind) md(Pl[[ind]], ind))

Plm.theta <- lapply(seq_along(Plm), function(ind) { # l
	ll <- ind - 1 
	sapply(seq_along(Plm[[ind]]), function(ind2) {
		mm <- ind2 - 1
	(-1)^mm *(1-cos(theta)^2)^(mm/2) * as.function(Plm[[ind]][[ind2]])( cos(theta)) 
	})
	})
	


Nlm <- lapply(seq_along(Plm.theta)-1, Nm, phi=phi) 

explm <- lapply(seq_along(Plm.theta)-1, Expm, phi=phi) 


Ylm <- lapply(seq_along(Plm.theta), function(.l)  Nlm[[.l]] * Plm.theta[[.l]] * explm[[.l]] ) 
Ylm
}


# sphericalHarmonics(1)
