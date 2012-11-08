Eint <- function(cd, k, r, phi=0, theta,  nmax = length(cd$cn), pt=Mie.pt(cos(theta),nmax)){
	
	kr <- k*r
	j <- besselJ(kr, 0:(nmax+1))
	Es.r <- 0
	Es.theta <- 0
	Es.phi <- 0
	
	for(n in 1:nmax){
		cn <- cd$cn[n]
		dn <- cd$dn[n]
		pn <- pt[n, 1]
		taun <- pt[n, 2]
		En <- (1i)^n* (2*n+1)/(n*(n+1))
		Es.r.n <- En*(-1i*dn*cos(phi)*n*(n+1)*sin(theta)*pn*j[n+1]/kr)
		Es.theta.n <- En*(cn*cos(phi)*pn*j[n+1] - 1i*dn*cos(phi)*taun * ((n*j[n] - (n+1)*j[n+2])/(2*n+1) + j[n+1]/kr))
		Es.phi.n <- En*(-cn*sin(phi)*taun*j[n+1] + 1i*dn*sin(phi)*pn * ((n*j[n] - (n+1)*j[n+2])/(2*n+1) + j[n+1]/kr))
		Es.r <- Es.r + Es.r.n
		Es.theta <- Es.theta + Es.theta.n
		Es.phi <- Es.phi + Es.phi.n
	}
	
	Esca <- c(Es.r=Es.r, Es.theta=Es.theta, Es.phi=Es.phi)
	Esca
}

Escat <- function(ab, k, r, phi=0, theta,  nmax = length(ab$an), pt=Mie.pt(cos(theta),nmax)){
	
	kr <- k*r
	hn <- hankel(0:(nmax+1), kr)

	Es.r <- 0
	Es.theta <- 0
	Es.phi <- 0
	
	for(n in 1:nmax){
		an <- ab$an[n]
		bn <- ab$bn[n]
		pn <- pt[n, 1]
		taun <- pt[n, 2]
		En <- (1i)^n* (2*n+1)/(n*(n+1))
		Es.r.n <- En*(1i*an*cos(phi)*n*(n+1)*sin(theta)*pn*hn[n+1]/kr)
		Es.theta.n <- En*(-bn*cos(phi)*pn*hn[n+1] + 1i*an*cos(phi)*taun * ((n*hn[n] - (n+1)*hn[n+2])/(2*n+1) + hn[n+1]/kr))
		Es.phi.n <- En*(bn*sin(phi)*taun*hn[n+1] - 1i*an*sin(phi)*pn * ((n*hn[n] - (n+1)*hn[n+2])/(2*n+1) + hn[n+1]/kr))
		Es.r <- Es.r + Es.r.n
		Es.theta <- Es.theta + Es.theta.n
		Es.phi <- Es.phi + Es.phi.n
	}
	
	Esca <- c(Es.r=Es.r, Es.theta=Es.theta, Es.phi=Es.phi)
	Esca
}

hankel <- function(n, x){
	sqx  <- sqrt(0.5*pi/x)
	nu <- (n+0.5)
	bx   <- besselJ(x,nu)*sqx
	yx   <- besselY(x,nu)*sqx
	
	hx  <-  bx+1i*yx
	hx
}

`Mie.ab` <-
function(m = 1,x = 1){

	z <- m * x
	nmax <- round(Re(2+x+4*x^(1/3)))
	nmx <- round(max(nmax,Mod(z)) + 16)
	n <- (1:nmax)
	nu <- (n + 0.5)

	sx <- sqrt(0.5 * pi * x)
	px <- sx * besselJ(x,nu)
	p1x <- c(sin(x), px[1:(nmax-1)])

	chx <- -sx * besselY(x,nu)
	ch1x <- c(cos(x), chx[1:(nmax-1)])

	gsx <- px-1i*chx
	gs1x <- p1x-1i*ch1x

	dnx <- rep(0i,nmx)
	
	# Computation of Dn(z) according to (4.89) of B+H (1983)
	for (j in seq(from=nmx,to=2,by=-1) ){  
	    dnx[j-1] <- j/z-1/(dnx[j]+j/z)
	}

	dn <- dnx[n]       # Dn(z), n=1 to nmax
	da <- dn/m + n/x
	db <- m * dn + n/x

	an <- (da * px - p1x) / (da * gsx - gs1x)
	bn <- (db * px - p1x) / (db * gsx - gs1x)

	result <- cbind(an, bn)
	return(as.data.frame(result))
}


`Mie.cd` <-
function(m = 1,x = 1){

	z <- m * x
	nmax <- round(Re(2+x+4*x^(1/3)))
	nmx <- round(max(nmax,Mod(z)) + 16)
	n <- (1:nmax)
	nu <- (n + 0.5)
	msq <- m^2

	cnx <- rep(0+0i, nmx)
	for (j in seq(nmx, 2, by=-1)){
	    cnx[j-1] <- j-z*z/(cnx[j]+j)
	}
	cnn <- cnx[n]
	
	sqx <- sqrt(0.5*pi/x) 
	sqz <- sqrt(0.5*pi/z)
	bx  <- besselJ(x, nu)*sqx
	bz  <- 1/besselJ(x, nu)/sqz
	yx <-  besselY(x, nu)*sqx
	hx <-  bx + 1i*yx

	b1x <- c(sin(x)/x, bx[1:(nmax-1)])
	y1x <- c(-cos(x)/x, yx[1:(nmax-1)])
	h1x <-  b1x+1i*y1x
	###
	
	ax  <-  x*b1x-n*bx
	ahx <-  x*h1x-n*hx

	nenn1 <- msq*ahx-hx*cnn

	nenn2 <- ahx-hx*cnn

	cn <- bz*(bx*ahx-hx*ax)/nenn2

	dn <- bz*m*(bx*ahx-hx*ax)/nenn1

	result <- cbind(cn, dn)
	return(as.data.frame(result))
}

`Mie.coated.ab` <-
function(m1 = 1, m2 =1,  x = 1, y=1){

m <- m2 / m1

# z <- m * x
# nmax <- round(Re(2+x+4*x^(1/3)))
# nmx <- round(max(nmax,Mod(z)) + 16)
# n <- (1:nmax)
# nu <- (n + 0.5)


u <- m1 * x
v <- m2 * x 
w <- m2 * y      # The arguments of Bessel Functions

nmax <- round(Re(2+y+4*y^(1/3)))     # The various nmax values
mx <- max(Mod(m1*y),Mod(m2*y))
nmx <- round(max(nmax,mx) + 16)
nmax1 <- nmax-1
n <- (1:nmax)

# Computation of Dn(z), z=u,v,w according to (4.89) of B+H (1983)

z <- u
dnx <- rep(0i,nmx)
for (j in seq(from=nmx,to=2,by=-1) ){  
    dnx[j-1] <- j/z-1/(dnx[j]+j/z)
}
dnu <- dnx[n]

z <- v
for (j in seq(from=nmx,to=2,by=-1) ){  
    dnx[j-1] <- j/z-1/(dnx[j]+j/z)
}
dnv <- dnx[n]

z <- w
for (j in seq(from=nmx,to=2,by=-1) ){  
    dnx[j-1] <- j/z-1/(dnx[j]+j/z)
}
dnw <- dnx[n]

# Computation of Psi, Chi and Gsi Functions and their derivatives

nu <- (n+0.5)

sv <- sqrt(0.5*pi*v)
pv <- sv * besselJ(v, nu)

sw <- sqrt(0.5*pi*w)
pw <- sw*besselJ(w,  nu)

sy <- sqrt(0.5*pi*y)
py <- sy*besselJ(y, nu) 

p1y <- c(sin(y), py[1:(nmax-1)])
chv <- -sv * besselY(v,nu)
chw <- -sw * besselY(w,nu)
chy <- -sy * besselY(y,nu)

ch1y <- c(cos(y), chy[1:(nmax-1)])

gsy <- py-1i*chy
gs1y <- p1y-1i*ch1y

# Computation of U, V, F Functions, avoiding products of Riccati-Bessel Fcts.

uu <- m*dnu-dnv
vv <- dnu/m-dnv

fv <- pv/chv
fw <- pw/chw

ku1 <- uu*fv/pw
kv1 <- vv*fv/pw

ku2 <- uu*(pw-chw*fv) + (pw/pv)/chv
kv2 <- vv*(pw-chw*fv)+(pw/pv)/chv

dns1 <- ku1/ku2
gns1 <- kv1/kv2

# Computation of Dn_Schlange, Gn_Schlange

dns <- dns1+dnw
gns <- gns1+dnw

a1 <- dns/m2+n/y
b1 <- m2*gns+n/y

# an and bn

an <- (py*a1-p1y)/(gsy*a1-gs1y)

bn <- (py*b1-p1y)/(gsy*b1-gs1y)
result <- cbind(an, bn)
return(as.data.frame(result))
}

`Mie.coated` <-
function(m1, m2, x, y){

if (x==y)     { # To reduce computing time
	result <- Mie(m1,y)
	return(result)
}                

if (x==0)     {      # To avoid a singularity at x=0
	result <- Mie(m2,y)
	return(result)
}

if (m1==m2)     { # To reduce computing time
	result <- Mie(m1,y)
	return(result)
}

# This is now the normal situation

	nmax <- round(2 + y + 4 * y^(1/3))
	n1 <- nmax-1

	n <- (1:nmax)
	cn <- 2*n+1
	c1n <- n * (n+2) / (n+1)
	c2n <- cn / n / (n+1)

	y2 <- y^2

	ab <- Mie.coated.ab(m1, m2, x, y)

	anp <- Re(ab$an)
	anpp <- Im(ab$an)

	bnp <- Re(ab$bn)
	bnpp <- Im(ab$bn)
	g1 <- matrix(0,nrow=nmax,ncol=4) 

	g1[1:n1,] <- cbind(anp[2:nmax], anpp[2:nmax], bnp[2:nmax], bnpp[2:nmax])

	dn <- cn * (anp+bnp)
	q <- sum(dn)

	qext <- 2 * q / y2

	en <- cn * (anp^2 + anpp^2 + bnp^2 + bnpp^2)

	q <- sum(en)
	qsca <- 2 * q / y2
	qabs <- qext - qsca

	result <- c(qext, qabs, qsca)
	return(result)
}

`Mie.coated.spectra` <-
function(lambda, core.n, inner.radius = 0.04 , 
	layer.thickness = 0.005, layer.n = 1.5, 
	medium = 1.333){
	
	if (length(layer.n) == 1 & length(core.n) == 1) {
		stop("both layer.n and core.epsilon are scalars, run Mie.coated instead")
	}
	if (length(layer.n) == 1 ) {
		layer.n <- rep(layer.n, length = length(core.n)) # case a scalar is given (non-dispersive)
	}
	if (length(core.n) == 1 ) {
		core.n <- rep(core.n, length = length(layer.n)) # case a scalar is given (non-dispersive)
	}
	
	outer.radius <- inner.radius + layer.thickness
	area <- pi * outer.radius^2
	k0 <- 2 * pi / lambda
	k <- k0 * medium
	
	x <- k * inner.radius
	y <- k * outer.radius
	
	m1 <-  core.n / medium
	m2 <-  layer.n / medium
	
	efficiencies <- mapply(Mie.coated,  m1=m1, m2 = m2, x = x, y = y)
	
	return(area*t(efficiencies)) # cross sections, this is a matrix that can be used with mapply to sweep over params
}

Mie.pt <- function (u,nmax){
#
# pi_n and tau_n, -1 <= u <= 1, n1 integer from 1 to nmax 
# angular functions used in Mie Theory
# Bohren and Huffman (1983), p. 94 - 95

p <- rep(1, 3*nmax)
tau <- rep(u, 3*nmax)
p[2] <- 3*u
tau[2] <- 3*cos(2*acos(u))

for (n1 in (3:nmax)){
	
    p1 <- (2*n1-1)/(n1-1)*p[n1-1]*u
    p2 <- n1/(n1-1)*p[n1-2]
    p[n1] <- p1-p2

    t1 <- n1*u*p[n1]
    t2 <- (n1+1)*p[n1-1]
    tau[n1] <- t1-t2
	
}

as.data.frame(cbind(p, tau))
}

Mie.S <- function(m, x, u){


# Computation of Mie Scattering functions S1 and S2
# for complex refractive index m=m'+im", 
# size parameter x=k0*a, and u=cos(scattering angle),
# where k0=vacuum wave number, a=sphere radius;
# s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
# C. Mätzler, May 2002

nmax <- round(2+x+4*x^(1/3))

ab <- Mie.ab(m,x)

an <- ab$an
bn <- ab$bn

pt <- Mie.pt(u,nmax)

pin <- pt$p
tin <- pt$tau

n <- 1:nmax
n2 <- (2*n+1)/(n*(n+1))

pin <- n2*pin
tin <- n2*tin

S1 <- (an*pin+bn*tin)
S2 <- (an*tin+bn*pin)

as.data.frame(cbind(S1, S2))
}

##' mie_spectrum
##'
##' spectrum for a sphere
##' @title mie_spectrum
##' @param m 
##' @param x 
##' @param u 
##' @return data.frame
##' @export
##' @author Baptiste Auguie (from C. Mätzler Matlab code)
`mie_spectrum` <-
  function( lambda, epsilon, radius = 0.04 , medium = 1.333, efficiency = FALSE , order='all'){
	
	area <- pi * (radius)^2
	k0 <- 2 * pi / lambda
	k <- k0 * medium
	x <- k * radius
	m <- sqrt(epsilon) / medium
        
	efficiencies <- mapply(Mie, m = m, x = x, order=order)
	
		if(efficiency == TRUE)
		return(t(efficiencies)) # um^2
		else
		return(area*t(efficiencies)) # cross sections
}


`Mie` <-
function(m, x, order='all'){

if (x == 0){                # To avoid a singularity at x=0
    result <- list(qext = 0, qsca = 0, qabs = 0)
	return(as.data.frame(result))
}
          # This is the normal situation

    nmax <- round(2 + x + 4 * x^(1/3))
    n1 <- nmax-1

    n   <- (1:nmax)
    cn  <- 2*n+1
    c1n <- n * (n+2) / (n+1)
    c2n <- cn / n / (n+1)

    x2 <- x^2

    ab <- Mie.ab(m,x)

    anp <- Re(ab$an) 
    anpp <- Im(ab$an)

    bnp <- Re(ab$bn) 
    bnpp <- Im(ab$bn)

    g1 <- matrix(0,nrow=nmax,ncol=4) #[0; 0; 0; 0]; # displaced numbers used for

	g1[1:n1,] <- cbind(anp[2:nmax], anpp[2:nmax], bnp[2:nmax], bnpp[2:nmax])

	# dn <- cn * (anp+bnp)
	dn <- switch(order, 
		'all' = cn * (anp+bnp), 
		'dipole' = cn * (anp), 
		'quadrupole' = cn * (anp)
		)

    # q <- sum(dn)
	qq <- switch(order, 
		'all' = sum(dn), 
		'dipole' = dn[1], 
		'quadrupole' = dn[2]
		)

    qext <- 2 * qq / x2 

    # en <- cn * (anp^2 + anpp^2 + bnp^2 + bnpp^2)
	en <- switch(order, 
		'all' = cn * (anp^2 + anpp^2 + bnp^2 + bnpp^2), 
		'dipole' = cn * (anp^2 + anpp^2 ), 
		'quadrupole' = cn * (anp^2 + anpp^2 )
		)

    # qq <- sum(en)
	qq <- switch(order, 
		'all' = sum(en), 
		'dipole' = en[1], 
		'quadrupole' = en[2]
		)

    qsca <- 2 * qq / x2

    qabs <- qext - qsca

result <- c(qext, qabs, qsca)
return(result)
}

