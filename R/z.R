function stZnAll=GenZnAll(nNmax, rho,sBessel)
% computes the three Zn(rho) auxiliary functions for the radial
% dependence of VSHs for n=1 to nNmax
% can be used for both regular VSHs (based on j(rho)) or
% irregular VSHs (based on h1(rho)).
%
% Parameters:
  % - nNmax: scalar integer
%          number of n in series 
% - rho:   column vector [R x 1] (no zero components allowed, even for regular VSH, for speed optimization)
%          arguments of the VSHs
% - sBessel: string defining the Bessel function to be used
%           sBessel='j' for or sBessel='h1'
%           


vsh_bessel <- function(rho, nmax, type = "h"){
  
  rho <- as.complex(rho)
  sq <- sqrt((pi/2)*rho)
  
  f <- if(type == "j")  sq * BesselJ(z = rho, nu = 0.5, nSeq = nmax+1) else
                        sq * BesselH(z = rho, nu = 0.5, nSeq = nmax+1, m = 1)
  
  psi <- bj[, -1L]
  xi <- bh[, -1L]
  norho <- outer(1/rho, seq_len(nmax))
  
  dpsi <- bj[, -(nmax+1)] - norho * psi
  dxi <- bh[, -(nmax+1)] - norho * xi
  
  list(psi=psi, xi=xi, dpsi=dpsi, dxi=dxi)
}

n=1:nNmax;
nm1=0:nNmax;
nu=nm1+0.5;

[numat, rhomat] = meshgrid(nu,rho);
if strcmp(sBessel,'h1')
f=besselh(numat, rhomat);
else
  if strcmp(sBessel,'j')
f=besselj(numat,rhomat);
else
  disp 'Error in GetZnAll: wrong sBessel string...'
end;
end;

% f is matrix [R x nNmax+1] of cylindrical Bessel
% Z_{n+0.5}(rho), n=0..nNmax

sq=sqrt((pi/2)./rho); % [R x 1]
f=bsxfun(@times,f,sq); % [R x nNmax+1]
% f is now matrix of spherical Bessel
% z_n(rho), n=0..nNmax or equivalently z_{n-1}(rho), n=1..nNmax+1

stZnAll.Z0=f(:,2:(nNmax+1));
stZnAll.Z1=bsxfun(@times,stZnAll.Z0, 1./ rho);

% Computes: Z2_n=z_{n-1} - n Z1_n
% check for loss of precision in sum
stZnAll.Z2=GenCheckSum2Mat(f(:,n), - bsxfun(@times,n,stZnAll.Z1),'Z2','GenZnAll');



