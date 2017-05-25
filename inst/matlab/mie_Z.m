function stZnAll=mie_Z(nNmax, rho,sBessel)
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
%
% Returns: stZnAll structure with 3 fields
%          containing matrices [R x nNmax]
% - stZnAll.Z0 is Z_n^0(rho)
% - stZnAll.Z1 is Z_n^1(rho)
% - stZnAll.Z2 is Z_n^2(rho)
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

if any(rho==0)
    disp 'Warning: rho=0 arguments not allowed in mie_Z...'
end

n=1:nNmax;
nm1=0:nNmax;
nu=nm1+0.5;

[numat, rhomat] = meshgrid(nu,rho);

switch sBessel
  case 'h1'
    f=besselh(numat, rhomat);
  case 'j'
    f=besselj(numat, rhomat);
  otherwise
    warning('Error in mie_Z: wrong sBessel string.')
end


% f is matrix [R x nNmax+1] of cylindrical Bessel
% Z_{n+0.5}(rho), n=0..nNmax

sq=sqrt((pi/2)./rho); % [R x 1]
f=bsxfun(@times,f,sq); % [R x nNmax+1]
% f is now matrix of spherical Bessel
% z_n(rho), n=0..nNmax or equivalently z_{n-1}(rho), n=1..nNmax+1

stZnAll.Z0=f(:,2:(nNmax+1));
stZnAll.Z1=bsxfun(@times,stZnAll.Z0, 1./ rho);

% Computes: Z2_n=z_{n-1} - n Z1_n
stZnAll.Z2=f(:,n) - bsxfun(@times,n,stZnAll.Z1);

end
