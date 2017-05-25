function stSuscep = mie_GDAB(nNmax,s,x)
% Calculates the Mie susceptibilities Gamma_n, Delta_n, A_n, and B_n for n=1:nNmax for a single sphere
%
% Parameters:
% - nNmax:  integer scalar
% - s:      column vector [L x 1]
%           wavelength-dependent relative refractive index (Eq. H.45)
% - x:      column vector [L x 1]
%           wavelength-dependent x=kM*a (Eq. H.45)
%
% Returns: structure with 4 fields
%          each a matrix [L x nNmax] with susceptibilities
%          defined in Eqs. H.12, H.13
% stSuscep.Gamma: Gamma_n
% stSuscep.Delta: Delta_n
% stSuscep.A:     A_n
% stSuscep.B:     B_n
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

z=s.*x;

% get psi_n(x), xi_n(x) and derivatives
stRBx=mie_RB(nNmax,x);
stRBz=mie_RB(nNmax,z);

% auxiliary functions
smat=repmat(s,1,nNmax);
PP1 = stRBz.psi .* stRBx.Dpsi;
PP2 = stRBx.psi .* stRBz.Dpsi;
PP3 = stRBz.psi .* stRBx.Dxi;
PP4 = stRBx.xi .* stRBz.Dpsi;

% numerators
NumGam = - PP1 + smat .* PP2;
NumDel =   PP2 - smat .* PP1;

% denominators
DenDelB = - PP4 + smat .* PP3;
DenGamA =   PP3 - smat .* PP4;

% calculates susceptibilities from Eqs. H.46, H.47, H. 51, H.52
% Make the cell to return
stSuscep.Gamma = NumGam ./ DenGamA;
stSuscep.Delta = NumDel ./ DenDelB;
stSuscep.A = i*smat./DenGamA;
stSuscep.B = i*smat./DenDelB;

end
