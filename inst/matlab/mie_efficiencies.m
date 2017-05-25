function stQ=mie_efficiencies(x,stGD)
% Calculates the Mie extinction, scattering and absorption coefficients for PWE from the Mie susceptibilities.
% The relevant formulae are given in  Eqs. H.76, H.77, H.78.
%
% Parameters:
% - x:      column vector [L x 1]
%           wavelength-dependent x=kM*a (Eq. H.45)
% - stGD:  structure with two fields, Gamma and Delta
%           each a matrix [L x nNmax]
%           with suceptibilities Gamma_n and Delta_n
%           can be obtained from function GenSuscepGDAB
%
% Returns: structure stQ with three fields
%          each a matrix [L x 1]
% - stQ.Qext: [L x 1] wavelength-dependent extinction coefficient
% - stQ.Qsca: [L x 1] wavelength-dependent scattering coefficient
% - stQ.Qabs: [L x 1] wavelength-dependent absorption coefficient
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

nNmax=size(stGD.Gamma,2);

n=transpose(1:nNmax); % [nNmax x 1]
cc = 2.*n+1; % [nNmax x 1]

% From Eq. H.76
Gamma2=abs(stGD.Gamma).^2;
Delta2=abs(stGD.Delta).^2;
scamat= Gamma2 + Delta2; % [L x nNmax]
stQ.Qsca= 2./((x).^2) .* (scamat * cc); % [L x 1]
% no checksum required since this is a sum of positive numbers

% Calculation of Qext
% This may result in loss-of-precission warnings for non-absorbing sphere,
% for which real(\Delta)=-|\Delta|^2 exactly
% These are not an issue since in this case Qext=Qsca, and only
% the calculation of Qsca should therefore be trusted.

% take real part
GammaR=real(stGD.Gamma); % [L x nNmax]
DeltaR=real(stGD.Delta); % [L x nNmax]
% From Eq. H.77
stQ.Qext= - 2./((x).^2) .* ((GammaR+DeltaR) * cc); % [L x 1]
% sum is on negative numbers only, no loss of precision

% From Eq. H.78
stQ.Qabs = stQ.Qext - stQ.Qsca; % [L x 1]

end