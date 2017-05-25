function stIncEabn1 = mie_ab_PWE(nNmax)
% Calculates a_{n,1}, and b_{n,1} for n=1:nNmax for the incident PWE field
% using Eqs. H.74 and H.75.
%
% Parameters:
% - nNmax:  integer scalar
% - s:      column vector [L x 1]
%           wavelength-dependent relative refractive index (Eq. H.45)
% - x:      column vector [L x 1]
%           wavelength-dependent x=kM*a (Eq. H.45)
%
% Returns: structure with 2 fields
%          each a vector [1 x nNmax]
% - stIncEabn1.an1: a_{n,1}
% - stIncEabn1.bn1: b_{n,1}
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

n=1:nNmax;
stIncEabn1.bn1 = i.^(n+1) .* sqrt(pi*(2*n+1));
stIncEabn1.an1 = stIncEabn1.bn1;

end
