function [Cs,Cx]=mie_CsCx(lambda,Ca,Cepsilon)
% Low-level function to calculate s_k and x_k for a sphere multilayer
% See p. 623 for definitions.
%
% Parameters:
% - lambda:    column vector [L x 1]
%              wavelength in nm
% - Ca:        cell of K scalars
%              a_k: radii of spherical interfaces (in nm) for k=1..K
% - Cepsilon:  cell of K+1 scalars or [L x 1] vectors
%              epsilon of media (possibly
%              wavelength-dependent) for k=0 (inside sphere)
%              to k=K (embedding medium).
% Returns:
% - Cs:        cell of K+1 vectors [L x 1]
%              s_k for k=1..N
% - Cx:        cell of K+1 vectors [L x 1]
%              x_k for k=1..N
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

nK=length(Ca); % number of spherical interfaces K
% Calculate x_k=k_k a_k and s_k=k_{k-1}/k_k, see p. 623
Cx=cell(1,nK);
Cs=cell(1,nK);
for kk=1:nK
    % Calculate x_k [L x 1]. Note that Cepsilon{jj} is epsilon in region jj+1
    Cx{kk}=2*pi* sqrt(Cepsilon{kk+1}) * Ca{kk} ./ lambda; % [L x 1]
    % Calculate s_k [L x 1]
    % 0*lambda is used to ensure that s_k is [L x 1] even when
    % both epsilon's are scalar
    Cs{kk}=sqrt(Cepsilon{kk}+0*lambda)./sqrt(Cepsilon{kk+1}); % [L x 1]
end
