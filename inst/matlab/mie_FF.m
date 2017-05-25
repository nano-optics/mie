function stMieRes=mie_FF(lambda,Cepsilon,Ca,nNmax)
% Solves the problem of PWE of a multilayer spherical system with Mie theory
% The parameters are given as arguments or may also be given in a
% structure stMP when called as mie_FF(stMP,sPlot,sCoeff)
% where stMP contains the parameter fields.
% Note that for PWE, we have for the Mie coefficients:
% |m|=1 only; a,c even in m; b,d odd in m.
% For example, the scattered field is entirely defined by
% c_{n,1} and d_{n,1}
% (c_{n,-1}=c_{n,1}, d_{n,-1}=-d_{n,1}, and
% all other coefficients are zero).
%
% Parameters:
% - nNmax:     integer [1 x 1]
%              number of n in series
% - Ca:        cell of integer {K x [1 x 1]}
%              K is number of interfaces (K=1 for single sphere)
%              radii of spherical interfaces (typically in NANOMETERS)
%              Ca{K} is outer interface (largest sphere)
% - lambda:    possibly column vector [L x 1]
%              wavelengths (typically in NANOMETERS)
% - Cepsilon:  cell of scalars or column vectors {K+1 x [L x 1]}
%              epsilon of media (possibly
%              wavelength-dependent) for k=0 (inside sphere)
%              to k=K (embedding medium).
%
% Returns: structure stMieRes with fields containing the parameters for
%          future reference and additional fields for the results
% - stMieRes.nNmax: scalar
% - stMieRes.a: scalar
% - stMieRes.lambda: [L x 1]
% - stMieRes.epsilonM: scalar or [L x 1]
% - stMieRes.x: [L x 1]
% - stMieRes.nK: scalar, number of interfaces
% - stMieRes.Ca: cell of {K scalars}
% - stMieRes.CepsilonM: cell of {K+1 x (scalar or [L x 1])}
% - stMieRes.Cx: cell of {K x [L x 1]}
% - stMieRes.Cs: cell of {K x [L x 1]}
% - stMieRes.Qext: [L x 1] wavelength-dependent extinction coefficient
% - stMieRes.Qsca: [L x 1] wavelength-dependent scattering coefficient
% - stMieRes.Qabs: [L x 1] wavelength-dependent absorption coefficient
% - stMieRes.MLocAve: [L x 1] wavelength-dependent average LFIEF
% - stMieRes.Gamma (if sCoeff='suscep' or 'both'): [L x nNmax]
% - stMieRes.Delta (if sCoeff='suscep' or 'both'): [L x nNmax]
% - stMieRes.CstMulGDAB (if sCoeff='suscep' or 'both'): [L x nNmax]
% - stMieRes.Cgammakn1 (if sCoeff='coeff' or 'both'): cell {K+1 x [L x nNmax]}
% - stMieRes.Cdeltakn1 (if sCoeff='coeff' or 'both'): cell {K+1 x [L x nNmax]}
% - stMieRes.Calphakn1 (if sCoeff='coeff' or 'both'): cell {K+1 x [L x nNmax]}
% - stMieRes.Cbetakn1 (if sCoeff='coeff' or 'both'): cell {K+1 x [L x nNmax]}
% - stMieRes.cn1 (if sCoeff='coeff' or 'both'): [L x nNmax]
% - stMieRes.dn1 (if sCoeff='coeff' or 'both'): [L x nNmax]
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information


nK=length(Ca); % number of spherical interfaces K

% Get s_k and x_k for k=1..nK
[Cs,Cx]=mie_CsCx(lambda,Ca,Cepsilon);

% Calculate susceptibilities Gamma^k_n, Delta^k_n, A^k_n, and B^k_n
% (defined in Eqs. H.108 and H.117)
% for all n,k and all lambda
CstMulGDAB=mie_GDAB_ml(nNmax,Cs,Cx);
% CstMulGDAB is a cell of K structures, each with the corresponding
% susceptibilities
% Gamma^k_n, Delta^k_n, A^k_n, and B^_n as matrices [L x Nmax]

% Calculate Q coeffs 
stMieRes=mie_efficiencies(Cx{nK},CstMulGDAB{nK});

% return results
stMieRes.nNmax=nNmax;
stMieRes.a=Ca{nK}; % for compatibility with single spheres
stMieRes.lambda=lambda;
stMieRes.epsilonM=Cepsilon{nK+1}; % for compatibility with single spheres
stMieRes.x=Cx{nK}; % for compatibility with single spheres

stMieRes.nK=nK;
stMieRes.Ca=Ca;
stMieRes.Cepsilon=Cepsilon;
stMieRes.Cx=Cx;
stMieRes.Cs=Cs;

stMieRes.Gamma=CstMulGDAB{nK}.Gamma; % for compatibility with single spheres
stMieRes.Delta=CstMulGDAB{nK}.Delta; % for compatibility with single spheres
stMieRes.CstMulGDAB=CstMulGDAB;

% calculate incident PW coefficient an1 and bn1 [1 x nNmax]
stIncEabn1=mie_ab_PWE(nNmax);
an1mat=repmat(stIncEabn1.an1,length(lambda),1); % [L x nNmax]
bn1mat=repmat(stIncEabn1.bn1,length(lambda),1); % [L x nNmax]
% calculate Mie coefficients using a downward recurrence (pp. 623,624)
% each cell of coeff corresponds to a given kk=0..nK
stMieRes.Cgammakn1=cell(1,nK+1);
stMieRes.Cdeltakn1=cell(1,nK+1);
stMieRes.Calphakn1=cell(1,nK+1);
stMieRes.Cbetakn1=cell(1,nK+1);
stMieRes.Cgammakn1{1}=0; % no scattered field inside the smallest sphere
stMieRes.Cdeltakn1{1}=0;
% Initialize recurrence
stMieRes.Calphakn1{nK+1}=an1mat;
stMieRes.Cbetakn1{nK+1}=bn1mat;
for kk=nK:-1:1
    % Eq. H.108 gamma^k=Gamma^k * alpha^k
    stMieRes.Cgammakn1{kk+1}=CstMulGDAB{kk}.Gamma .* stMieRes.Calphakn1{kk+1};
    stMieRes.Cdeltakn1{kk+1}=CstMulGDAB{kk}.Delta .* stMieRes.Cbetakn1{kk+1};
    % Eq. H.117 alpha^{k-1}=A^k alpha^k
    stMieRes.Calphakn1{kk}=CstMulGDAB{kk}.A .* stMieRes.Calphakn1{kk+1};
    stMieRes.Cbetakn1{kk}=CstMulGDAB{kk}.B .* stMieRes.Cbetakn1{kk+1};
end
stMieRes.cn1=stMieRes.Cgammakn1{nK+1}; % for compatibility with single spheres
stMieRes.dn1=stMieRes.Cdeltakn1{nK+1}; % for compatibility with single spheres

    
end
