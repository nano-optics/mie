function CstMulSuscep=mie_GDAB_ml(nNmax,Cs,Cx)
% calculates the Mie susceptibilities for a multilayer sphere problem
% i.e. Gamma^k_n, Delta^k_n, A^k_n, and B^k_n for n=1:nNmax and k=1..nK.
% See Section H.5.2 for more details.
%
% Parameters:
% - nNmax:  integer scalar
% - Cs:     cell K column vectors [L x 1]
%           wavelength-dependent relative refractive index for each
%           interface s_k = k_{k-1}/k_k (see p. 623)
% - x:      cell of K column vectors [L x 1]
%           wavelength-dependent x_k=k_k*a_k (see p. 623)
%
% Returns: a cell of K structures, each with 4 fields 
%          with susceptibilities
%          defined in Eqs. H.108 and H.117
% CstMulSuscep{k}.Gamma: Gamma^k_n [L x nNmax]
% CstMulSuscep{k}.Delta: Delta^k_n [L x nNmax]
% CstMulSuscep{k}.A: A^k_n [L x nNmax]
% CstMulSuscep{k}.B: B^k_n [L x nNmax]
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

% Extract number of interfaces nK
nK=length(Cx);
nNbLambda=length(Cx{1});

% Calculate the necessary Riccati-Bessel functions
CstRBz=cell(1,nK);
CstRBx=cell(1,nK);
for kk=1:nK
    zk=Cs{kk}.*Cx{kk}; % [L x 1]
    xk=Cx{kk};  % [L x 1]

    % get psi_n(sk*xk), xi_n(sk*xk) and derivatives
    CstRBz{kk}=mie_RB(nNmax,zk);

    % get psi_n(xk), xi_n(xk) and derivatives
    CstRBx{kk}=mie_RB(nNmax,xk);
end

CstMulSuscep=cell(1,nK);
% upward recurrence (see page 623) for Gamma and Delta
% A and B can be calculated at the same time
for kk=1:nK
    if kk==1 % implement initial condition inside loop for convenience
        Gamkm1=zeros(nNbLambda,nNmax); % [L x N]
        Delkm1=Gamkm1; % [L x N]
    else
        % Get Gamma^{k-1} and Delta^{k-1}
        Gamkm1=CstMulSuscep{kk-1}.Gamma; % [L x N]
        Delkm1=CstMulSuscep{kk-1}.Delta; % [L x N]
    end
    
    smat=repmat(Cs{kk},1,nNmax);

    % auxialiary functions for Gamma and A
    PP1 = CstRBz{kk}.psi + Gamkm1.*CstRBz{kk}.xi;
    PP2 = CstRBz{kk}.Dpsi + Gamkm1.*CstRBz{kk}.Dxi;
    % Eq. H.112 
    Dnk= smat .* CstRBx{kk}.xi .* PP2 - PP1 .* CstRBx{kk}.Dxi;
    % Eqs. H.110 and H.111 for Gamma
    CstMulSuscep{kk}.Gamma= (PP1 .* CstRBx{kk}.Dpsi - smat .* CstRBx{kk}.psi .* PP2) ./ Dnk;
    % Eqs. H.118, H.119, H.120 for A
    CstMulSuscep{kk}.A= -1i*smat./ Dnk;
    
    % Same method for Delta and B
    PP1 = CstRBz{kk}.psi + Delkm1.*CstRBz{kk}.xi;
    PP2 = CstRBz{kk}.Dpsi + Delkm1.*CstRBz{kk}.Dxi;
    % Eq. H.115
    Dnk= smat .* CstRBx{kk}.Dxi .* PP1 - PP2 .* CstRBx{kk}.xi;
    % Eqs. H.113 and H.114
    CstMulSuscep{kk}.Delta= (PP2 .* CstRBx{kk}.psi - smat .* CstRBx{kk}.Dpsi .* PP1) ./ Dnk;
    % Eqs. H.121, H.122, H.123 for B
    CstMulSuscep{kk}.B= 1i*smat./ Dnk;

end

end