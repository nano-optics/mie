function stPinTaun=mie_pitau(nNmax, theta)
% Computes angle functions pi_n(cos(theta)) and tau_n(cos(theta)) for n=1..nNmax
%
% Parameters:
% - nNmax: scalar integer [1 x 1]
% - theta: column vector [T x 1]
%          with theta (in radians)
%          all theta's must be between 0 and pi
%
% Returns: structure with fields:
% - pin: matrix [T x nNmax] with pi_n(cos(theta)) n=1..nNmax
% - taun: matrix [T x nNmax] with tau_n(cos(theta)) n=1..nNmax
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

theta = theta(:);
if any(theta < 0)
    disp 'Warning: theta must be >0 in PwePinTaun...';
end

nrows=length(theta);
mu=cos(theta);
% Initialize recurrence pi_0=0, pi_1=1
% pinm1 contains pi_{n-1} n=1..nNmax+1
pinm1=zeros(nrows,nNmax+1);
pinm1(:,2)=ones(nrows,1);

% Get pi_2 to pi_nNmax by recurrence (see SPlaC guide)
% pi_n is pinm1(:,n+1)
for n=2:(nNmax)
    pinm1(:,n+1)=(2*n-1)/(n-1)*mu.*pinm1(:,n)-n/(n-1)*pinm1(:,n-1);
end

% return pi_n matrix (except n=0)
stPinTaun.pin=pinm1(:,2:(nNmax+1));

% return tau_n matrix
nmat=repmat(1:nNmax,nrows,1);
mumat=repmat(mu,1,nNmax);
stPinTaun.taun=nmat.*mumat.*(stPinTaun.pin) - (nmat+1).*pinm1(:,1:nNmax);

end
