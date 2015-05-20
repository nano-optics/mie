theta <- seq(0,pi,length=10)
ntheta <- length(theta)
nmax <- 5
pitau <- pitaun(nmax, theta)





lambda = 400:402'
nMedium = 1.45
r0 = 50
theta = linspace(0, pi, 10)'


wavelength = 400:402
nMedium = 1.45
r0 = 50
theta <- seq(0, pi, length=10)

incident_field <- function(wavelength,nMedium,r0,theta){
  nL <- length(wavelength)
  nT <- length(theta)
  k <- 2*pi*nMedium / wavelength
  kr <- r0 * outer(k, cos(theta))
  # exp(ikM z) is [L x T], obtained by matrix product of [L x 1] by [1 x T]
  phasefact <- exp(1i * kr)
  matsin <- matrix(sin(theta), ncol=nT, nrow=nL, byrow = TRUE) 
  matcos <- matrix(cos(theta), ncol=nT, nrow=nL, byrow = TRUE) 
  phasefact * matsin 
  # Results are all [L x T] matrices
  # They result from Eq. H.16 for e_x and H.76 for E_inc
  list(Ecr= phasefact * matsin ,
       Ect = phasefact * matcos,
       Esf = - phasefact)
}

region <- "outside"
# get regular field from known plane wave expression of incident field
E.reg <- incident_field(wavelength,nMedium,r0,theta)

# calculate irregular VSH expansion (coeffs c,d)
E.irr <- PweEgenThetaAllPhi(lambda,epsilon,stAbcdn1.cn1,stAbcdn1.dn1,r0,theta,'h1',stPinTaun);


Ecr2theta=abs(stE.Ecr).^2;
Ect2theta=abs(stE.Ect).^2;
Esf2theta=abs(stE.Esf).^2;

if nNbTheta==1
disp 'PweEFaverages: only one theta, all averages will be set to zero.'
dtheta=0;
else
dtheta=pi/(nNbTheta-1);
end
sintcolnorm=(dtheta/2)*transpose(sin(theta)); % [T x 1]

% computes surface-averages (integrating
% <f>=0.5*int(f(t)*sin(t)dt using a simple sums in the rectangle approximation)
% Sums are computed as matrix [LR x T]* vector [T x 1] multiplication
% All results are [LR x 1]

% average LFIEF
stE.MLocPerpAve=1/2* Ecr2theta * sintcolnorm;
stE.MLocParaAve=1/2* (Ect2theta +Esf2theta) * sintcolnorm;
stE.MLocAve=stE.MLocPerpAve+stE.MLocParaAve;

% average SERS EF (E4 approximation)
stE.F0E4Ave= ...
(3/8* ( ((Ecr2theta+Ect2theta).^2+Esf2theta.^2) * sintcolnorm) + ...
1/4* ( ((Ecr2theta+Ect2theta).*Esf2theta) * sintcolnorm) );
% average SERS EF (E4 approximation - perpendicular component only)
stE.F0E4PerpAve= 3/8* (Ecr2theta.^2) * sintcolnorm;
% average SERS EF (E4 approximation - parallel component only)
stE.F0E4ParaAve= 1/8* ( (3*(Ect2theta.^2+Esf2theta.^2) + 2*Ect2theta.*Esf2theta) * sintcolnorm );
