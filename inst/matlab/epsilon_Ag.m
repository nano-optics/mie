function epsAg = epsilon_Ag(lambda)
%EPSILON_AG Silver dielectric function
%
% Uses the analytical expression given in Eq. (E.1).
% The exp(-i omega t) convention is assumed.
%
% PARAMETERS:
% - lambda: scalar, vector, or matrix
%           wavelengths in NANOMETERS (nm)
%
% RETURNS: epsilon(lambda) as a complex column vector
%
% DEPENDS: none
%
% FAMILY: user_level, epsilon, utility
%

lambda = lambda(:); % ensure column
eps_infty = 4.0;
lambda_p = 282.0;
mu_p = 17000.0;
epsAg = eps_infty *(1-1./(lambda_p.^2 *( (1./lambda).^2 + 1i./(mu_p.*lambda))));
end
