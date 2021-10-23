function y = adsorption_potential(isotherm, minlnP, maxlnP, tol)

if nargin < 4 || isempty(tol)
    tol = 1e-5;
end

if maxlnP < minlnP + tol
    y = 0;
else
    y = trapz(minlnP:tol:maxlnP, isotherm(minlnP:tol:maxlnP));
end
end
