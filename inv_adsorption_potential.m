function lnP0 = inv_adsorption_potential(Psi, isotherm, minlnP, tol, options)

if nargin < 4 || isempty(tol)
    tol = 1e-4;
end

if nargin < 5 || isempty(options)
    options = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 600, 'MaxIter', 100, 'TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');
end

func = @(y)(Psi - adsorption_potential(y, isotherm, minlnP, tol));
[lnP0, fval, exitflag, output] = fzero(func, 1e3*minlnP, options);
end
