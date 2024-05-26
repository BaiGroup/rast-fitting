function lnP0 = inv_adsorption_potential(Psi, ads_pot, options)

if nargin < 3 || isempty(options)
    options = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 600, 'MaxIter', 100, 'TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');
end

func = @(y)(Psi - ads_pot(y));
[lnP0, fval, exitflag, output] = fsolve(func, 7, options);
end
