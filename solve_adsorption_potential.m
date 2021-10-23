function lnP0 = solve_adsorption_potential(Omega, isotherm, minlnP, tol, options)

if nargin < 4 || isempty(tol)
    tol = 1e-5;
end

if nargin < 5 || isempty(options)
    options = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 600, 'MaxIter', 100, 'Display', 'off');
end

lnP0 = minlnP;
N = length(minlnP);
for i = 1:N
    func = @(y)(Omega - adsorption_potential(isotherm{i}, minlnP(i), y, tol));
    [lnP0(i), fval, exitflag, output, jacobian] = fsolve(func, minlnP(i), options);
end
end
