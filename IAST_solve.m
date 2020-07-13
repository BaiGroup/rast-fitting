function y = IAST_solve(M, S)
% IAST solver for N-component systems
% M(i, j): i-th partial pressure for component j
% S{j}(i, 1:2): component j, i-th pressure & loading

tol = 6e-6;
method = 'linear';
options_0 = optimset('FinDiffType', 'central', 'MaxFunEvals', 600, 'MaxIter', 400, 'FunValCheck', 'on', 'Display', 'iter');
options = optimset(options_0, 'TolFun', 1e-2, 'TolX', 1e-2);

[M1, M2] = size(M);
minlnP = zeros(M2, 1);

for i = 1 : M2  % components
    lnP = log(S{i}(:, 1));
    N = S{i}(:, 2);
    isotherm_pp = interp1(lnP, N, method, 'pp');
    isotherm{i} = @(x)ppval(isotherm_pp, x);
    [minlnP(i), fval, exitflag, output, jacobian] = fsolve(isotherm{i}, min(lnP), options_0);
end

lnP_mixture = log(M);

N_tot_mixture = zeros(M1, 1);
Z = zeros(M1, M2);
x0 = zeros(2*M2, 1);
F = zeros(M1, 2*M2-1);

for i = 1 : M1  % mixture partial pressures
    i
    x0(M2+1:2*M2) = 1/M2;
    x0(1:M2) = lnP_mixture(i, :)' - log(x0(M2+1:2*M2));
    factor_scale = 10*(M2-1) / (1-x0(M2));
    x0(M2+1:2*M2-1) = x0(M2+1:2*M2-1) * factor_scale;
    D_scale = diag([ones(1, M2), ones(1, M2-1)/factor_scale]);
    func = @(x)AST_func(D_scale*x, isotherm, minlnP, lnP_mixture(i, :), tol);
    
    [x, fval, exitflag, output, jacobian] = fsolve(func, x0(1:2*M2-1), options);
    F(i, :) = fval;
    Z(i, 1:M2-1) = x(M2+1:2*M2-1) / factor_scale;
    Z(i, M2) = 1 - sum(Z(i, 1:M2-1));
    for j = 1 : M2
        N_tot_mixture(i) = N_tot_mixture(i) + Z(i, j) / isotherm{j}(x(j));
    end
    N_tot_mixture(i) = 1 / N_tot_mixture(i);
end

y = N_tot_mixture .* Z;

end