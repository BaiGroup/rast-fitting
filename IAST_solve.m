function [Q_predicted, x, F] = IAST_solve(M, S, isotherm, minlnP, EoS, options, mode, x0)
% IAST solver for N-component systems
% M(i, j): i-th partial pressure for component j
% S{j}(i, 1:2): component j, i-th pressure & loading
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N])
% options: argument for fsolve() or lsqnonlin()
% mode: 1 uses fsolve() and P*y1/P_i^0 - x_i*gamma_i = 0
%       2 uses lsqnonlin() and ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i),
%       gamma_i cannot be negative during optimizations
%       0 uses MultiStart() for mode 2
% x0: initial guess for lnP_i, z_i

[ndata, N] = size(M);

if nargin < 3 || isempty(isotherm) || isempty(minlnP)
    [isotherm, minlnP, maxlnP] = fit_isotherm(S);
end

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 800, 'MaxIter', 100, 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'iter');

if nargin < 5 || isempty(EoS)
    EoS = @(x)[];  % return [] when calling IAST_func()
end

if nargin < 5 || isempty(options)
    options = options_0;
else
    options = optimset(options_0, options);
end

if nargin < 7
    mode = 1;
end

lnP_mixture = log(M);

if nargin < 8
    z = ones(ndata, N-1) * 1/N;
    x0 = IAST_init(lnP_mixture, z, EoS(z));
end

Q_tot = zeros(ndata, 1);
Z = zeros(ndata, N);
F = zeros(ndata, 2*N-1);
x = zeros(ndata, 2*N-1);

for i = 1 : ndata  % mixture partial pressures
    i
    func = @(x)IAST_func(x, isotherm, minlnP, lnP_mixture(i, :), EoS(x(N+1:2*N-1)'), [], mode);
    if mode == 1
        [x1,fval,exitflag,output,jacobian] = fsolve(func, x0(i, :), options);
    else
        lb = [minlnP, zeros(1, N-1)];
        if nargin < 3 || isempty(isotherm) || isempty(minlnP)
            ub = [maxlnP, ones(1, N-1)];
        else
            ub = [Inf*ones(1, N), ones(1, N-1)];
        end
        if mode == 2
            [x1,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(func, x0(i, :), lb, ub, options);
        else
            prob = createOptimProblem('lsqnonlin', 'objective', func, 'x0', x0(i, :), ...
                'lb', lb, 'ub', ub, 'options', options);
            ms = MultiStart('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'UseParallel', 'always');
            [x1, fval, exitflag, output, allmins] = run(ms, prob, 10);
        end
    end

    x(i, :) = x1;
    F(i, :) = func(x1);
    Z(i, 1:N-1) = x1(N+1:2*N-1);
    Z(i, N) = 1 - sum(Z(i, 1:N-1));
    for j = 1 : N
        Q_tot(i) = Q_tot(i) + Z(i, j) / isotherm{j}(x1(j));
    end
    Q_tot(i) = 1 / Q_tot(i);
end

Q_predicted = Q_tot .* Z;

end