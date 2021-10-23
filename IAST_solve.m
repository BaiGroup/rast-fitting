function [Q_predicted, x, err, lnP0, psi, exitflags] = IAST_solve(M, S, isotherm, minlnP, EoS, options, mode, x0, tol, EoS_deriv)
% IAST solver for N-component systems
% M(i, j): i-th partial pressure for component j
% S{j}(i, 1:2): component j, i-th pressure & loading
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N-1, Omega])
% options: argument for fsolve() or lsqnonlin()
% mode: 1 or 3 uses fsolve() and P*y1/P_i^0 - x_i*gamma_i = 0
%       2 or 4 uses lsqnonlin() and ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i),
%       gamma_i cannot be negative during optimizations
%       -2 uses MultiStart() for mode 2
% x0: initial guess for z_i, lnP_i

[ndata, N] = size(M);

if nargin < 4 || isempty(isotherm) || isempty(minlnP)
    [isotherm, minlnP, maxlnP] = fit_isotherm(S);
end

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 800, 'MaxIter', 100, 'TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');

if nargin < 5 || isempty(EoS)
    EoS = @(x)ones(1, N);
end

if nargin < 5 || isempty(options)
    options = options_0;
else
    options = optimset(options_0, options);
end

if nargin < 7 || isempty(mode)
    mode = 1;
end

lnP_mixture = log(M);

if nargin < 8 || isempty(x0)
    z = ones(ndata, N-1) * 1/N;
    x0 = IAST_init(lnP_mixture, z, EoS([z,1e3*ones(ndata,1)]), mode);
end

if nargin < 9 || isempty(tol)
    tol = 1e-5;
end

if nargin < 10 || isempty(EoS_deriv)
    % central difference to calculate numerical derivative
    EoS_deriv = @(x)(sum((log(EoS([x(1:end-1), x(end)+tol]))-log(EoS([x(1:end-1), x(end)-tol])))/2/tol.*x(1:end-1)));
end

Q_tot = zeros(ndata, 1);
Z = zeros(ndata, N);
lnP0 = zeros(ndata, N);
psi = zeros(ndata, N);
exitflags = zeros(ndata, 1);
if mode == 1 || mode == 2 || mode == -2
    err = zeros(ndata, 2*N-1);
    x = zeros(ndata, 2*N-1);
elseif mode == 3 || mode == 4
    err = zeros(ndata, N);
    x = zeros(ndata, N);
else
    error('IAST_solve:UnknownMode','Mode parameter outside known choices.')
end

for i = 1 : ndata  % mixture partial pressures
    if ~strcmp(optimget(options, 'Display'), 'off')
        fprintf('\n======Data point %d ======\n', i);
    end
    func = @(x)IAST_func(x, isotherm, minlnP, lnP_mixture(i, :), EoS, [], mode);
    if mode == 1 || mode == 3
        [x1,fval,exitflags(i),output,jacobian] = fsolve(func, x0(i, :), options);
    elseif mode == 2 || mode == 4 || mode == -2
        if mode == 4
            lb = [zeros(1, N-1), 0];
            ub = [ones(1, N-1), Inf];
        else
            lb = [zeros(1, N-1), minlnP];
            if nargin < 4 || isempty(isotherm) || isempty(minlnP)
                ub = [ones(1, N-1), maxlnP];
            else
                ub = [ones(1, N-1), Inf*ones(1, N)];
            end
        end
        if mode == 2 || mode == 4
            [x1,resnorm,residual,exitflags(i),output,lambda,jacobian] = lsqnonlin(func, x0(i, :), lb, ub, options);
        elseif mode == -2
            prob = createOptimProblem('lsqnonlin', 'objective', func, 'x0', x0(i, :), ...
                'lb', lb, 'ub', ub, 'options', options);
            ms = MultiStart('Display', 'iter', 'TolFun', 1e-3, 'TolX', 1e-4, 'UseParallel', 'always');
            [x1, fval, exitflags(i), output, allmins] = run(ms, prob, 10);
        end
    end

    x(i, :) = x1;
    [err(i, :), lnP0(i, :), psi(i, :)] = func(x1);
    if strcmp(optimget(options, 'Display'), 'final')
        disp('err(i, :): ')
        disp(err(i, :));
    end
    Z(i, 1:N-1) = x1(1:N-1);
    Z(i, N) = 1 - sum(Z(i, 1:N-1));
    for j = 1 : N
        Q_tot(i) = Q_tot(i) + Z(i, j) / isotherm{j}(lnP0(i, j));
    end
    Q_tot(i) = 1 / (Q_tot(i) + EoS_deriv([Z(i, 1:N-1), psi(i, N)]));
end

Q_predicted = Q_tot .* Z;

end