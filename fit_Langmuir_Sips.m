function [isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_Langmuir_Sips(S, M, options)
% Fit numerical interpolants for N single-component isotherms
%   S{j}(i, 1:2): component j, i-th pressure (P) & loading (Q)
%   method is passed to interp1 and defaults to 'linear'. It is used in
%     interpolating single-component isotherms.
%   options: argument for lsqnonlin()
%
% OUTPUT:
%   isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
%     minlnP(i), 1<=i<=N: lnP at which Q is zero

if nargin < 2 || isempty(M)
    M = 1;
end

if nargin < 3 || isempty(options)
    options = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 200000, 'MaxIter', 10000, 'Display', 'iter');
end

N = length(S);
minlnP = zeros(1, N);
maxlnP = zeros(1, N);
isotherm = cell(1, N);
ads_pot = cell(1, N);
inv_ads_pot = cell(1, N);
for i = 1 : N  % components
    lnP = log(S{i}(:, 1));
    maxlnP(i) = max(lnP);
    minlnP(i) = -Inf;
    x0 = ones(M, 3);
    func = @(x)(relative_error_safe(Langmuir_Sips(lnP, x(:, 1), x(:, 2), x(:, 3)), S{i}(:,2)));
    [param, resnorm, residual, exitflags, output, lambda, jacobian] = lsqnonlin(func, x0, [], [], options);
    param
    isotherm{i} = @(x)Langmuir_Sips(x, param(:, 1), param(:, 2), param(:, 3));
    ads_pot{i} = @(x)Langmuir_Sips_ads_pot(x, param(:, 1), param(:, 2), param(:, 3));
    if M == 1
        inv_ads_pot{i} = @(x)Langmuir_Sips_inv_ads_pot(x, param(1), param(2), param(3));
    else
        inv_ads_pot{i} = [];
    end
end

end