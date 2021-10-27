function [isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_piecewise_polynomial(S, method, options)
% Fit numerical interpolants for N single-component isotherms
%   S{j}(i, 1:2): component j, i-th pressure (P) & loading (Q)
%   method is passed to interp1 and defaults to 'linear'. It is used in
%     interpolating single-component isotherms.
%   options: argument for fsolve()
%
% OUTPUT:
%   isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
%     minlnP(i), 1<=i<=N: lnP at which Q is zero

if nargin < 2 || isempty(method)
    method = 'linear';
end

if nargin < 3
    options = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 600, 'MaxIter', 100, 'Display', 'off');
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
    Q = S{i}(:, 2);
    isotherm_pp = interp1(lnP, Q, method, 'pp');
    isotherm{i} = @(x)ppval(isotherm_pp, x);
    [minlnP(i), fval, exitflag, output, jacobian] = fsolve(isotherm{i}, min(lnP), options);
    ads_pot{i} = @(y)piecewise_polynomial_ads_pot(y, isotherm_pp, minlnP(i));
    inv_ads_pot{i} = [];
end

end