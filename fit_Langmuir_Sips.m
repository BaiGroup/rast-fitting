function [isotherm, ads_pot, inv_ads_pot] = fit_Langmuir_Sips(S, options)
% Fit numerical interpolants for N single-component isotherms
%   S{j}(i, 1:2): component j, i-th pressure (P) & loading (Q)
%   method is passed to interp1 and defaults to 'linear'. It is used in
%     interpolating single-component isotherms.
%   options: argument for fsolve()
%
% OUTPUT:
%   isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
%     minlnP(i), 1<=i<=N: lnP at which Q is zero

if nargin < 2
    options = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 600, 'MaxIter', 100, 'Display', 'off');
end

N = length(S);
isotherm = cell(1, N);
ads_pot = cell(1, N);
inv_ads_pot = cell(1, N);
for i = 1 : N  % components
    x0 = ones(1, 3);
    func = @(x)(Langmuir_Sips(log(S{i}(:,1)), x(:, 1), x(:, 2), x(:, 3)) - S{i}(:,2));
    [param,fval,exitflag,output,jacobian] = fsolve(func, x0, options);
    isotherm{i} = @(x)Langmuir_Sips(x, param(1), param(2), param(3));
    ads_pot{i} = @(x)Langmuir_Sips_ads_pot(x, param(1), param(2), param(3));
    inv_ads_pot{i} = @(x)Langmuir_Sips_inv_ads_pot(x, param(1), param(2), param(3));
end

end