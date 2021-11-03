% AST_func Adsorbed Solution Theory (AST) equations for N-component systems
%
% Syntax
%   y = AST_func(x, isotherm, minlnP, logP, tol)
%   y = AST_func(x, isotherm, minlnP, logP, EOS, tol)
%
% Description
%   x(i), 1<=i<=N-1: absorbed-phase mole fraction for the i-th component
%   Mode 1 or 2 or 102:
%     x(N-1+i), 1<=i<=N: ln p_i^0
%   Mode -1 or -2:
%     x(N), adsorption potential \Psi = \frac{\pi A}{RT}
%   isotherm(i), 1<=i<=N: function handle for isotherm of the i-th component, Q(lnP)
%   minlnP(i), 1<=i<=N: lower limit of integration domain,
%        i.e., ln(pressure) at which the amount of adsorption is zero
%   lnP(i), 1<=i<=N: ln(bulk phase pressure)
%   tol: precision
%   EoS: A function handle that takes N mole fractions and adsorption
%       potential and returns N activity coefficients. EoS is optional and
%       defaults to 1 (ideal adsorbed solution)
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N-1, Psi])
%
%   y(i), 1<=i<=N:
%     mode 1 or -1 uses P*y1/P_i^0 - x_i*gamma_i = 0
%     2 or -2 uses ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i), gamma_i cannot
%     be negative during optimizations
%   y(N+i), 1<=i<=N-1: Raoult's law Pi_i^0 (p_i^0) = Pi_N^0 (p_N^0)
%                  evaluated using trapezoidal expression
%
%   See Bai et al., Langmuir 28 (2012) 15566 for reference to equations
%   and symbols.
%
% See also IAST_solve, trapz

function [err, lnP0, psi] = IAST_func(x, lnP, varargin)
if nargin < 2 || rem(nargin,2) ~= 0
    error('IAST_func:Number of arguments incorrect');
end

N = numel(lnP);

names          = {'isotherm'; 'minlnP'; 'EoS'; 'mode'; 'tol'; 'ads_pot'; 'inv_ads_pot'};
default_values = {        [];       [];    [];      1;  1e-5;        []; []};
opt_args = process_variable_arguments(names, default_values, varargin);
isotherm = opt_args.('isotherm');
minlnP = opt_args.('minlnP');
EoS = opt_args.('EoS');
mode = opt_args.('mode');
tol = opt_args.('tol');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');

if isempty(isotherm)
    isotherm = cell(1, N);
end
if isempty(ads_pot)
    ads_pot = cell(1, N);
end
if isempty(inv_ads_pot)
    inv_ads_pot = cell(1, N);
end

z = [x(1:N-1), 1-sum(x(1:N-1))];

if isempty(EoS)
    EoS = @(y)ones(1, N);
end

if mode == -1 || mode == -2
    err = zeros(1, N);
    lnP0 = zeros(1, N);
    psi = ones(1, N) * x(N);
    for i = 1:N
        if isempty(inv_ads_pot{i})
            if isempty(ads_pot{i})
                if isempty(isotherm{i}) || isempty(minlnP)
                    error('IAST_func:isotherm and minlnP must be provided without inv_ads_pot');
                end
                ads_pot{i} = @(y)adsorption_potential(y, isotherm{i}, minlnP(i), tol);
            end
            inv_ads_pot{i} = @(y)inv_adsorption_potential(y, ads_pot{i});
        end
        lnP0(i) = inv_ads_pot{i}(psi(i));
    end
elseif mode == 1 || mode == 2 || mode == 102
    for i = 1:N
        if isempty(ads_pot{i})
            if isempty(isotherm{i}) || isempty(minlnP)
                error('IAST_func:isotherm and minlnP must be provided without ads_pot');
            end
            ads_pot{i} = @(y)adsorption_potential(y, isotherm{i}, minlnP(i), tol);
        end
    end
    
    err = zeros(1, 2*N-1);
    lnP0 = x(N:2*N-1);
    psi = zeros(1, N);
    psi(N) = ads_pot{N}(lnP0(N));
    for i = 1:N-1
        psi(i) = ads_pot{i}(lnP0(i));
        err(N+i) = psi(N) - psi(i);
    end
else
    error('IAST_func:UnknownMode','Mode parameter outside known choices.')
end

gamma = EoS([x(1:N-1), mean(psi)]);
if mode == 1 || mode == -1
    err(1:N) = exp(lnP-lnP0) - z.*gamma;
elseif mode == 2 || mode == -2 || mode == 102
    err(1:N) = lnP - lnP0 - log(z.*gamma);
end
% err(~isfinite(err)) = 1e6;  % protect z.*gamma being too small (mode 2) or gamma too large (mode 1)
end