function [err, lnP0, psi] = IAST_func(x, isotherm, minlnP, lnP, EoS, tol, mode)
% AST equations for N-component systems
% x(i), 1<=i<=N-1: z_i for component i in the absorbed phase
% Mode 1 or 2 or -2:
%   x(N-1+i), 1<=i<=N: ln p_i^0
% Mode 3 or 4:
%   x(N), adsorption potential \Omega = \frac{\pi A}{RT}
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% lnP(i), 1<=i<=N: bulk phase pressure
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N-1, Omega])
% tol: precision
% mode: 1 or 3 uses P*y1/P_i^0 - x_i*gamma_i = 0
%       2 or 4 uses ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i), gamma_i cannot
%       be negative during optimizations

N = length(minlnP);
z = [x(1:N-1), 1-sum(x(1:N-1))];

if nargin<5 || isempty(EoS)
    EoS = @(x)ones(1, N);
end

if nargin < 6 || isempty(tol)
    tol = 1e-3;
end

if nargin < 7 || isempty(mode)
    mode = 1;
end

if mode == 3 || mode == 4
    psi = ones(N) * x(N);
    lnP0 = solve_adsorption_potential(psi(N), isotherm, minlnP, tol);
elseif mode == 1 || mode == 2 || mode == -2
    lnP0 = x(N:2*N-1);
    err = zeros(1, 2*N-1);
    psi = zeros(1, N);
    funcN = isotherm{N};
    psi(N) = adsorption_potential(funcN, minlnP(N), lnP0(N), tol);
    for i = 1:N-1
        psi(i) = adsorption_potential(isotherm{i}, minlnP(i), lnP0(i), tol);
        err(N+i) = psi(N) - psi(i);
    end
else
    error('IAST_func:UnknownMode','Mode parameter outside known choices.')
end

gamma = EoS([x(1:N-1), psi(N)]);
if mode == 1 || mode == 3
    err(1:N)=exp(lnP-lnP0)-z.*gamma;
elseif mode == 2 || mode == 4 || mode == -2
    err(1:N)=lnP-lnP0-log(z.*gamma);
end
end