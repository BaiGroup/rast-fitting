function [err, Q_predicted] = RAST_func_per_point(x, isotherm, minlnP, lnP, Q, mode)
% Real Adsorbed Solution Theory constraint function
% x(i), 1<=i<=N: ln p_i^0
% x(N+i), 1<=i<=N-1: z_i for component i in the absorbed phase
% x(2N-1+i), 1<=i<=N: i-th activitiy coefficient
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% lnP(i), 1<=i<=N: bulk phase pressure
% Q(i), 1<=i<=N: loading
% mode: 1 uses P*y1/P_i^0 - x_i*gamma_i = 0
%       2 uses ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i), gamma_i cannot be
%       negative during optimizations

N = length(minlnP);

if nargin < 6
    mode = 1;
end

IAST_err = IAST_func(x(1:2*N-1), isotherm, minlnP, lnP, x(2*N:end), [], mode);

z = x(N+1:2*N-1);
z = [z, 1-sum(z)];
Q_predicted = 0;
for j=1:N
    Q_predicted = Q_predicted + z(j) / isotherm{j}(x(j));
end
Q_predicted = 1 / Q_predicted * z;

err = [IAST_err, ((Q_predicted-Q)./Q)];

end