function q = Langmuir_Sips(lnP, lnK, q_sat, v)
% M-site Langmuir_Sips models for single-component systems
% lnP(i), 1<=i<=N: ln(pressures)
% lnK(j), 1<=j<=M: ln(Henry's constants)
% q_sat(j), 1<=j<=M: saturation loadings
% v(j), 1<=j<M: Sip exponents

M = length(lnK);

if nargin < 4
    v = ones(M, 1);
elseif M ~= length(v)
    error('Langmuir_Sips:NonEqualParameters', 'K and v must have same number of parameters');
end

if M ~= length(q_sat)
    error('Langmuir_Sips:NonEqualParameters', 'K and q_sat must have same number of parameters');
end

q = 0;
for i = 1:M
    t = exp(lnK(i) + lnP.*v(i));
    q = q + q_sat(i).*t./(1+t);
end
end 