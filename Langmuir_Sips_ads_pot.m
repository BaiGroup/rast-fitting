function y = Langmuir_Sips_ads_pot(maxlnP, lnK, q_sat, v)
% M-site Langmuir_Sips models for single-component systems
% maxlnP(i), 1<=i<=N: upper limits of integration
% lnK(j), 1<=j<=M: ln(Henry's constants)
% q_sat(j), 1<=j<=M: saturation loadings
% v(j), 1<=j<M: Sip exponents

M = length(lnK);

if nargin < 4 || isempty(v)
    v = ones(M, 1);
elseif M ~= length(v)
    error('Langmuir_Sips_ads_pot:NonEqualParameters', 'K and v must have same number of parameters');
end

if M ~= length(q_sat)
    error('Langmuir_Sips_ads_pot:NonEqualParameters', 'K and q_sat must have same number of parameters');
end

y = 0;
for i = 1:M
    y = y + q_sat(i) ./ v(i) .* log(1 + exp(lnK(i) + maxlnP.*v(i)));
end

end
