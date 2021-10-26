function lnP = Langmuir_Sips_inv_ads_pot(Psi, lnK, q_sat, v)
% 1-site Langmuir_Sips models for single-component systems
% p(i), 1<=i<=N: pressures
% lnK: ln(Henry's constants)
% q_sat: saturation loadings
% v: Sip exponents

if nargin < 4 || isempty(v)
    v = 1;
end

if length(lnK) ~= 1 || length(q_sat) ~= 1 || length(v) ~= 1
    error('Langmuir_Sips_inv_ads_pot:TooManyParameters', 'Langmuir_Sips_inv_ads_pot currently supports only 1-site Langmuir-Sips model');
end

lnP = (log(exp(Psi.*v./q_sat) - 1) - lnK) ./ v;

end
