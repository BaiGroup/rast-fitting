function lnP = Langmuir_Sips_inv_ads_pot(Psi, lnK, q_sat, v)
% 1-site Langmuir_Sips models for single-component systems
% p(i), 1<=i<=N: pressures
% lnK: ln(Henry's constants)
% q_sat: saturation loadings
% v: Sip exponents

if nargin < 4
    v = 1;
end

lnP = (log(exp(Psi.*v./q_sat) - 1) - lnK) ./ v;

end
