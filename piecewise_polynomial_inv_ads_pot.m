function lnP = piecewise_polynomial_inv_ads_pot(psi, isotherm_pp, minlnP)
% Analytical integration of piecewise polynomials to compute adsorption
% potential. Only works for linear functions for now due to extrapolations
% and the use of the quadratic formula.

% coeffs(j, 1) * [x - breaks(j)] + coeffs(j, 2) for interval [breaks(j), breaks(j+1)]
[breaks, coeffs, num_interval, num_coeffs] = unmkpp(isotherm_pp);

if num_coeffs ~= 2
    error('piecewise_polynomial_inv_ads_pot: Only works for linear functions for now');
end

% Extend in both directions so that we can do extrapolations
coeffs = [coeffs(1, 1), coeffs(1,1)*(minlnP-breaks(1))+coeffs(1,2); coeffs];
breaks = [minlnP, breaks];

% integration constants
prefactor = 1./(num_coeffs+1-[1:num_coeffs]);
exponents = num_coeffs:-1:1;

y = psi;
for j = 1:num_interval
    delta = sum(coeffs(j, :) .* prefactor .* ((breaks(j+1) - breaks(j)).^exponents));
    if delta >= y || j == num_interval
        a = coeffs(j, 1) / 2;
        b = coeffs(j, 2) - coeffs(j, 1) * breaks(j);
        c = coeffs(j, 1) * breaks(j)^2 / 2 - coeffs(j, 2) * breaks(j) - y;
        lnP = (-b + sqrt(abs(b^2 - 4*a*c)))/2/a;
        break;
    else
        y = y - delta;
    end
end

end