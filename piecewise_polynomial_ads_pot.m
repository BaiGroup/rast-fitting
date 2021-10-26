function y = piecewise_polynomial_ads_pot(maxlnP, isotherm_pp, minlnP)
% Analytical integration of piecewise polynomials to compute adsorption
% potential

% coeffs(j, 1) * [x - breaks(j)] + coeffs(j, 2) for interval [breaks(j), breaks(j+1)]
[breaks, coeffs, num_interval, num_coeffs] = unmkpp(isotherm_pp);

% Extend in both directions so that we can do extrapolations
coeffs = [coeffs(1, 1), coeffs(1,1)*(minlnP-breaks(1))+coeffs(1,2); coeffs; coeffs(end, 1), coeffs(end, 1)*(breaks(end)-breaks(end-1))+coeffs(end, 2)];
breaks = [minlnP, breaks, maxlnP];

% integration constants
prefactor = 1./(num_coeffs+1-[1:num_coeffs]);
exponents = num_coeffs:-1:1;

y = 0;
for j = 1:num_interval+2
    if breaks(j) >= maxlnP
        break;
    elseif breaks(j+1) < maxlnP
        upper_limit = breaks(j+1);
    else
        upper_limit = maxlnP;
    end
    y = y + sum(coeffs(j, :) .* prefactor .* ((upper_limit - breaks(j)).^exponents));
end

end