function [gamma, psi, x2, lnpi0] = activity_from_MC(M, S_ads_pot, S_inv_ads_pot, changefit)
% Calculate activity coefficients from binary fixed composition MC data
% INPUT:
%   M{j}(i, 1:2): component j, i-th pressure (P) & loading (Q)
%   S_ads_pot{j}: function handles for j-th component pure isotherm
% OUTPUT:
%   gamma(i, 1:2): i-th activity coefficient for i-th pressure point

% set defaults
arguments
    M cell
    S_ads_pot cell = {}
    S_inv_ads_pot cell = {}
    changefit string = 'off'
end

N = length(M); 
ndata = length(M{1});
psi = ones(ndata, 1);
lnpi0 = ones(ndata, N);

if M{1}(1,1)+M{2}(1,1) == M{1}(end,1)+M{2}(end,1) % fixed pressure campaign
    disp('FIXED PRESSURE')

    % first integrate pure component isotherm to mixture pressure
    totalP = M{1}(1,1)+M{2}(1,1);
    totalQ = M{1}(:,2)+M{2}(:,2);
    y1 = M{1}(:,1)./totalP;
    term1 = S_ads_pot{2}(log(totalP));

    % get piece-wise for loading = n_i = f(y_i)
    M_fixP_a = {[y1, M{1}(:,2)]};
    M_fixP_b = {[1-y1, M{2}(:,2)]};
    [term2_isotherm_a, min_y1_term2_a, max_y1_term2_a, M_ads_pot_term2_a] = fit_piecewise_polynomial(M_fixP_a);
    [term2_isotherm_b, min_y1_term2_b, max_y1_term2_b, M_ads_pot_term2_b] = fit_piecewise_polynomial(M_fixP_b);

    % solve for psi using both terms
    for n = 1:ndata
        fprintf('\n======Data point %d ======\n', n);
        term2_a = M_ads_pot_term2_a{1}(log(y1(n)));
        term2_b =  M_ads_pot_term2_b{1}(0) - M_ads_pot_term2_b{1}(log(1-y1(n)));
        psi(n) = term1 + term2_a + term2_b;
        fprintf('For state point %.0f, totalP is %.0f, psi is %.2f\n', n, totalP, psi(n))

    % using psi, solve for sorption pressure pi_0
        for c = 1:N
            if isempty(S_inv_ads_pot{c})
                lnpi0(n,c) = inv_adsorption_potential(psi(n), S_ads_pot{c});
            else
                lnpi0(n,c) = S_inv_ads_pot{c}(psi(n));
            end
        end
    end

else % fixed composition campaign
    disp('FIXED COMP')
    M = {sort(M{1}); sort(M{2})}; % so plot isn't criss/cross
    totalM = {[M{1}(:,1)+M{2}(:,1), M{1}(:,2)+M{2}(:,2)]};
    totalP = totalM{1}(:,1); totalQ = totalM{1}(:,2);

    % fit function for integration
    % solve for psi by integrating binary isotherm to mixture pressure
    if strcmp(changefit, 'off')
        disp('RUNNING PW FIT')
        [M_isotherm, M_minlnP, M_maxlnP, M_ads_pot, M_inv_ads_pot] = fit_piecewise_polynomial(totalM);
        % if piecewise, use:
        for i = 1:length(totalP)
            psi(i) = M_ads_pot{1}(log(totalP(i)));
        end
    else
        disp('RUNNING DS LANGMUIR FIT')
        [M_isotherm, M_minlnP, M_maxlnP, M_ads_pot, M_inv_ads_pot] = fit_Langmuir(totalM,2);
        psi = M_ads_pot{1}(log(totalP));
    end

    % using psi, solve for sorption pressure pi_0
    for n = 1:ndata
        % fprintf('For state point %.0f, psi is %.2f\n', n, psi(n))
        for c = 1:N
            if isempty(S_inv_ads_pot{c})
                lnpi0(n,c) = inv_adsorption_potential(psi(n), S_ads_pot{c});
            else
                lnpi0(n,c) = S_inv_ads_pot{c}(psi(n));
            end
        end
    end
end

% using raoult's law, solve for activity
pi0 = exp(lnpi0);
x2 = M{2}(:,2)./totalQ;
gamma = [M{1}(:,1)./(pi0(:,1).*(1-x2)), M{2}(:,1)./(pi0(:,2).*(x2))];
end