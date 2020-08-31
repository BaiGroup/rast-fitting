%% Load data 
clc

% single comp water properties: zeo name, uc mid_N mid_mu widom_mu
zeo_prop = [importdata('iast_water.txt')]';
names = zeo_prop.textdata;
uc = zeo_prop.data(:,1); 
n_sat = zeo_prop.data(:,2).*2./uc;
mid_n = mean(zeo_prop.data(:,[4,6]),2)./uc; 
mid_mu = mean(zeo_prop.data(:,[3,5]),2); 
widom_mu = zeo_prop.data(:,5);
% T = 323; R=8.314; mu_ig = -8605.407; psat=[12260, 29545]; vol=0.018233e-3;
% mid_p = psat(1) * exp((mid_mu - mu_ig)/T); 
% avgp = mean(mid_p(3,:));
% avgn = mean(mid_n(3,:));
% semilogx(mid_p(3,:), mid_n(3,:), 'ro', avgp, avgn, 'bo')

% single comp etoh data
etoh_loading = [importdata('iast_etoh.txt')]'; % Zx12
etoh_fug = [95.11, 127.2, 170, 235, 324.9, 449.1, 620.9, 999.6, 1609, ...
    2591, 4319, 7200]';

% binary data
binary_data = importdata('iast_binary.txt'); % 2*Zx11, first row is water
%% %% water/ethanol
% using single comp to predict binary isotherm
% S={[water(:,1),water(:,2)], [ethanol(:,1),ethanol(:,2)]};
% M=[ethanol_water(:,6)*1000, ethanol_water(:,3)*1000];

% S: cell array, pressure+loading columns 1-2 for each component (j)
% S{j}(i, 1:2): component j, i-th pressure & loading
for i = 1:2
    [fitted_p, fitted_loading, all_p, all_n] = getwateriso(uc(i), n_sat(i), mid_n(i), mid_mu(i), widom_mu(i), 1e-6);
    S={[fitted_p', fitted_loading'], [etoh_fug, etoh_loading(:,i)/uc(i)]};
    % M: partial pressures in mixture where col = comp and row = pressure
    % partial pressure = fugacity
    fug_mix_water = [18549, 18068, 17994, 17945, 17630, 16664, 15382,...
        13224, 9731.9, 6531.9, 965.87]';
    fug_mix_etOH = [140.29, 401.12, 1299.4, 4837.1, 10364, 15655,...
        16870, 20355, 23307, 25318, 27455]';
    M=[fug_mix_water, fug_mix_etOH];
    Q=IAST_solve(M,S);
    
    % plotting
    figure
    semilogx(M(:,2), binary_data((i*2),:)./uc(i), 'rs',...
        M(:,2), Q(:,2), 'r*', M(:,2), binary_data((i*2-1),:)./uc(i),...
        'bs', M(:,2), Q(:,1), 'b*')
    xlabel('{P} [Pa]') 
    ylabel('{N} [molec/uc]')
    ylim([0 inf])
    legend('sim\_EtOH','IAST\_EtOH','sim\_H2O','IAST\_H2O','Location','NorthEastOutside')
    title(sprintf('IAST %s Fitting', char(names(i,:))))
end

%% water isotherm fitting
% just graph water isotherm
clc
for i = 1:3
%   M = linspace(2, 10, 5);
%   K = logspace(-12, -8, 5);
%   kFixedM = zeros(1, length(M));
%   MFixedK = zeros(1, length(K));
    [fitted_p, fitted_loading, all_p, all_n] = getwateriso(uc(i), n_sat(i), mid_n(i), mid_mu(i), widom_mu(i), 1e-6);
    figure
    semilogx(fitted_p, fitted_loading, '--', all_p, all_n, 'ko')
    hold on
    [fitted_p, fitted_loading, all_p, all_n] = getwateriso(uc(i), n_sat(i), mid_n(i), mid_mu(i), widom_mu(i), 1);
    semilogx(fitted_p, fitted_loading, '*')
    legend('Scaled Fitting', 'Simulation','Unscaled Fitting', 'location','northeast')
    title(sprintf('%s Water Isotherm', char(names(i,:))))
    xlabel('P (Pa)')
    ylabel('N (molec/uc)')
%     legendcell = cellstr(num2str(M', 'M=%-d'));
%     legendcell(end+1, :) = cellstr('simulation');
%     validation = [importdata('isotherm-TIP4P-MFI-1-323K.txt')];
%     x = validation.data(:,2); y = validation.data(:,3);
%     uitable('Data', kFixedM', 'ColumnName', 'Fixed M, K = ', 'Position', [700 100 250 150]);
end

function [water_p, water_loading, all_p, all_n] = getwateriso(uc, n_sat, mid_n, mid_mu, widom_mu, scale)
% NIST: psat=12260  vol=0.018233e-3; 
psat_sim = [18.14 18.08 18.48 17.96].*1000; psat = mean(psat_sim);
vol_sim = 1./[0.974351 0.976658 0.976596 0.976051].*1e-6.*18.01528;
vol = mean(vol_sim); T = 323; R=8.314; mu_ig = -8605.407;

% get henry's constant first: k_H = exp(-(mu_ex-mu_ig)/kbT) / uc
k_h = exp(-(widom_mu-mu_ig)/T)/uc; % PUT /UC BACK

% sips model: q_e = (q_lf*(K*c_e)^m_lf) / (1 + K*c_e)^m_lf
% parameters: q_lf = saturation loading via N_mid iterative script, 
% K = equil constant for solid, m_lf = heterogenous parameter
% convert mu to pressure: p = psat*fugacity = psat*exp((mu-mu_ig)/T)
mid_p = psat * exp((mid_mu - mu_ig)/T);
% low pressure: pick p (log scale) so q = k*p is small
pick_p = logspace(7.5, 9, 8); % for liquid [=] Pa
tr_low_p = psat * exp((pick_p - psat)*(vol/(R*T))); % for vapor
tr_low_n = k_h*tr_low_p;
filter = find(tr_low_n < 0.2);
low_p = pick_p(filter);
low_n = tr_low_n(filter);
all_p = [low_p, mid_p]; 
scaled_all_p = all_p.*1e-6; % for fitting, need to scale Pa -> MPa
all_n = [low_n, mid_n];

% regression fitting: model, initial guess, fit
% para1 = sat loading, para2 = K, para3 = M
% N = nsat*K*P^M / (1 + K*P^M)
fun = @(para, all_p) n_sat.*(para(1)*scale).*all_p.^para(2)./(1+(para(1)*scale)*all_p.^para(2));
%fun = @(para, all_p) n_sat*para(1)*all_p.^M./(1+para(1)*all_p.^M);
%fun = @(para, all_p) n_sat*K*all_p.^para(1)./(1+K*all_p.^para(1));
options = optimoptions('lsqcurvefit', 'Display','iter', 'MaxFunEvals', 300, 'StepTolerance', 1e-4, 'FunctionTolerance', 1e-5);
format shortEng;
para0 = [1/mid_p(1) 0.5]; lb = []; ub = []; % OG GUESS = 1/MIDP
para = lsqcurvefit(fun, para0, scaled_all_p, all_n, lb, ub, options)
water_p = logspace(7, 11);
water_p_scaled = water_p*1e-6; % for fitting, P [=] MPa
water_loading = fun(para, water_p_scaled);
end
