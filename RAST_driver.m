%% MFI-1/298K data
load('MFI-1-298K.mat');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=exp(1.50116)*1e3;  % T=298K
V_H2O=18.015e-6/0.991665;  % g/mL -> m^3/mol
water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/298)*p_sat_H2O;

%% MFI-1/323K data
load('MFI-1-323K.mat');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=18.165e3;  % T=323K
V_H2O=18.46e-6;  % m^3/mol
% water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/323)*p_sat_H2O;

%% Preparation
set(0,'DefaultTextInterpreter','latex');
options=optimset('Display','iter');

S={[water(:,1),water(:,2)],[ethanol(:,1),ethanol(:,2)]};
M={[ethanol_water(:,1),ethanol_water(:,3)],[ethanol_water(:,2),ethanol_water(:,4)]};
Q=[ethanol_water(:,3),ethanol_water(:,4)];
z=Q./sum(Q,2);
N = length(M);
ndata = size(M{1}, 1);

% Fitting single-component isotherms
% [isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_Langmuir_Sips(S(1));
[isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_piecewise_polynomial(S(1));
[isotherm(2), minlnP(2), maxlnP(2), ads_pot(2), inv_ads_pot(2)] = fit_piecewise_polynomial(S(2));

%% Plotting single-component isotherms
lnP_plot = linspace(min(log(S{1}(:,1)))*0.9, max(log(S{1}(:,1)))*1.1, 100);
Q_plot = isotherm{1}(lnP_plot);
semilogx(S{1}(:,1), S{1}(:,2), 'o', 'DisplayName', 'H2O');
hold on; semilogx(exp(lnP_plot), Q_plot, '-', 'DisplayName', 'H2O fitted');

figure;
lnP_plot = linspace(min(log(S{2}(:,1)))*0.9, max(log(S{2}(:,1)))*1.1, 100);
Q_plot = isotherm{2}(lnP_plot);
semilogx(S{2}(:,1), S{2}(:,2), 'o', 'DisplayName', 'EtOH')
hold on; semilogx(exp(lnP_plot), Q_plot, '-', 'DisplayName', 'EtOH fitted');

%% Mode 1 or 2 or 102
x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 1, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules, 'x0', x0);
gamma = x(:, N+1:end);
coeff = fit_activity_model(@Margules_ads, 3, [z(:,1), psi(:, 1)], gamma)
% C_ub = 0 is important because C can be arbitrarily varied without
% affecting the results
coeff2 = fit_activity_model(@Margules, 3, [z(:,1), psi(:, 1)], gamma, 'C_ub', 0)

%% Plot activity coefficients
z_plot = linspace(0, 1, 100)';
gamma_plot = Margules(coeff2, z_plot);
plot(z(:,2),gamma(:,2),'ro',z(:,2),gamma(:,1),'bo',1-z_plot,gamma_plot(:,2),'m-',1-z_plot,gamma_plot(:,1),'c-');
xlabel('$x_{\textrm{EtOH}}$'); ylabel('$\gamma$');
legend('sim\_EtOH','sim\_H2O','Margules\_EtOH','Margules\_H2O','Location','NorthEastOutside');

figure;
plot(z(:,2),log(gamma(:,2)),'ro',z(:,2),log(gamma(:,1)),'bo',1-z_plot,log(gamma_plot(:,2)),'m-',1-z_plot,log(gamma_plot(:,1)),'c-');
xlabel('$x_{\textrm{EtOH}}$'); ylabel('$\ln\gamma$');
legend('sim\_EtOH','sim\_H2O','Margules\_EtOH','Margules\_H2O','Location','NorthEastOutside');

hold on
gamma_plot_0_5 = Margules_ads(coeff, [z_plot, 0.5*ones(size(z_plot))]);
gamma_plot_5 = Margules_ads(coeff, [z_plot, 5*ones(size(z_plot))]);
gamma_plot_50 = Margules_ads(coeff, [z_plot, 50*ones(size(z_plot))]);
gamma_plot_100 = Margules_ads(coeff, [z_plot, 100*ones(size(z_plot))]);
gamma_plot_130 = Margules_ads(coeff, [z_plot, 130*ones(size(z_plot))]);
plot(1-z_plot,log(gamma_plot_0_5(:,2)),'m--','DisplayName','0.5 EtOH');
plot(1-z_plot,log(gamma_plot_0_5(:,1)),'c--','DisplayName','0.5 H2O');
plot(1-z_plot,log(gamma_plot_5(:,2)),'m:','DisplayName','5 EtOH');
plot(1-z_plot,log(gamma_plot_5(:,1)),'c:','DisplayName','5 H2O');
plot(1-z_plot,log(gamma_plot_50(:,2)),'m^','DisplayName','50 EtOH');
plot(1-z_plot,log(gamma_plot_50(:,1)),'c^','DisplayName','50 H2O');
plot(1-z_plot,log(gamma_plot_100(:,2)),'m+','DisplayName','100 EtOH');
plot(1-z_plot,log(gamma_plot_100(:,1)),'c+','DisplayName','100 H2O');
plot(1-z_plot,log(gamma_plot_130(:,2)),'ms','DisplayName','130 EtOH');
plot(1-z_plot,log(gamma_plot_130(:,1)),'cs','DisplayName','130 H2O');

%% Mode 3
x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 3, 'C_ub', +Inf, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv);
coeff=x(ndata*(2*N-1)+1:end);

%% Mode 4
x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 4, 'C_ub', 500, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv);
coeff=x(ndata*N+1:end);

%% Mode 4 with no \Psi-dependence
% Set C_ub = 0, EoS = @Margules, N_EoS_param = 3, EoS_deriv = @(x,y)0
% C_ub = 0 is important because C can be arbitrarily varied without
% affecting the results
x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 4, 'C_ub', 0, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @(x,y)0);
coeff=x(ndata*N+1:end);

%% Mode 5
x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 5, 'C_ub', +Inf, 'g_lb', -1e2, 'g_ub', 1e2, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv);
coeff = x;

%% Mode 6 or 7
x0 = [];
[gamma, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 6, 'C_ub', +Inf, 'g_lb', -1e2, 'g_ub', 1e2, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv);
coeff = x;

%% RAST fitting with fitted EoS
% x0 = [z(:,1:end-1), lnP0];
x0 = [];
[err, Q_predicted, x, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(coeff, isotherm, M, @Margules_ads, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @Margules_ads_deriv, 'x0', x0, 'options', options);

%% Plotting isotherms
semilogx(M{2}(:,1),M{2}(:,2),'rs',M{2}(:,1),M{1}(:,2),'bo',M{2}(:,1),Q_predicted(:,2),'md',M{2}(:,1),Q_predicted(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$Q$ [mol/kg]');
legend('sim\_EtOH','sim\_H2O','RAST\_EtOH','RAST\_H2O','Location','NorthEastOutside');

figure;
loglog(M{2}(:,1),(M{2}(:,2)./M{1}(:,2))./(M{2}(:,1)./M{1}(:,1)),'rs',M{2}(:,1),(Q_predicted(:,2)./Q_predicted(:,1))./(M{2}(:,1)./M{1}(:,1)),'bo');
xlabel('$p$ [Pa]'); ylabel('$S$');
legend('sim','RAST','Location','NorthEastOutside');

%% Activity coefficients
gamma1 = Margules_ads(coeff, [x(:,1), psi_IAST(:,1)])
gamma2 = Margules_ads(coeff, [x(:,1), psi_IAST(:,2)])
gamma1 - gamma2
q_ex1 = Margules_ads_deriv(coeff, [x(:, 1), psi_IAST(:,1)])
q_ex2 = Margules_ads_deriv(coeff, [x(:, 1), psi_IAST(:,2)])
q_ex1 - q_ex2
Q_predicted1 = 1./(1./Q_predicted - q_ex1)
Q_predicted2 = 1./(1./Q_predicted - q_ex2)
Q_predicted - Q_predicted1
Q_predicted - Q_predicted2