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

%% LTA-0/323K data
load('LTA-0-323K.mat');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=18.165e3;  % T=323K
V_H2O=18.46e-6;  % m^3/mol
% water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/323)*p_sat_H2O;

%% FAU-2/323K data
load('FAU-2-323K.mat');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=18.165e3;  % T=323K
V_H2O=18.46e-6;  % m^3/mol
% water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/323)*p_sat_H2O;

%% Preparation
set(0,'DefaultTextInterpreter','latex');
options=optimset('Display','iter');

S={[water(:,1),water(:,2)],[ethanol(:,1),ethanol(:,2)]};
M={[ethanol_water(:,6)*1000,ethanol_water(:,5)],[ethanol_water(:,3)*1000,ethanol_water(:,2)]};
Q=[ethanol_water(:,5),ethanol_water(:,2)];
z=Q./sum(Q,2);
N = length(M);
ndata = size(M{1}, 1);

% Fitting single-component isotherms
[isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_Langmuir_Sips(S(1));
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

%% Mode 1 or 2
% x0 = [9.54659916507525,4.45822504852441,2.15015220983416,1.03285147519273;10.7092938249809,8.49842481784590,2.46727400686816,1.01342556124652;10.8598606883720,8.94395738832644,2.92892137665078,0.862798610037701;10.8013182779368,8.77124933794977,2.81583193277256,0.869093701483178;10.8660845799886,8.96227898971391,1.00984319158690,0.910873219810828;10.7824636368829,8.71549854469242,2.57892239958312,0.949060084134431;10.9331442608880,9.15925089128384,2.88408028020013,0.710472042871600;10.8518734104143,8.92044227889966,2.80653371721469,0.876490902667201;10.6354780724902,8.27801822024613,2.13091117211549,1.06051664853207;9.46847853923580,3.96055003075170,2.02507717540426,1.07826261415205;10.5764030026928,8.09877533994039,2.13590409238473,0.991620827971923;9.41428645344238,3.51197080883258,1.74561738840598,1.16350873586547;9.83768171256791,5.69098927267615,2.24765635120982,0.910258115420987];  % MFI-1, 298K
% x0 = [9.19074037382604,5.09736861287326,3.03271843050511,2.27915890872639;9.22930125758334,6.06111397875104,4.36054648137393,1.57614414145470;9.28511251961734,7.49244696500350,4.29082822781169,1.18543463456380;9.75903510736948,10.1366045075341,2.57077962210034,0.296151159422557;9.58400894142692,9.37039980212221,3.37120938170920,1.36585454540269;9.81803704095643,10.3845083678447,2.68370034335134,0.713789089929751;9.80180682142322,10.3169912086654,2.74225698564950,0.796137246747622;9.84022393484934,10.4765016692427,2.44234340204493,0.797246040196091;9.82150072042286,10.3989549906469,2.23261466089940,0.918720382520571;9.79039123665228,10.2690571146685,2.17338756602035,1.04612143815037;9.78967192063356,10.2660260838746,1.00092129902735,1.00201175753042];  % FAU-2, 323K
x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 1, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules, 'x0', x0);
gamma = x(:, N+1:end);
coeff2 = fit_activity_model(@Margules, 2, z(:,1), gamma)
coeff = fit_activity_model(@Margules_ads, 3, [z(:,1), psi(:, 1)], gamma)

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
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 3, 'C_ub', 1e-3, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv);
coeff=x(ndata*(2*N-1)+1:end);

%% Mode 4
% x0 = [9.56932133693453,4.58066201273206,10.7082849511984,8.49498775552621,10.8136012238689,8.80710044581727,10.7775612220731,8.70082664812700,10.8398067691370,8.88440817552311,10.7573893791015,8.64146621420108,10.8190808116371,8.82396736941342,10.8138929966375,8.80816989694733,10.6028098033075,8.17959060389455,9.47753850237952,4.01646853415337,10.6078313303063,8.19611051552091,9.42404124399680,3.57813075474896,9.85707396913871,5.75919521857211,3.89521912479900,1.15354151246764,0.0352689984934897]; % MFI-1, 298K
x0 = [9.21176887922679,5.20543227719186,9.23720988318593,6.17303125509601,9.28722972181409,7.51946624358724,9.67478076237356,9.77457287521604,9.73071309545655,10.0158456845029,9.78217417546879,10.2343714235391,9.79667418222772,10.2952652244416,9.82099728592047,10.3967243392332,9.83977853762424,10.4745099104085,9.83059983986975,10.4365418255004,9.79232267420267,10.2770168186431,5.83384765614761,-1.19886573258603,0.00994314902757180];  % FAU-2, 323K
% x0=[9.83936297560183,4.97292910849719,9.84866598014940,6.09052993415287,9.90264670393793,7.40188150999622,10.1082476807207,8.48418100739620,10.2532490607674,9.24794117177083,10.3212287711267,9.65680802703133,10.3326892389181,9.72836846051617,10.3643190469011,9.93127734182888,10.3923370573292,10.1181347206316,10.3931617016872,10.1254748503674,10.4212127900103,10.3197024200992,1.06936343759260,2.96878064702337];  % LTA-0, 323K
% x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 4, 'C_ub', 1e-2, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv);
coeff=x(ndata*N+1:end);

%% Mode 5
% x0 = [34.1512916754137,0.797662310980325,0.00999386469064226];  % MFI-1, 298K
% x0 = [5.83384765614761,-1.19886573258603,0.00994314902757180];  % FAU-2, 323K
x0 = [];
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, 'mode', 5, 'C_ub', 1e-2, 'g_lb', -1e2, 'g_ub', 1e2, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv);
coeff = x;

%% RAST fitting with fitted EoS
x0 = [z(:,1:end-1), lnP0];
[err, Q_predicted, x, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(coeff, isotherm, M, @Margules_ads, 'mode', 1, 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @Margules_ads_deriv, 'x0', x0, 'options', options);

%% Plotting isotherms
semilogx(M{2}(:,1),M{2}(:,2),'rs',M{2}(:,1),M{1}(:,2),'bo',M{2}(:,1),Q_predicted(:,2),'md',M{2}(:,1),Q_predicted(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [molec/uc]');
legend('sim\_EtOH','sim\_H2O','RAST\_EtOH','RAST\_H2O','Location','NorthEastOutside');

%% Activity coefficients
gamma1 = Margules_ads(coeff, [z(:,1:end-1), psi(:,1)])
gamma2 = Margules_ads(coeff, [z(:,1:end-1), psi(:,2)])
gamma1 - gamma2
q_ex1 = Margules_ads_deriv(coeff, [z(:, 1:end-1), psi(:,1)])
q_ex2 = Margules_ads_deriv(coeff, [z(:, 1:end-1), psi(:,2)])
q_ex1 - q_ex2
Q_predicted1 = 1./(1./Q_predicted - q_ex1)
Q_predicted2 = 1./(1./Q_predicted - q_ex2)
Q_predicted - Q_predicted1
Q_predicted - Q_predicted2