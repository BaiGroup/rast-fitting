%% MFI-1/298K data
load('MFI-1-298K.mat');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=exp(1.50116)*1e3;  % T=298K
V_H2O=18.015e-6/0.991665;  % g/mL -> m^3/mol
water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/298)*p_sat_H2O;

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

S={[water(:,1),water(:,2)],[ethanol(:,1),ethanol(:,2)]};
M={[ethanol_water(:,6)*1000,ethanol_water(:,5)],[ethanol_water(:,3)*1000,ethanol_water(:,2)]};
Q=[ethanol_water(:,5),ethanol_water(:,2)];
z=Q./sum(Q,2);
N = length(M);
ndata = size(M{1}, 1);

%% Mode 1 or 2
% x0=[18.8249587745800,4.50258482656915,0.164809895587217,0.000188549760651445,1;19.9579736270755,8.51319839836400,0.0379042919906339,0.000228786934793868,1;20.0472128712751,8.79663421666627,0.00676836636075732,0.000288671378447038,1;19.9964587073367,8.63546802334574,0.0290109776526397,0.000242517827250023,1;20.0697841752496,8.86825616941019,0.00310915132873038,0.000298248261582597,1;20.0051269236248,8.66300293044885,0.0208613387774240,0.000257273524427223,1;20.0537857294982,8.81751143898200,0.00182818729579778,0.000301756776734216,1;20.0448286971051,8.78906535283560,0.0123430162495032,0.000275269102518928,1;19.9010742436970,8.33226180074198,0.0430738362338025,0.000221811933763693,1;18.7142998100242,4.03561965294625,0.176835048578078,0.000195704557112266,1;19.8256097787257,8.08940325353482,0.0542809105169930,0.000208921389758377,1;17.9696610841657,3.70410321681428,0.251113098944458,0.000294422887795130,1;19.0953420159099,5.62231497601938,0.132314655258204,0.000178280018723405,1];
[Q_predicted,x]=RAST_solve(M,S,@Margules,[],1,x0);
% [Q_predicted,x, Q_IAST, x_IAST, F_IAST]=RAST_solve(M,S,@Margules,[],1);
gamma=x(:,2*N:3*N-1);
coeff=fit_activity_model(@Margules, 2, z(:,1), gamma)

%% Plot activity coefficients
z_predicted=Q_predicted./sum(Q_predicted,2);
z_plot = linspace(0,1,100)';
gamma_plot = Margules(coeff,z_plot);
plot(z_predicted(:,2),gamma(:,2),'rs',z_predicted(:,2),gamma(:,1),'bo',1-z_plot,gamma_plot(:,2),'m-',1-z_plot,gamma_plot(:,1),'c-');
xlabel('$x_{\textrm{EtOH}}$'); ylabel('$\gamma$');
legend('sim\_EtOH','sim\_H2O','Margules\_EtOH','Margules\_H2O','Location','NorthEastOutside');

figure;
plot(z_predicted(:,2),log(gamma(:,2)),'rs',z_predicted(:,2),log(gamma(:,1)),'bo',1-z_plot,log(gamma_plot(:,2)),'m-',1-z_plot,log(gamma_plot(:,1)),'c-');
xlabel('$x_{\textrm{EtOH}}$'); ylabel('$\ln\gamma$');
legend('sim\_EtOH','sim\_H2O','Margules\_EtOH','Margules\_H2O','Location','NorthEastOutside');


%% Mode 3
x0=[18.8222120392256,4.49015668974766,19.9575875060930,8.51197469499720,20.0472146402693,8.79663433021242,19.9962472562403,8.63481055728581,20.0697859465553,8.86826691429560,20.0050364583461,8.66271118508367,20.0537931017180,8.81751957437988,20.0448082554387,8.78900000160829,19.9005667496345,8.33063418479486,18.7104388249224,4.02292883214846,19.8247939770017,8.08673779841672,17.9272852086422,3.69494022863267,19.0924045485707,5.61167415159729,0.154365399591610,0.0367265883087434,0.00676849253558309,0.0283721961779323,0.00311950847263224,0.0205756531872349,0.00183634919874207,0.0122782821345066,0.0415153651625246,0.166321926889219,0.0517567989118417,0.244219557585561,0.123032530947916,-20.4119003897631,-8.09496768250881];
[Q_predicted,x, Q_IAST, x_IAST, F_IAST]=RAST_solve(M,S,@Margules,[],3,x0);
% [Q_predicted,x]=RAST_solve(M,S,@Margules,[],3,[],2);
coeff=x(ndata*(2*N-1)+1:end);

%% Mode 4
x0=[18.8249587745800,4.50258482656915,19.9579736270755,8.51319839836400,20.0472128712751,8.79663421666627,19.9964587073367,8.63546802334574,20.0697841752496,8.86825616941019,20.0051269236248,8.66300293044885,20.0537857294982,8.81751143898200,20.0448286971051,8.78906535283560,19.9010742436970,8.33226180074198,18.7142998100242,4.03561965294625,19.8256097787257,8.08940325353482,17.9696610841657,3.70410321681428,19.0953420159099,5.62231497601938,-20.848654436889824,-8.088954640210344];  % MFI-1, 298K
x0=[9.83936297560183,4.97292910849719,9.84866598014940,6.09052993415287,9.90264670393793,7.40188150999622,10.1082476807207,8.48418100739620,10.2532490607674,9.24794117177083,10.3212287711267,9.65680802703133,10.3326892389181,9.72836846051617,10.3643190469011,9.93127734182888,10.3923370573292,10.1181347206316,10.3931617016872,10.1254748503674,10.4212127900103,10.3197024200992,1.06936343759260,2.96878064702337];  % LTA-0, 323K
[Q_predicted,x, Q_IAST, x_IAST, F_IAST]=RAST_solve(M,S,@Margules,[],4,x0);
% [Q_predicted,x]=RAST_solve(M,S,@Margules,[],4,[],2);
coeff=x(ndata*N+1:end);

%% Mode 5
x0=[4.80 -11.84];
[Q_predicted,x]=RAST_solve(M,S,@Margules,[],5,x0);
% [Q_predicted,x]=RAST_solve(M,S,@Margules,[],5,[],2);
coeff=x;

%% Plotting isotherms
semilogx(M{2}(:,1),ethanol_water(:,2),'rs',M{2}(:,1),ethanol_water(:,5),'bo',M{2}(:,1),Q_predicted(:,2),'md',M{2}(:,1),Q_predicted(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [molec/uc]');
legend('sim\_EtOH','sim\_H2O','RAST\_EtOH','RAST\_H2O','Location','NorthEastOutside');

%% MultiStart for solving AST equations with provided EoS
N = length(M);
ndata = size(M{1}, 1);
P = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    P(:, i) = M{i}(:, 1);
    Q(:, i) = M{i}(:, 2);
end

EoS_with_coeff = @(y)EoS(x, y);

prob=createOptimProblem('fmincon','objective',func,'x0',k','lb',k_lb,'ub',k_ub,'options',opts);
ms=MultiStart('Display','iter','TolFun',1e-6,'TolX',1e-6,'UseParallel','always');
[k,fv,ef,out,allmins]=run(ms,prob,1000);

[Q_predicted, x, F] = IAST_solve(P, [], isotherm, minlnP, EoS_with_coeff, [], 2, x0);

%% Activity coefficients
z_predicted=Q_predicted./sum(Q_predicted,2);
gamma=Margules(coeff,z_predicted(:,1:end-1))