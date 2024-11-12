%% set-up
clc;
clear all; close all;

% run("startup.m")
% fixed_comp = 0.5;
% y2=fixed_comp;
% zeo_label{1} = 'BEA-1';
fixed_comp = [0.1 0.3 0.5 0.7 0.9]; %y2
zeo_label = importdata('data/zeo185.txt');

for i = [42 75 113]

    % if strcmp(zeo_label{i}, 'MFI-1') | strcmp(zeo_label{i}, 'DDR-1') | ...
    %         strcmp(zeo_label{i}, 'CHA-1') | strcmp(zeo_label{i}, 'BEA-1') ...
    %         | strcmp(zeo_label{i}, 'FAU-2')
    %     continue
    % end

    fprintf('STARTING %s\n', zeo_label{i})

    % fixed comp campaign
    x1_all = cell(1, length(fixed_comp));
    gMC_all = cell(1, length(fixed_comp));
    psiMC_all = cell(1, length(fixed_comp));
    totalP_all = cell(1, length(fixed_comp));
    lnpi0_all = cell(1, length(fixed_comp));
    M1_all = cell(1, length(fixed_comp));
    M2_all = cell(1, length(fixed_comp));
    conf95_all = cell(1, length(fixed_comp));

    for c = 1:length(fixed_comp)
        y2=fixed_comp(c);
        % load unary + Campaign B data
        % align binary data
        align_data = [12, 14, 15, 15, 16]; % in order, y2 = 0.1 --> 0.9
        start_data = align_data(end) - align_data(c) + 1;
        water = importdata(sprintf('data/indsim/%s/indsim_water_avg.txt', zeo_label{i}));
        if strcmp(zeo_label{i}, 'MFI-1') | strcmp(zeo_label{i}, 'DDR-1')
            etoh = importdata(sprintf('data/indsim/%s/indsim_etoh_avg.txt', zeo_label{i}));
            p_etoh = etoh.data(:,1);
            q_etoh = etoh.data(:,2);
            conf_etoh = etoh.data(:,3);
            if y2 == 0.5
                etoh_water = importdata(sprintf('data/indsim/%s/indsim_CampB_avg.txt', zeo_label{i}));
                M={[etoh_water.data(start_data:end,1)-etoh_water.data(start_data:end,end), etoh_water.data(start_data:end,3)], ...
                    [etoh_water.data(start_data:end,end), etoh_water.data(start_data:end,4)]};
                conf_binary = [etoh_water.data(start_data:end,5), etoh_water.data(start_data:end,6)];
            else
                etoh_water = importdata(sprintf('data/indsim/%s/extsim_CampB_%.1f_avg.txt', zeo_label{i}, y2));
                M={[etoh_water.data(start_data:end,3)-etoh_water.data(start_data:end,end), etoh_water.data(start_data:end,4)], ...
                    [etoh_water.data(start_data:end,end), etoh_water.data(start_data:end,5)]};
                conf_binary = [etoh_water.data(start_data:end,6), etoh_water.data(start_data:end,7)];
            end
            conf95 = {water.data(:,4), conf_etoh, conf_binary(:,1), conf_binary(:,2)};
        else
            etoh = importdata(sprintf('data/indsim/%s/sim_etoh.txt', zeo_label{i}));
            etoh_water = importdata(sprintf('data/indsim/%s/extsim_CampB_%.1f_avg.txt', zeo_label{i}, y2));
            M={[etoh_water.data(start_data:end,3)-etoh_water.data(start_data:end,end), etoh_water.data(start_data:end,4)], ...
                [etoh_water.data(start_data:end,end), etoh_water.data(start_data:end,5)]};
            conf_binary = [etoh_water.data(start_data:end,6), etoh_water.data(start_data:end,7)];
            p_etoh = etoh.data(:,2);
            q_etoh = etoh.data(:,3);
            conf95 = {water.data(:,4), [], conf_binary(:,1), conf_binary(:,2)};
        end
        
        [p_sorted, sort_idx] = sort(p_etoh);
        S={[water.data(:,2), water.data(:,3)],[p_sorted, q_etoh(sort_idx)]};
        [isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_piecewise_polynomial(S);
        N = length(M); ndata = size(M{1}, 1);

        % calculate activities at each state point
        [gammaMC, psiMC, x2_campB, lnpi0] = activity_from_MC(M, ads_pot, inv_ads_pot);
        totalP = M{1}(:,1)+M{2}(:,1);

        % save all
        x1_all{c} = 1-x2_campB;
        gMC_all{c} = gammaMC;
        psiMC_all{c} = psiMC;
        totalP_all{c} = totalP;
        lnpi0_all{c} = lnpi0;
        M1_all{c} = M{1}; M2_all{c} = M{2};
        conf95_all{c} = conf95;
    end

    % multiple comps
    x1_all_combo = cat(1, x1_all{:});
    gMC_all_combo = cat(1, gMC_all{:});
    psiMC_all_combo = cat(1, psiMC_all{:});
    totalP_all_combo = cat(1,totalP_all{:});
    lnpi0_all_combo = cat(1,lnpi0_all{:});
    M1_all_combo = cat(1,M1_all{:}); M2_all_combo = cat(1,M2_all{:});

    % % save MC data to text file
    % fid = fopen('isotherm.txt', 'w');
    % fprintf(fid, '#%s\n#Unary Water\n#f_H20[Pa] Q[mol/kg]\n', zeo_label{i});
    % formatSpec = '%0.8e %0.4g\n'; fprintf(fid, formatSpec, S{1}');
    % fprintf(fid, '#Unary Ethanol\n#f_EtOH[Pa] Q[mol/kg]\n');
    % formatSpec = '%0.4e %0.4g\n'; fprintf(fid, formatSpec, S{2}');
    % fprintf(fid, '#Binary\n#f_H20[Pa] f_EtOH[Pa] Q_H20[mol/kg] Q_EtOH[mol/kg]\n');
    % formatSpec = '%0.2e %0.2e %0.4g %0.4g\n'; fprintf(fid, formatSpec, [M1_all_combo(:,1), M2_all_combo(:,1), M1_all_combo(:,2), M2_all_combo(:,2)]');
    % fclose(fid); movefile('isotherm.txt', sprintf('data/%s/indsim', zeo_label{i}));
    % fprintf('Saved to text file\n')

    fprintf('DONE WITH %s CAMPAIGN B GAMMAS\n', zeo_label{i})

    % get binary data and gammas from campaign A
    binary_campA = importdata(sprintf('data/indsim/%s/extsim_CampA_avg.txt', zeo_label{i}));
    M_campA{1}=[binary_campA.data(:,2), binary_campA.data(:,4)];
    M_campA{2}=[binary_campA.data(:,3), binary_campA.data(:,5)];
    totalP_campA = M_campA{1}(:,1)+M_campA{2}(:,1);
    [gammaMC_campA, psiMC_campA, x2_campA, lnpi0_campA] = activity_from_MC(M_campA, ads_pot, inv_ads_pot);
    fprintf('DONE WITH %s CAMPAIGN A GAMMAS\n', zeo_label{i})

    % save MC data to text file
    fid = fopen('isotherm_CampA.txt', 'w');
    fprintf(fid, '#%s\n#Unary Water\n#f_H20[Pa] Q[mol/kg]\n', zeo_label{i});
    formatSpec = '%0.9g %0.9g\n'; fprintf(fid, formatSpec, S{1}');
    fprintf(fid, '#Unary Ethanol\n#f_EtOH[Pa] Q[mol/kg]\n');
    formatSpec = '%0.9g %0.9g\n'; fprintf(fid, formatSpec, S{2}');
    fprintf(fid, '#Binary\n#f_H20[Pa] f_EtOH[Pa] Q_H20[mol/kg] Q_EtOH[mol/kg]\n');
    formatSpec = '%0.9g %0.9g %0.9g %0.9g\n'; fprintf(fid, formatSpec, [M_campA{1}(:,1), M_campA{2}(:,1), M_campA{1}(:,2), M_campA{2}(:,2)]');
    fclose(fid); movefile('isotherm_CampA.txt', sprintf('data/%s/indsim', zeo_label{i}));
    fprintf('Saved to text file\n')

    % combine campaign A and B data for fitting
    x1_campAB_combo = cat(1, x1_all_combo, 1-x2_campA);
    gMC_campAB_combo = cat(1, gMC_all_combo, gammaMC_campA);
    psiMC_campAB_combo = cat(1, psiMC_all_combo, psiMC_campA);
    totalP_campAB_combo = cat(1, totalP_all_combo, totalP_campA);
    lnpi0_campAB_combo = cat(1,lnpi0_all_combo, lnpi0_campA);
    M1_campAB_combo = cat(1, M1_all_combo, M_campA{1});
    M2_campAB_combo = cat(1, M2_all_combo, M_campA{2});
    fprintf('DONE WITH %s COMBINING GAMMAS\n', zeo_label{i})
    % save(sprintf('data/%s/indsim.mat', zeo_label{i}));

    % mode 1

    gamma_inf = isinf(gMC_campAB_combo);
    gMC_fit = gMC_campAB_combo; gMC_fit(gamma_inf,:)=[];
    x1_fit = x1_campAB_combo; x1_fit(gamma_inf, :) = [];
    psiMC_fit = psiMC_campAB_combo; psiMC_fit(gamma_inf, :) = [];

    [coeffpsi, errpsi, residualpsi] = fit_activity_model(@Margules_ads, 3, [x1_fit, psiMC_fit], gMC_fit);
    [coeff, err, residual] = fit_activity_model(@Margules, 3,  [x1_fit, psiMC_fit], gMC_fit, 'C_ub', 0);
    err_mo1 = err/numel(residual);
    errpsi_mo1 = errpsi/numel(residual);
    % 
    % fprintf('DONE FITTING %s MC GAMMAS (with psi): coeff %.2f %.2f %.2f\n', zeo_label{i}, coeffpsi)
    % fprintf('DONE FITTING %s MC GAMMAS (no psi): coeff %.2f %.2f %.2f\n', zeo_label{i}, coeff)
    % 
    % txt file
    fid = fopen('Method1-2.txt', 'w');
    fprintf(fid, '#%s\n#Method 1 Margules: A1 = %0.4g, A2 = %0.4g, C = %0.4g, err_regression = %0.3f\n', zeo_label{i}, coeffpsi, errpsi_mo1);
    fprintf(fid, '#Method 2 Margules: A1 = %0.4g, A2 = %0.4g, C = %0.0f, err_regression = %0.3f\n', coeff, err_mo1);
    fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] gamma_H20 gamma_EtOH lnP0_H2O lnP0_EtOH\n');
    formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g\n';
    fprintf(fid, formatSpec, [M1_campAB_combo(:,1), M2_campAB_combo(:,1), gMC_campAB_combo, lnpi0_campAB_combo]');
    fclose(fid); movefile('Method1-2.txt', sprintf('data/%s/indsim', zeo_label{i}));
    fprintf('Saved to text file\n')

    save(sprintf('data/%s/indsim.mat', zeo_label{i}));

end

%% mode 4

clear all; load(sprintf('data/%s/indsim.mat', 'USI-0'));

% x0=[];
% for idx_x0 = 1:N*length(lnpi0_campAB_combo)
%     if ~mod(idx_x0, 2) == 1
%         x0(idx_x0) = lnpi0_campAB_combo(i, 2);
%     else
%         x0(idx_x0) = lnpi0_campAB_combo(i, 1);
%     end
% end
for row = 1:length(lnpi0_campAB_combo)
    x0(N*row) = lnpi0_campAB_combo(row, 2);
    x0(N*row-1) = lnpi0_campAB_combo(row, 1);
end
x0=x0';
for idx_xappend = 1:length(coeffpsi)
    x0(end+1) = coeffpsi(idx_xappend);
end
x0(end-2:end)=[1 0.5 10];
% x0=[cat(1, [lnpi0_campAB_combo(:,1); lnpi0_campAB_combo(:,2)], coeffpsi')];
options=optimset('Display','iter', 'TolFun', 4e-5, 'TolX', 1e-5);
[Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve({M1_campAB_combo, M2_campAB_combo}, S, 'mode', 4, 'C_ub', 1e2, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv, 'options', options);
coeffpsi_mo4=x(length(M1_campAB_combo)*N+1:end);
% x_initmo4 = x;

fprintf('DONE WITH %s MARG FITTING (mode 4, psi): coeff %.2f %.2f %.2f\n', zeo_label{i}, coeffpsi_mo4)

x0(end-2:end)=[1 1 0];
% get coeffs: mode 4, no psi dependence
options=optimset('Display','iter', 'TolFun', 4e-5, 'TolX', 1e-5);
[Q_predicted2, x2, err2, lnP02, psi2, Q_IAST2, x_IAST2, err_IAST2, lnP0_IAST2, psi_IAST2] = RAST_solve({M1_campAB_combo, M2_campAB_combo}, S, 'mode', 4, 'C_ub', 0, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @(x,y)0);
coeff_mo4=x2(length(M1_campAB_combo)*N+1:end);

fprintf('DONE WITH %s MARG FITTING (mode 4, no psi): coeff %.2f %.2f %.2f\n', zeo_label{i}, coeff_mo4)
% save(sprintf('data/%s/indsim_campAB_coeffmo4.mat', zeo_label{i}));

%% txt files
save(sprintf('data/%s/indsim_campAB_coeffmo4.mat', zeo_label{i}));

fid = fopen('Method3-4.txt', 'w'); 
fprintf(fid, '#%s\n#Method 3 Margules: A1 = %0.4g, A2 = %0.4g, C = %0.4g, err_regression = %0.3f\n', zeo_label{i}, coeffpsi_mo4, err);
fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] psi_H2O[mol/kg] psi_EtOH lnp0_H2O lnp0_EtOH\n');
formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g\n';
fprintf(fid, formatSpec, [M1_campAB_combo(:,1), M2_campAB_combo(:,1), psi, lnP0]');
fprintf(fid, '#Method 4 Margules: A1 = %0.4g, A2 = %0.4g, C = %0.0f, err_regression = %0.3f\n', coeff_mo4, err2);
fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] psi_H2O[mol/kg] psi_EtOH lnp0_H2O lnp0_EtOH\n');
formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g\n';
fprintf(fid, formatSpec, [M1_campAB_combo(:,1), M2_campAB_combo(:,1), psi2, lnP02]');
fclose(fid); movefile('Method3-4.txt', sprintf('data/%s/indsim', zeo_label{i}));
fprintf('Saved %s data to text file\n', zeo_label{i}); pause(2);

%% save mode1+4+bulk+IAST predictions to text file

% mode = 4;
% y2 = 0.5;
% k = find(fixed_comp==y2);
% x0=[];
% x0=[x1_all{k}, lnpi0_all{k}];
% options=optimset('Display','iter');
% i=97; % zeo_label{1} = 'BEA-1';
% zeo_label = importdata('data/zeo185.txt');
% load(sprintf('data/%s/indsim_campAB_coeffmo4.mat', zeo_label{i}));

for k = 5

    % mode = 1;
    y2=fixed_comp(k);
    % x0=x;
    % x0=[x1_all{k}, lnpi0_all{k}];
    % options=optimset('Display','iter');

    % % mode 1 predictions, with psi
    % [err, RAST_backcalc, x, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(coeffpsi, isotherm, {M1_all{k}, M2_all{k}}, @Margules_ads, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @Margules_ads_deriv, 'x0', x0, 'options', options);
    % % extract psi and lnpi0 from err_IAST + sum
    % err_psi = sum(err_IAST(:, N+1));
    % err_lnpi0 = sum(err_IAST(:, 1:N), 1);
    % % save to text file
    % fid = fopen(sprintf('%s-mode%.0fPSI-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), 'w');
    % fprintf(fid, '#%s\n#Margules: A1 = %0.4g, A2 = %0.4g, C = %0.4g, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, coeffpsi, err, err_lnpi0, err_psi);
    % fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
    % formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
    % fprintf(fid, formatSpec, [M1_all{k}(:,1), M2_all{k}(:,1), RAST_backcalc, psi_IAST, lnP0_IAST]');
    % fclose(fid); movefile(sprintf('%s-mode%.0fPSI-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), sprintf('data/%s/indsim', zeo_label{i}));
    % fprintf('Saved to text file\n')

    % % mode 1 predictions, without psi
    % [err, RAST_backcalc, x, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(coeff, isotherm, {M1_all{k}, M2_all{k}}, @Margules, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @(x,y)0, 'x0', x0, 'options', options);
    % % extract psi and lnpi0 from err_IAST + sum
    % err_psi = sum(err_IAST(:, N+1));
    % err_lnpi0 = sum(err_IAST(:, 1:N), 1);
    % % save to text file
    % fid = fopen(sprintf('%s-mode%.0f-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), 'w');
    % fprintf(fid, '#%s\n#Margules: A1 = %0.4g, A2 = %0.4g, C = 0, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, coeff(1:2), err, err_lnpi0, err_psi);
    % fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
    % formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
    % fprintf(fid, formatSpec, [M1_all{k}(:,1), M2_all{k}(:,1), RAST_backcalc, psi_IAST, lnP0_IAST]');
    % fclose(fid); movefile(sprintf('%s-mode%.0f-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), sprintf('data/%s/indsim', zeo_label{i}));
    % fprintf('Saved to text file\n')

    mode=4;
    % mode 4 predictions, with psi
    % [err, RAST_backcalc, x, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(coeffpsi_mo4, isotherm, {M1_all{k}, M2_all{k}}, @Margules_ads, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @Margules_ads_deriv, 'x0', x0, 'options', options);
    % % extract psi and lnpi0 from err_IAST + sum
    % err_psi = sum(err_IAST(:, N+1));
    % err_lnpi0 = sum(err_IAST(:, 1:N), 1);
    % % save to text file
    % fid = fopen(sprintf('%s-mode%.0fPSI-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), 'w');
    % fprintf(fid, '#%s\n#Margules: A1 = %0.4g, A2 = %0.4g, C = %0.4g, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, coeffpsi_mo4, err, err_lnpi0, err_psi);
    % fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
    % formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
    % fprintf(fid, formatSpec, [M1_all{k}(:,1), M2_all{k}(:,1), RAST_backcalc, psi_IAST, lnP0_IAST]');
    % fclose(fid); movefile(sprintf('%s-mode%.0fPSI-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), sprintf('data/%s/indsim', zeo_label{i}));
    % fprintf('Saved to text file\n')

    % mode 4 predictions, without psi
    [err, RAST_backcalc, x, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(coeff_mo4, isotherm, {M1_all{k}, M2_all{k}}, @Margules, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @(x,y)0, 'x0', x0, 'options', options);
    % extract psi and lnpi0 from err_IAST + sum
    err_psi = sum(err_IAST(:, N+1));
    err_lnpi0 = sum(err_IAST(:, 1:N), 1);
    % save to text file
    fid = fopen(sprintf('%s-mode%.0f-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), 'w');
    fprintf(fid, '#%s\n#Margules: A1 = %0.4g, A2 = %0.4g, C = 0, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, coeff_mo4(1:2), err, err_lnpi0, err_psi);
    fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
    formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
    fprintf(fid, formatSpec, [M1_all{k}(:,1), M2_all{k}(:,1), RAST_backcalc, psi_IAST, lnP0_IAST]');
    fclose(fid); movefile(sprintf('%s-mode%.0f-RASTpred-y2-%.1f.txt', zeo_label{i}, mode, y2), sprintf('data/%s/indsim', zeo_label{i}));
    fprintf('Saved to text file\n')
    
    % % IAST predictions
    % x0=[]; options=[];
    % [Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = IAST_solve([M1_all{k}(:, 1), M2_all{k}(:,1)], S, 'options', optimset('Display', 'iter'), 'mode', 1, 'x0', x0);
    % % % extract psi and lnpi0 from err_IAST + sum
    % err_psi = sum(err_IAST(:, N+1));
    % err_lnpi0 = sum(err_IAST(:, 1:N), 1);
    % % get prediction error
    % diff_sq = (relative_error_safe(Q_IAST, [M1_all{k}(:,2), M2_all{k}(:,2)])).^2;
    % err = sum(diff_sq, 'all') / numel(diff_sq);
    % % save to text file
    % fid = fopen(sprintf('%s-IASTpred-y2-%.1f.txt', zeo_label{i}, y2), 'w');
    % fprintf(fid, '#%s\n#IAST, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, err, err_lnpi0, err_psi);
    % fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
    % formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
    % fprintf(fid, formatSpec, [M1_all{k}(:,1), M2_all{k}(:,1), Q_IAST, psi_IAST, lnP0_IAST]');
    % fclose(fid); movefile(sprintf('%s-IASTpred-y2-%.1f.txt', zeo_label{i}, y2), sprintf('data/%s/indsim', zeo_label{i}));
    % fprintf('Saved to text file\n')


    % bulk predictions
    % coeffbulk=[1.6022 0.7947];
    % [err, QRAST_bulk, x, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(coeffbulk, isotherm, {M1_all{k}, M2_all{k}}, @Margules, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @(x,y)0, 'x0', x0, 'options', options);
    % err_psi = sum(err_IAST(:, N+1));
    % err_lnpi0 = sum(err_IAST(:, 1:N), 1);
    % % save to text file
    % fid = fopen(sprintf('%s-BULKpred-y2-%.1f.txt', zeo_label{i}, y2), 'w');
    % fprintf(fid, '#%s\n#Margules: A1 = %0.4g, A2 = %0.4g, C = 0, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, coeffbulk, err, err_lnpi0, err_psi);
    % fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
    % formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
    % fprintf(fid, formatSpec, [M1_all{k}(:,1), M2_all{k}(:,1), QRAST_bulk, psi_IAST, lnP0_IAST]');
    % fclose(fid); movefile(sprintf('%s-BULKpred-y2-%.1f.txt', zeo_label{i}, y2), sprintf('data/%s/indsim', zeo_label{i}));
    % fprintf('Saved to text file\n')

    fprintf(sprintf('DONE WITH RAST PREDICTIONS at y2 = %0.1f\n', y2)); pause(2)

end

%% plot

y2=0.9;
k = find(fixed_comp==y2);

data_mo4PSI = importdata(sprintf('data/%s/indsim/%s-mode4PSI-RASTpred-y2-%0.1f.txt', zeo_label{i}, zeo_label{i}, y2));
data_mo4PSI = data_mo4PSI.textdata(4:end, :);
QRAST_mo4PSI = str2double(data_mo4PSI(:,[3,4]));

data_mo4 = importdata(sprintf('data/%s/indsim/%s-mode4-RASTpred-y2-%0.1f.txt', zeo_label{i}, zeo_label{i}, y2));
data_mo4 = data_mo4.textdata(4:end, :);
QRAST_mo4 = str2double(data_mo4(:,[3,4]));

data_mo1PSI = importdata(sprintf('data/%s/indsim/%s-mode1PSI-RASTpred-y2-%0.1f.txt', zeo_label{i}, zeo_label{i}, y2));
data_mo1PSI = data_mo1PSI.textdata(4:end, :);
QRAST_mo1PSI = str2double(data_mo1PSI(:,[3,4]));

data_mo1 = importdata(sprintf('data/%s/indsim/%s-mode1-RASTpred-y2-%0.1f.txt', zeo_label{i}, zeo_label{i}, y2));
data_mo1 = data_mo1.textdata(4:end, :);
QRAST_mo1 = str2double(data_mo1(:,[3,4]));

data_IAST = importdata(sprintf('data/%s/indsim/%s-IASTpred-y2-%0.1f.txt', zeo_label{i}, zeo_label{i}, y2));
data_IAST = data_IAST.textdata(4:end, :);
Q_IAST = str2double(data_IAST(:,[3,4]));

data_bulk = importdata(sprintf('data/%s/indsim/%s-BULKpred-y2-%0.1f.txt', zeo_label{i}, zeo_label{i}, y2));
data_bulk = data_bulk.textdata(4:end, :);
QRAST_bulk = str2double(data_bulk(:,[3,4]));

% RAST pred vs. logP
close all; figure;
run('startup.m');
fid=sprintf('data/%s/indsim/Q-lgf-fitCampAB_y2_%.1f', zeo_label{i}, y2);
% rast_compare = {[M1_all{k}(:,2), M2_all{k}(:,2)], QRAST_mo4, QRAST_mo4PSI};
rast_compare = {[M1_all{k}(:,2), M2_all{k}(:,2)], QRAST_mo1PSI, QRAST_mo4PSI, Q_IAST, QRAST_bulk};
plot_manu({M1_all{k}(:,1)+M2_all{k}(:,1)}, rast_compare, {'fig-QRAST-lgf', fixed_comp(k)}, zeo_label{i});
% plot_manu({M1_all{k}(:,1)+M2_all{k}(:,1)}, rast_compare, {'fig-QRAST-lgf-cconstr', fixed_comp(k)}, zeo_label{i});
% % savefig(gcf, sprintf('%s.fig', fid))

%% RAST selectivity vs. logP
x2_MC = M2_all{k}(:,2)./(M1_all{k}(:,2)+M2_all{k}(:,2));
x2_mo1PSI = QRAST_mo1PSI(:,2)./sum(QRAST_mo1PSI, 2);
x2_mo4 = QRAST_mo4(:,2)./sum(QRAST_mo4, 2);
x2_mo4PSI = QRAST_mo4PSI(:,2)./sum(QRAST_mo4PSI, 2);
x2_IAST = Q_IAST(:,2)./sum(Q_IAST, 2);
x2_bulk = QRAST_bulk(:,2)./sum(QRAST_bulk, 2);

S_MC = (x2_MC./(1-x2_MC))./(y2/(1-y2));
S_mo1PSI = (x2_mo1PSI./(1-x2_mo1PSI))./(y2/(1-y2));
S_mo4 = (x2_mo4./(1-x2_mo4))./(y2/(1-y2));
S_mo4PSI = (x2_mo4PSI./(1-x2_mo4PSI))./(y2/(1-y2));
S_IAST = (x2_IAST./(1-x2_IAST))./(y2/(1-y2));
S_bulk = (x2_bulk./(1-x2_bulk))./(y2/(1-y2));

figure;
run('startup.m');
fid=sprintf('data/%s/indsim/S-lgf-fitCampAB_y2_%.1f', zeo_label{i}, y2);
sel_compare = {S_MC, S_mo1PSI, S_mo4PSI, S_IAST, S_bulk};
plot_manu({M1_all{k}(:,1)+M2_all{k}(:,1)}, sel_compare, {'fig-S-lgf', fixed_comp(k)}, zeo_label{i});
% sel_compare = {S_MC, S_mo4, S_mo4PSI};
% plot_manu({M1_all{k}(:,1)+M2_all{k}(:,1)}, sel_compare, {'fig-S-lgf-cconstr', fixed_comp(k)}, zeo_label{i});
% savefig(gcf, sprintf('%s.fig', fid))

fprintf(sprintf('DONE WITH PLOTTING Q-lgf + S-lgf at y2 = %0.1f\n', y2));

%% REPRO FIGS: gamma vs. lgf

y2=0.9;
k = find(fixed_comp==y2);

% gMarg_psi = Margules_ads(coeffpsi, [x1_all{k}, psiMC_all{k}]);
% gMarg = Margules(coeff, [x1_all{k}, psiMC_all{k}]);
% gbulk=Margules(coeffbulk, [x1_all{k}, psiMC_all{k}]);


% get gammas from Marg mode4 fit, with and without psi
gmo4_psi = Margules_ads(coeffpsi_mo4, [x1_all{k}, psiMC_all{k}]);
gmo4 = Margules(coeff_mo4, [x1_all{k}, psiMC_all{k}]);

% gamma vs. logP
close all; figure
run('startup.m')
fid=sprintf('data/%s/indsim/gamma_lgf_y2_%.1f',  zeo_label{i}, y2);
% gamma_compare = {gMC_all{k}, gMarg_psi, gMarg, gmo4_psi, gmo4};
% plot_manu({totalP_all{k}}, gamma_compare, {'fig-gamma-lgf', y2}, zeo_label{i})
gamma_compare = {gMC_all{k}, gmo4, gmo4_psi};
plot_manu({totalP_all{k}}, gamma_compare, {'fig-gamma-lgf-cconstr', y2}, zeo_label{i})
savefig(gcf, sprintf('%s.fig', fid))

%% plot ind. MC data w/error bars + compare RUN1 data

i=1; zeo_label{i} = 'DDR-1';
load(sprintf('data/%s/indsim.mat', zeo_label{i}));

blue = {'#add8e6', '#00bfff', '#1e90ff', '#4682b4', '#0000cd', '#add8e6', '#00bfff'};
red = {'#f08080', '#cd5c5c', '#ff6347', '#ffa500', '#8b0000', '#f08080', '#cd5c5c'};

% RUN1 data
close all; figure;
run('startup.m')
idx = 3;
y2 = 0.5;
start_data = 2;
k = find(fixed_comp==y2);
etoh_water_run1 = importdata(sprintf('data/indsim/%s/indsim_CampB_RUN%d_avg.txt', zeo_label{i}, 1));
M={[etoh_water_run1.data(start_data:end,1)-etoh_water_run1.data(start_data:end,end), etoh_water_run1.data(start_data:end,3)], ...
    [etoh_water_run1.data(start_data:end,end), etoh_water_run1.data(start_data:end,4)]};
% conf_binary_run1 = [etoh_water_run1.data(start_data:end,6), etoh_water_run1.data(start_data:end,7)];
% conf95_run1 = {conf_binary_run1(:,1), conf_binary_run1(:,2)};
semilogx(M{1}(:,1)+M{2}(:,1), M{1}(:,2), 'o', 'Color', blue{idx+2}, 'MarkerFaceColor', blue{idx+2}); hold on
semilogx(M{1}(:,1)+M{2}(:,1), M{2}(:,2), 'o', 'Color', red{idx+2}, 'MarkerFaceColor', red{idx+2});
% open(sprintf('data/%s/Q-lgf-fitCampAB_y2_0.5.fig', zeo_label{i}));

% add new data w/error bars. conf95_all = water, etoh, binary1, binary2
fid=sprintf('data/%s/indsim/Q-lgf-err-wprev-y2_%.1f',  zeo_label{i}, y2);
errorbar(totalP_all{k}, M1_all{k}(:,2), conf95_all{k}{3}, 'o', 'Color', blue{idx}, 'MarkerFaceColor', blue{idx}, 'CapSize', 10, 'LineWidth', 1); hold on
errorbar(totalP_all{k}, M2_all{k}(:,2), conf95_all{k}{4}, 'o', 'Color', red{idx}, 'MarkerFaceColor', red{idx}, 'CapSize', 10, 'LineWidth', 1);
xlabel('{\it f}_{total} (Pa)'); ylabel('Q (mol/kg)');
set(gca,'XScale','log');
xlim([1e2 1e6]); ylim([0 5]);
xticks([1e2 1e3 1e4 1e5 1e6]);
yticks([0 1 2 3 4 5 6]);
XAxis.MinorTick = 'on';
YAxis.MinorTick = 'on';
YAxis.MinorTickValues = 0:0.5:6;
% savefig(gcf, sprintf('%s.fig', fid))

%% REPRO FIGS: lgQ vs. lgf

% load(sprintf('data/%s/indsim_campAB_coeffmo4.mat', zeo_label{i}));

% y(1:2) MC unary data, etoh rs, water bo, empty
% y(3:4) unary piece-wise fits, solid lines, start at minlnP
% y(5:6) MC binary data, same as above, filled markers
% y(7:8) MC total binary, filled black diamonds + pw fit solid line
totalM = {[M{1}(:,1)+M{2}(:,1), M{2}(:,2)+M{2}(:,2)]};
[totalM_isotherm, totalM_minlnP, totalM_maxlnP, totalM_ads_pot] = fit_piecewise_polynomial(totalM);
xb=linspace(totalM_minlnP(1), totalM_maxlnP(1), 100)';
xp=[linspace(pure_minlnP(1), pure_maxlnP(1), 100)', linspace(pure_minlnP(2), pure_maxlnP(2), 100)'];

% % y(9) total binary pw fit, add lowP, black dotted line
% totalM_withlowP = {[M_withlowP{k}{1}(:,1)+M_withlowP{k}{2}(:,1), M_withlowP{k}{1}(:,2)+M_withlowP{k}{2}(:,2)]};
% [M_withlowP_isotherm, M_withlowP_minlnP, M_withlowP_maxlnP, M_withlowP_adspot] = fit_piecewise_polynomial(totalM_withlowP);
% x_withlowP = linspace(M_withlowP_minlnP, M_withlowP_maxlnP, 100)';
%
% % y(10) total binary dual-Langmuir fit, add lowP, dashed
% [totalM_dulang_isotherm, totalM_dulang_minlnP, totalM_dulang_maxlnP, totalM_dulang_ads_pot, totalM_dulang_inv_ads_pot] = fit_Langmuir(totalM_withlowP, 2);
% % totalM_x_dulang = linspace(M_withlowP_minlnP, totalM_dulang_maxlnP, 100)';

% % sans lowP water
% S1_sanslow=S{1}(6:end,:);
% [S_sanslow_isotherm, S_sanslow_minlnP, S_sanslow_maxlnP, S_sanslow_adspot] = fit_piecewise_polynomial({S1_sanslow});
% x_sanslow = linspace(S_sanslow_minlnP(1), S_sanslow_maxlnP(1), 100)';

% % y(11:12) unary water/etoh dual-Lang fit, blue/red dashed line
% [dulang_isotherm, dulang_minlnP, dulang_maxlnP, dulang_ads_pot, dulang_inv_ads_pot] = fit_Langmuir({S1_sanslow, S{2}}, 2);
% x_dulang = linspace(S_sanslow_minlnP(1), dulang_maxlnP(1), 100)';
% etoh_x_dulang = log(y2.*exp(x_withlowP))';
% etoh_dulang_isotherm=dulang_isotherm{2};

close all; figure;
run('startup.m');
fid=sprintf('data/%s/indsim/lgQ_lgf_y2_%.1f', zeo_label{i}, y2);
%STOPPED EDITING HERE
iso_compare = {S1_sanslow(:,2), S{2}(:,2), ...
    S_sanslow_isotherm{1}(x_sanslow), pure_isotherm{2}(xp(:,2)), ...
    M1_all{k}(:,2),  M2_all{k}(:,2), ...
    totalM{1}(:,2), totalM_isotherm{1}(xb), ...
    M_withlowP_isotherm{1}(x_withlowP), ...
    totalM_dulang_isotherm{1}(x_withlowP), ...
    dulang_isotherm{1}(x_dulang), ...
    etoh_dulang_isotherm(etoh_x_dulang)};
plot_manu({S1_sanslow(:,1), S{2}(:,1), ...
    exp(x_sanslow), exp(xp(:,2)), ...
    M1_all{k}(:,1), M2_all{k}(:,1), ...
    y2.*totalM{1}(:,1), y2.*exp(xb), ...
    y2.*exp(x_withlowP), y2.*exp(x_withlowP), ...
    exp(x_dulang),  exp(etoh_x_dulang)}, ...
    iso_compare, {'fig-lgQ-lgf', fixed_comp(k)}, zeo_label{i});
savefig(gcf, sprintf('%s.fig', fid))

% fig: logPSI-logf
clear water_adspot etoh_adspot total_adspot total_adspot2 ...
    water_adspot2 etoh_adspot2  totalM_ads_dulang water_ads_dulang ...
    etoh_ads_dulang

% for j = 1:length(lnpi0_all{k}(:,1))
%     water_adspot2(j)=pure_ads_pot{1}(lnpi0_all{k}(j,1));
%     etoh_adspot2(j)=pure_ads_pot{2}(lnpi0_all{k}(j,2));
% end
for j = 1:length(xp(:,1))
    % water_adspot(j)=pure_ads_pot{1}(xp(j,1));
    etoh_adspot(j)=pure_ads_pot{2}(xp(j,2));
    % etoh_adspot2(j)=pure_ads_pot{2}(xp(j,1));
end
for j = 1:length(xb)
    total_adspot2(j)=totalM_ads_pot{1}(xb(j));
end
for j = 1:length(x_sanslow(:,1))
    water_adspot_sanslow(j)=S_sanslow_adspot{1}(x_sanslow(j,1));
end
for j = 1:length(x_withlowP)
    totalM_ads_withlowP(j)=M_withlowP_adspot{1}(x_withlowP(j));
end
totalM_ads_dulang=totalM_dulang_ads_pot{1}(x_withlowP);
water_ads_dulang=dulang_ads_pot{1}(x_dulang);
etoh_ads_dulang=dulang_ads_pot{2}(etoh_x_dulang);


close all; figure;
run('startup.m');
fid=sprintf('data/%s/lgPSI_lgf_y2_%.1f', zeo_label{i}, y2);
adspot_compare = {water_adspot_sanslow, etoh_adspot, total_adspot2, totalM_ads_withlowP, totalM_ads_dulang, water_ads_dulang, etoh_ads_dulang};
plot_manu({exp(x_sanslow), exp(xp(:,2)), y2.*exp(xb), y2.*exp(x_withlowP), exp(x_dulang), exp(etoh_x_dulang)}, adspot_compare, {'fig-lgPSI-lgf', fixed_comp(k)}, zeo_label{i});
savefig(gcf, sprintf('%s.fig', fid))

%% fig: Q-lgf-CampA

clear all;
i=1; zeo_label{i} = 'DDR-1';
load(sprintf('data/%s/indsim_campAB_coeffmo4.mat', zeo_label{i}));
% M_campA = {flip(M_campA{1}, 1), flip(M_campA{2}, 1)};

x0=[]; 
options=optimset('Display','iter');
[Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = IAST_solve([M_campA{1}(:, 1), M_campA{2}(:,1)], [], 'isotherm', isotherm, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'options', options, 'mode', 1);
% extract psi and lnpi0 from err_IAST + sum
err_psi = sum(err_IAST(:, N+1));
err_lnpi0 = sum(err_IAST(:, 1:N), 1);
% get prediction error
diff_sq = (relative_error_safe(Q_IAST, [M_campA{1}(:,2), M_campA{2}(:,2)])).^2;
err = sum(diff_sq, 'all') / numel(diff_sq);
% save to text file
fid = fopen(sprintf('%s-IASTpred-CampA.txt', zeo_label{i}), 'w');
fprintf(fid, '#%s\n#IAST, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, err, err_lnpi0, err_psi);
fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
fprintf(fid, formatSpec, [M_campA{1}(:,1), M_campA{2}(:,1), Q_IAST, psi_IAST, lnP0_IAST]');
fclose(fid); movefile(sprintf('%s-IASTpred-CampA.txt', zeo_label{i}), sprintf('data/%s/indsim', zeo_label{i}));
fprintf('Saved to text file\n')

% bulk predictions
coeffbulk=[1.6022 0.7947];
[err_bulk, QRAST_bulk, x, err_IAST1, lnP0_IAST1, psi_IAST1] = RAST_func_IAST_solve(coeffbulk, isotherm, M_campA, @Margules, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @(x,y)0, 'x0', x0, 'options', options);
% extract psi and lnpi0 from err_IAST + sum
err_psi = sum(err_IAST1(:, N+1));
err_lnpi0 = sum(err_IAST1(:, 1:N), 1);
% save to text file
fid = fopen(sprintf('%s-BULKpred-CampA.txt', zeo_label{i}), 'w');
fprintf(fid, '#%s\n#Margules: A1 = %0.4g, A2 = %0.4g, C = 0, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', zeo_label{i}, coeffbulk, err_bulk, err_lnpi0, err_psi);
fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
fprintf(fid, formatSpec, [M_campA{1}(:,1), M_campA{2}(:,1), QRAST_bulk, psi_IAST1, lnP0_IAST1]');
fclose(fid); movefile(sprintf('%s-BULKpred-CampA.txt', zeo_label{i}), sprintf('data/%s/indsim', zeo_label{i}));
fprintf('Saved to text file\n')

% mode 4 with psi
x0=[];
[Q_predicted, x_init, err, lnP0, psi, Q_IAST2, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M_campA, S, 'mode', 4, 'C_ub', 1e2, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', @Margules_ads, 'x0', x0, 'N_EoS_param', 3, 'EoS_deriv', @Margules_ads_deriv, 'options', options);
x0=[Q_predicted(:,2)./sum(Q_predicted,2), lnP0];
[err_mo4, QRAST_mo4PSI, x, err_IAST2, lnP0_IAST2, psi_IAST2] = RAST_func_IAST_solve(coeffpsi_mo4, isotherm, M_campA, @Margules_ads, 'mode', 1, 'x_lb', [0], 'x_ub', [1], 'tol', 1e-6, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @Margules_ads_deriv, 'x0', x0, 'options', options);
err_psi = sum(err_IAST2(:, N+1));
err_lnpi0 = sum(err_IAST2(:, 1:N), 1);
% save to text file
mode=4;
fid = fopen(sprintf('%s-mode%.0fPSI-RASTpred-CampA.txt', zeo_label{i}, mode), 'w');
fprintf(fid, '#%s\n#Margules: A1 = %0.4g, A2 = %0.4g, C = %0.4g, err_Qp = %0.4g, err_lnp0_H20 = %0.4g, err_lnp0_EtOH = %0.4g, err_psi = %0.4g\n', ...
    zeo_label{i}, coeffpsi_mo4, err_mo4, err_lnpi0, err_psi);
fprintf(fid, '#f_H20[Pa] f_EtOH[Pa] Qp_H20[mol/kg] Qp_EtOH[mol/kg] psi_H20[mol/kg] psi_EtOH[mol/kg] lnp0_H20 lnp0_EtOH\n');
formatSpec = '%0.2e %0.2e %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g\n';
fprintf(fid, formatSpec, [M_campA{1}(:,1), M_campA{2}(:,1), QRAST_mo4PSI, psi_IAST2, lnP0_IAST2]');
fclose(fid); movefile(sprintf('%s-mode%.0fPSI-RASTpred-CampA.txt', zeo_label{i}, mode), sprintf('data/%s/indsim', zeo_label{i}));
fprintf('Saved to text file\n')

%%
close all; figure;
run('startup.m');
fid=sprintf('data/%s/indsim/QCampA_lgf', zeo_label{i});
Q_compare = {M_campA{1}(:,2), M_campA{2}(:,2), QRAST_mo1PSI, QRAST_mo4PSI, QRAST_bulk, Q_IAST};
f_compare = {M_campA{2}(:,1)};
plot_manu(f_compare, Q_compare, {'fig-QCampA-lgf', []}, zeo_label{i});
% savefig(gcf, sprintf('%s.fig', fid))

% S-lgf-CampA

% M_campA = {flip(M_campA{1}, 1), flip(M_campA{2}, 1)};

y2 = M_campA{2}(:,1)./(M_campA{1}(:,1)+M_campA{2}(:,1));
x2_MC = M_campA{2}(:,2)./(M_campA{1}(:,2)+M_campA{2}(:,2));
x2_mo1 = QRAST_mo1PSI(:,2)./sum(QRAST_mo1PSI, 2);
x2_mo4 = QRAST_mo4PSI(:,2)./sum(QRAST_mo4PSI, 2);
x2_IAST = Q_IAST(:,2)./sum(Q_IAST, 2);
x2_bulk = QRAST_bulk(:,2)./sum(QRAST_bulk, 2);

S_MC = (x2_MC./(1-x2_MC))./(y2./(1-y2));
S_mo1 = (x2_mo1./(1-x2_mo1))./(y2./(1-y2));
S_mo4 = (x2_mo4./(1-x2_mo4))./(y2./(1-y2));
S_IAST = (x2_IAST./(1-x2_IAST))./(y2./(1-y2));
S_bulk = (x2_bulk./(1-x2_bulk))./(y2./(1-y2));

figure;
run('startup.m');
fid=sprintf('data/%s/indsim/S-lgf-CampA', zeo_label{i});
sel_compare = {S_MC, S_mo1, S_mo4, S_IAST, S_bulk};
plot_manu({M_campA{2}(:,1)}, sel_compare, {'fig-S-lgf', []}, zeo_label{i});
xlabel('{\it f}_{EtOH} (Pa)'); ylabel('{S}_{EtOH}');
xlim([1e2 1e5]); ylim([1 1e3]);
xticks([1e2 1e3 1e4 1e5]);
savefig(gcf, sprintf('%s.fig', fid))