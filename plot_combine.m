% combine plots

%% six-panel
clc; close all; clear all;
plot_zeo = {'MFI-1', 'DDR-1'};

% manufig='S-lgf-fitCampAB';
manufig='lgPSI_lgf';
% manufig='gamma_lgf_withlowP';
fid=sprintf('data/manu/fig_rev1/%s', manufig);

fid1=sprintf('data/%s/indsim/%s_y2_%.1f.fig',  plot_zeo{1}, manufig, 0.1);
fid2=sprintf('data/%s/indsim/%s_y2_%.1f.fig',  plot_zeo{1}, manufig, 0.5);
fid3=sprintf('data/%s/indsim/%s_y2_%.1f.fig',  plot_zeo{1}, manufig, 0.9);
fid4=sprintf('data/%s/indsim/%s_y2_%.1f.fig',  plot_zeo{2}, manufig, 0.1);
fid5=sprintf('data/%s/indsim/%s_y2_%.1f.fig',  plot_zeo{2}, manufig, 0.5);
fid6=sprintf('data/%s/indsim/%s_y2_%.1f.fig',  plot_zeo{2}, manufig, 0.9);

% fid1=sprintf('data/%s/indsim/%s_f%.0f.fig',  plot_zeo{1}, manufig, 900);
% fid2=sprintf('data/%s/indsim/%s_f%.0f.fig',  plot_zeo{1}, manufig, 30000);
% fid3=sprintf('data/%s/indsim/%s_f%.0f.fig',  plot_zeo{1}, manufig, 300000);
% fid4=sprintf('data/%s/indsim/%s_f%.0f.fig',  plot_zeo{2}, manufig, 900);
% fid5=sprintf('data/%s/indsim/%s_f%.0f.fig',  plot_zeo{2}, manufig, 30000);
% fid6=sprintf('data/%s/indsim/%s_f%.0f.fig',  plot_zeo{2}, manufig, 300000);

figure;
tcl=tiledlayout(2,3, 'TileSpacing', 'none');
h1 = openfig(fid1); ax1 = h1.CurrentAxes; set(ax1,'Xticklabel',[], 'Xlabel', [], 'Ylabel', []);
ax1.Parent=tcl; ax1.Layout.Tile=1;
h2 = openfig(fid2); ax2 = h2.CurrentAxes; set(ax2,'Xticklabel',[], 'Yticklabel',[], 'Xlabel', [], 'Ylabel', []);
ax2.Parent=tcl; ax2.Layout.Tile=2;
h3 = openfig(fid3); ax3 = h3.CurrentAxes; set(ax3,'Xticklabel',[], 'Yticklabel',[], 'Xlabel', [], 'Ylabel', []);
ax3.Parent=tcl; ax3.Layout.Tile=3;
h4 = openfig(fid4, 'reuse'); ax4 = h4.CurrentAxes; set(ax4,'Xlabel', [], 'Ylabel', []);
ax4.Parent=tcl; ax4.Layout.Tile=4;
h5 = openfig(fid5, 'reuse'); ax5 = h5.CurrentAxes; set(ax5,'Yticklabel',[], 'Xlabel', [], 'Ylabel', []);
ax5.Parent=tcl; ax5.Layout.Tile=5;
h6 = openfig(fid6, 'reuse'); ax6 = h6.CurrentAxes; set(ax6,'Yticklabel',[], 'Xlabel', [], 'Ylabel', []);
ax6.Parent=tcl; ax6.Layout.Tile=6;

xlabel(tcl, '{\it f} (Pa)', 'FontSize', 22);
% xlabel(tcl, '{\it f}_{total} (Pa)', 'FontSize', 22);
% xlabel(tcl, '{\it x}_{EtOH}', 'FontSize', 22);
% ylabel(tcl, 'S_{ads}', 'FontSize', 22);
% ylabel(tcl, 'Q (mol/kg)', 'FontSize', 22);
% ylabel(tcl, '\gamma', 'FontSize', 22);
ylabel(tcl, '\Psi', 'FontSize', 22);
% savefig(gcf, sprintf('%s.fig', fid))

%% three-panel

clc; close all; clear all;
plot_zeo = 'BEA-1';

% manufig='gamma_lgf';
% manufig='Q-lgf-fitCampAB';
manufig='S-lgf-fitCampAB';
fid=sprintf('data/manu/fig_rev1/%s_%s', manufig, plot_zeo);

fid1=sprintf('data/%s/indsim/%s_y2_%.1f.fig', plot_zeo, manufig, 0.1);
fid2=sprintf('data/%s/indsim/%s_y2_%.1f.fig', plot_zeo, manufig, 0.5);
fid3=sprintf('data/%s/indsim/%s_y2_%.1f.fig', plot_zeo, manufig, 0.9);

figure;
tcl=tiledlayout(1,3, 'TileSpacing', 'none');
h1 = openfig(fid1); ax1 = h1.CurrentAxes; set(ax1, 'Xlabel', [], 'Ylabel', []); %, 'Ylim', [0 12]);
ax1.Parent=tcl; ax1.Layout.Tile=1;
h2 = openfig(fid2); ax2 = h2.CurrentAxes; set(ax2, 'Yticklabel',[], 'Xlabel', [], 'Ylabel', []); %, 'Ylim', [0 12]);
ax2.Parent=tcl; ax2.Layout.Tile=2;
h3 = openfig(fid3); ax3 = h3.CurrentAxes; set(ax3, 'Yticklabel',[], 'Xlabel', [], 'Ylabel', []); %, 'Ylim', [0 12]);
ax3.Parent=tcl; ax3.Layout.Tile=3;

xlabel(tcl, '{\it f}_{total} (Pa)', 'FontSize', 22);
% ylabel(tcl, '\gamma', 'FontSize', 22);
% ylabel(tcl, 'Q (mol/kg)', 'FontSize', 22);
ylabel(tcl, 'S_{ads}', 'FontSize', 22);
savefig(gcf, sprintf('%s.fig', fid))

%% four-panel

clc; close all; clear all;
plot_zeo = {'DDR-1', 'MFI-1', 'BEA-1', 'LTA-2'};

% manufig='QCampA_lgf';
manufig='S-lgf-CampA';
fid=sprintf('data/manu/fig_rev1/%s', manufig);

fid1=sprintf('data/%s/indsim/%s.fig', plot_zeo{1}, manufig);
fid2=sprintf('data/%s/indsim/%s.fig', plot_zeo{2}, manufig);
fid3=sprintf('data/%s/indsim/%s.fig', plot_zeo{3}, manufig);
fid4=sprintf('data/%s/indsim/%s.fig', plot_zeo{4}, manufig);

figure;
tcl=tiledlayout(2,2, 'TileSpacing', 'none');
h1 = openfig(fid1); ax1 = h1.CurrentAxes; set(ax1,'Xticklabel',[], 'Xlabel', [], 'Xlabel', [], 'Ylabel', [], 'XLim', [50 1e5]);
ax1.Parent=tcl; ax1.Layout.Tile=1;
h2 = openfig(fid2); ax2 = h2.CurrentAxes; set(ax2,'Xticklabel',[], 'Yticklabel',[], 'Xlabel', [], 'Ylabel', [],'XLim', [50 1e5]);
ax2.Parent=tcl; ax2.Layout.Tile=2;
h3 = openfig(fid3); ax3 = h3.CurrentAxes; set(ax3,'Xlabel', [], 'Ylabel', [], 'XLim', [50 1e5]);
ax3.Parent=tcl; ax3.Layout.Tile=3;
h4 = openfig(fid4); ax4 = h4.CurrentAxes; set(ax4, 'Yticklabel',[], 'Xlabel', [], 'Ylabel', [], 'XLim', [50 1e5]);
ax4.Parent=tcl; ax4.Layout.Tile=4;

xlabel(tcl, '{\it f}_{EtOH} (Pa)', 'FontSize', 22);
% ylabel(tcl, 'Q (mol/kg)', 'FontSize', 22);
% xlabel(tcl, '{\it f}_{total} (Pa)', 'FontSize', 22);
ylabel(tcl, 'S_{EtOH}', 'FontSize', 22);
savefig(gcf, sprintf('%s.fig', fid))

%% three-panel
clc; close all; clear all;
plot_combo = {'c_pld', 'gamma_lcd', 'gamma_pld'};

manufig='gamma-C-pld';
fid=sprintf('data/manu/%s', manufig);

fid1=sprintf('data/manu/%s', plot_combo{1});
fid2=sprintf('data/manu/%s', plot_combo{2});
fid3=sprintf('data/manu/%s', plot_combo{3});

figure;
tcl=tiledlayout(1,3, 'TileSpacing', 'tight');
h1 = openfig(fid1); ax1 = h1.CurrentAxes; set(ax1);
ax1.Parent=tcl; ax1.Layout.Tile=1;
h2 = openfig(fid2); ax2 = h2.CurrentAxes; set(ax2);
ax2.Parent=tcl; ax2.Layout.Tile=2;
h3 = openfig(fid3); ax3 = h3.CurrentAxes; set(ax3);
ax3.Parent=tcl; ax3.Layout.Tile=3;

savefig(gcf, sprintf('%s.fig', fid))

%% three-panel type 2

clc; close all; clear all;
plot_combo = {'c_pld', 'gamma_lcd', 'gamma_pld'};

manufig='gamma-C-pld';
fid=sprintf('data/manu/%s', manufig);

fid1=sprintf('data/manu/%s', plot_combo{1});
fid2=sprintf('data/manu/%s', plot_combo{2});
fid3=sprintf('data/manu/%s', plot_combo{3});

figure;
tcl = tiledlayout(2,2, 'TileSpacing', 'tight');
h1 = openfig(fid1); ax1 = h1.CurrentAxes; set(ax1);
ax1.Parent=tcl;  ax1.Layout.Tile=[2 1]; %ax1.Layout.Tile=1;
h2 = openfig(fid2); ax2 = h2.CurrentAxes; set(ax2);
ax2.Parent=tcl;  ax2.Layout.Tile=3;%ax2.Layout.Tile=2;
h3 = openfig(fid3); ax3 = h3.CurrentAxes; set(ax3);
ax3.Parent=tcl;  ax3.Layout.Tile=4; %ax3.Layout.Tile=3;

% savefig(gcf, sprintf('%s.fig', fid))

% nexttile
% contour(X,Y,Z)


