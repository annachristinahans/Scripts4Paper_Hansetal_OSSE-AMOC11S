%% Figure 5 - Comparing density reconstruction for Mooring method

% --- loading data ---
nmid = 'VIK';

% dT/dp sur0
file_name = ['wbr1ebr1_vst1eof0sur0_idr1_' nmid];
vst1brs1 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd1 = 'Initial, Gradient';
file_name = ['wbr2ebr2_vst1eof0sur0_idr1_' nmid];
vst1brs2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd2 = 'Current, Gradient';
file_name = ['wbr3ebr3_vst1eof0sur0_idr1_' nmid];
vst1brs3 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd3 = 'Enhanced, Gradient';

% T'' sur0
file_name = ['wbr1ebr1_vst2eof0sur0_idr1_' nmid];
vst2brs1 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd4 = 'Initial, Linear';
file_name = ['wbr2ebr2_vst2eof0sur0_idr1_' nmid];
vst2brs2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd5 = 'Current, Linear';

% EOF sur0
file_name = ['wbr1ebr1_vst3eof1sur0_idr1_' nmid];
vst3brs1 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd6 = 'Initial, EOF';
file_name = ['wbr2ebr2_vst3eof2sur0_idr1_' nmid];
vst3brs2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
file_name = ['wbr2ebr2_vst3eof1sur0_idr1_' nmid];
vst3brs2_eb = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd7 = 'Current, EOF';

time_mly = vst1brs1.time_mly;
dep_grid = vst1brs1.dep_grid;

% --- computing density errors ---
vst1brs1.wb_rho_diff = vst1brs1.wb_rho_grid_orig - vst1brs1.wb_rho_grid;
vst1brs1.eb_rho_diff = vst1brs1.eb_rho_grid_orig - vst1brs1.eb_rho_grid;
vst1brs2.wb_rho_diff = vst1brs2.wb_rho_grid_orig - vst1brs2.wb_rho_grid;
vst1brs2.eb_rho_diff = vst1brs2.eb_rho_grid_orig - vst1brs2.eb_rho_grid;
vst1brs3.wb_rho_diff = vst1brs3.wb_rho_grid_orig - vst1brs3.wb_rho_grid;
vst1brs3.eb_rho_diff = vst1brs3.eb_rho_grid_orig - vst1brs3.eb_rho_grid;

vst2brs1.wb_rho_diff = vst2brs1.wb_rho_grid_orig - vst2brs1.wb_rho_grid;
vst2brs1.eb_rho_diff = vst2brs1.eb_rho_grid_orig - vst2brs1.eb_rho_grid;
vst2brs2.wb_rho_diff = vst2brs2.wb_rho_grid_orig - vst2brs2.wb_rho_grid;
vst2brs2.eb_rho_diff = vst2brs2.eb_rho_grid_orig - vst2brs2.eb_rho_grid;

vst3brs1.wb_rho_diff = vst3brs1.wb_rho_grid_orig - vst3brs1.wb_rho_grid;
vst3brs1.eb_rho_diff = vst3brs1.eb_rho_grid_orig - vst3brs1.eb_rho_grid;
vst3brs2.wb_rho_diff = vst3brs2.wb_rho_grid_orig - vst3brs2.wb_rho_grid;
vst3brs2.eb_rho_diff = vst3brs2_eb.eb_rho_grid_orig - vst3brs2_eb.eb_rho_grid;

% --- defining MC positions for plot ---
% define instrument positions
dep_VIK = 1e+3 .* [0.0030    0.0095    0.0164    0.0239    0.0322    0.0415    0.0519    0.0639    0.0776    0.0936    0.1123    0.1343   0.1603    0.1911    0.2276    0.2709    0.3220    0.3821    0.4524    0.5340    0.6279    0.7347    0.8551    0.9892    1.1369    1.2977    1.4709    1.6555    1.8504    2.0544    2.2665    2.4854    2.7101    2.9398    3.1736    3.4108    3.6507    3.8929    4.1370    4.3827    4.6295    4.8773    5.1259    5.3752    5.6250    5.8751];
mc_wb4_i = dep_VIK([1 20 21 29]); % WB4 initial: surface, 500m, 650m, 1900m
mc_wb4_c = dep_VIK([1 10 14 16 18 20 21 29]); % WB4 current: surface, 100m, 200m, 300m, 400m, 500m, 650m, 1900m
mc_wb4_e = dep_VIK([1 7 10 14 16 18 20 21 23 25 29]); % WB4 enhanced: surface, 50m, 100m, 200m, 300m, 400m, 500m, 650m, 850m, 1100m, 1900m

mc_eb1_i = dep_VIK([1 17 20 22 24 26]); % EB1 initial: surface, 300m, 500m, 700m, 950m, 1200m
mc_eb1_i(end) = 1234; % adjust as bottom cell is actually thinner
mc_eb1_c = dep_VIK([1 13 17 20 22 24 26]); % EB1 current: surface, 150m, 300m, 500m, 700m, 950m, 1200m
mc_eb1_c(end) = 1234;
mc_eb1_e = dep_VIK([1 4 10 14 16 18 20 22 24 26]); % EB1 enhanced: surface, 20m, 100m, 200m, 300m, 400m, 500m, 700m, 950m, 1200m
mc_eb1_e(end) = 1234;

mc_size = 30;

% --- plotting ---
figure()
tiledlayout(2,3,'TileSpacing','compact','Padding','compact')

col1 = [224 0 101]./255;    % red, define colors for interpolation methods (vst)
col2 = [15 160 235]./255;   % blue
col3 = [224 201 0]./255;    % yellow

nexttile % WB4 Setup comp
plot(mean(abs(vst1brs1.wb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col1)
hold on
plot(mean(abs(vst1brs2.wb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col1)
plot(mean(abs(vst1brs3.wb_rho_diff),2),dep_grid,':','LineWidth',2.5,'Color',col1)
plot(repmat(0.138,1,length(mc_wb4_i)),mc_wb4_i,'k.','Markersize',mc_size)
plot(repmat(0.144,1,length(mc_wb4_c)),mc_wb4_c,'k.','Markersize',mc_size)
plot(repmat(0.15,1,length(mc_wb4_e)),mc_wb4_e,'k.','Markersize',mc_size)
axis ij
ylim([min(dep_grid) max(dep_grid)])
ylabel({'WB4 mooring', 'Depth [m]'})
legend(lgd1,lgd2,lgd3,'Location','southeast')
legend('boxoff')
set(gca,'FontSize',20)
text(0, 1.01, '(a)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile % WB4 Initial
h1(3) = plot(mean(abs(vst3brs1.wb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col3);
hold on
h1(2) = plot(mean(abs(vst2brs1.wb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col2);
h1(1) = plot(mean(abs(vst1brs1.wb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col1);
plot(repmat(0.25,1,length(mc_wb4_i)),mc_wb4_i,'k.','Markersize',mc_size)
axis ij
ylim([min(dep_grid) max(dep_grid)])
xlabel('   ')
legend(h1,lgd1,lgd4,lgd6,'Location','southeast')
legend('boxoff')
set(gca,'FontSize',20)
text(0, 1.01, '(b)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile % WB4 Current
h2(3) = plot(mean(abs(vst3brs2.wb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col3);
hold on
h2(2) = plot(mean(abs(vst2brs2.wb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col2);
h2(1) = plot(mean(abs(vst1brs2.wb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col1);
plot(repmat(0.115,1,length(mc_wb4_c)),mc_wb4_c,'k.','Markersize',mc_size)
axis ij
ylim([min(dep_grid) max(dep_grid)])
legend(h2,lgd2,lgd5,lgd7,'Location','southeast')
legend('boxoff')
set(gca,'FontSize',20)
text(0, 1.01, '(c)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile % EB1 Setup comp
plot(mean(abs(vst3brs1.eb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col3)
hold on
plot(mean(abs(vst3brs2.eb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col3)
plot(mean(abs(vst1brs3.eb_rho_diff),2),dep_grid,':','LineWidth',2.5,'Color',col1)
plot(repmat(0.37,1,length(mc_eb1_i)),mc_eb1_i,'k.','Markersize',mc_size)
plot(repmat(0.385,1,length(mc_eb1_c)),mc_eb1_c,'k.','Markersize',mc_size)
plot(repmat(0.4,1,length(mc_eb1_e)),mc_eb1_e,'k.','Markersize',mc_size)
axis ij
ylim([min(dep_grid) max(dep_grid)])
legend(lgd6,lgd7,lgd3,'Location','southeast')
legend('boxoff')
ylabel({'EB1 mooring', 'Depth [m]'})
set(gca,'FontSize',20)
text(0, 1.01, '(d)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile % EB1 Initial
plot(mean(abs(vst3brs1.eb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col3)
hold on
plot(mean(abs(vst2brs1.eb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col2)
plot(mean(abs(vst1brs1.eb_rho_diff),2),dep_grid,'--','LineWidth',2.5,'Color',col1)
plot(repmat(1,1,length(mc_eb1_i)),mc_eb1_i,'k.','Markersize',mc_size)
axis ij
ylim([min(dep_grid) max(dep_grid)])
xlabel('Mean absolute density error [kg m^-^3]')
set(gca,'FontSize',20)
text(0, 1.01, '(e)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile % EB1 Current
plot(mean(abs(vst3brs2.eb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col3)
hold on
plot(mean(abs(vst2brs2.eb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col2)
plot(mean(abs(vst1brs2.eb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col1)
plot(repmat(0.8,1,length(mc_eb1_c)),mc_eb1_c,'k.','Markersize',mc_size)
axis ij
ylim([min(dep_grid) max(dep_grid)])
set(gca,'FontSize',20)
text(0, 1.01, '(f)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

set(gcf,'Position',[1 76 1440 720])
