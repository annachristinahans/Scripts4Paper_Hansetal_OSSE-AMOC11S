%% Figures Appendix

%% - Comparing surface options for Mooring method

% --- loading ---
nmid = 'VIK';

% dT/dp
file_name = ['wbr2ebr2_vst1eof0sur0_idr1_' nmid];
vst1brs2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd1 = 'Gradient, m-sur #0';
file_name = ['wbr2ebr2_vst1eof0sur1_idr1_' nmid];
vst1brs2sur1 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd2 = 'Gradient, m-sur #1';
file_name = ['wbr2ebr2_vst1eof0sur3_idr1_' nmid];
vst1brs2sur3 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd3 = 'Gradient, m-sur #3';

% T''
file_name = ['wbr2ebr2_vst2eof0sur0_idr1_' nmid];
vst2brs2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd4 = 'Linear, m-sur #0';
file_name = ['wbr2ebr2_vst2eof0sur1_idr1_' nmid];
vst2brs2sur1 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd5 = 'Linear, m-sur #1';
file_name = ['wbr2ebr2_vst2eof0sur2_idr1_' nmid];
vst2brs2sur2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd6 = 'Linear, m-sur #2';

% EOF
file_name = ['wbr2ebr2_vst3eof2sur0_idr1_' nmid];
vst3brs2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
file_name = ['wbr2ebr2_vst3eof1sur0_idr1_' nmid];
vst3brs2_eb = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd7 = 'EOF, m-sur #0';
file_name = ['wbr2ebr2_vst3eof2sur1_idr1_' nmid];
vst3brs2sur1 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
file_name = ['wbr2ebr2_vst3eof1sur1_idr1_' nmid];
vst3brs2sur1_eb = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_density_' file_name '.mat']);
lgd8 = 'EOF, m-sur #1';

time_mly = vst1brs2.time_mly;
dep_grid = vst1brs2.dep_grid;

% --- computing density errors ---
vst1brs2.wb_rho_diff = vst1brs2.wb_rho_grid_orig - vst1brs2.wb_rho_grid;
vst1brs2.eb_rho_diff = vst1brs2.eb_rho_grid_orig - vst1brs2.eb_rho_grid;
vst1brs2sur1.wb_rho_diff = vst1brs2sur1.wb_rho_grid_orig - vst1brs2sur1.wb_rho_grid;
vst1brs2sur1.eb_rho_diff = vst1brs2sur1.eb_rho_grid_orig - vst1brs2sur1.eb_rho_grid;
vst1brs2sur3.wb_rho_diff = vst1brs2sur3.wb_rho_grid_orig - vst1brs2sur3.wb_rho_grid;
vst1brs2sur3.eb_rho_diff = vst1brs2sur3.eb_rho_grid_orig - vst1brs2sur3.eb_rho_grid;

vst2brs2.wb_rho_diff = vst2brs2.wb_rho_grid_orig - vst2brs2.wb_rho_grid;
vst2brs2.eb_rho_diff = vst2brs2.eb_rho_grid_orig - vst2brs2.eb_rho_grid;
vst2brs2sur1.wb_rho_diff = vst2brs2sur1.wb_rho_grid_orig - vst2brs2sur1.wb_rho_grid;
vst2brs2sur1.eb_rho_diff = vst2brs2sur1.eb_rho_grid_orig - vst2brs2sur1.eb_rho_grid;
vst2brs2sur2.wb_rho_diff = vst2brs2sur2.wb_rho_grid_orig - vst2brs2sur2.wb_rho_grid;
vst2brs2sur2.eb_rho_diff = vst2brs2sur2.eb_rho_grid_orig - vst2brs2sur2.eb_rho_grid;

vst3brs2.wb_rho_diff = vst3brs2.wb_rho_grid_orig - vst3brs2.wb_rho_grid;
vst3brs2.eb_rho_diff = vst3brs2_eb.eb_rho_grid_orig - vst3brs2_eb.eb_rho_grid;
vst3brs2sur1.wb_rho_diff = vst3brs2sur1.wb_rho_grid_orig - vst3brs2sur1.wb_rho_grid;
vst3brs2sur1.eb_rho_diff = vst3brs2sur1_eb.eb_rho_grid_orig - vst3brs2sur1_eb.eb_rho_grid;

% --- plotting ---
figure()
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

col1 = [224 0 101]./255;    % red, define colors for interpolation methods (vst)
col2 = [15 160 235]./255;   % blue
col3 = [224 201 0]./255;    % yellow

nexttile % WB4 sur comp
h(7) = plot(mean(abs(vst3brs2.wb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col3);
hold on
h(8) = plot(mean(abs(vst3brs2sur1.wb_rho_diff),2),dep_grid,'-.o','LineWidth',2.5,'Color',col3,'MarkerSize',10);

h(4) = plot(mean(abs(vst2brs2.wb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col2);
h(5) = plot(mean(abs(vst2brs2sur1.wb_rho_diff),2),dep_grid,'-.o','LineWidth',2.5,'Color',col2,'MarkerSize',10);
h(6) = plot(mean(abs(vst2brs2sur2.wb_rho_diff),2),dep_grid,'-.*','LineWidth',2.5,'Color',col2);

h(1) = plot(mean(abs(vst1brs2.wb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col1);
h(2) = plot(mean(abs(vst1brs2sur1.wb_rho_diff),2),dep_grid,'-.o','LineWidth',2.5,'Color',col1,'MarkerSize',10);
h(3) = plot(mean(abs(vst1brs2sur3.wb_rho_diff),2),dep_grid,'-..','LineWidth',2.5,'Color',col1);

axis ij
ylim([0 150])
ylabel('Depth [m]')
title('WB4 mooring','FontWeight','normal')
xlabel('Mean absolute density error [kg m^-^3]')
legend(h,lgd1,lgd2,lgd3,lgd4,lgd5,lgd6,lgd7,lgd8,'Location','southeast')
legend('boxoff')
set(gca,'FontSize',20)
text(0, 1.01, '(a)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile % EB1 sur comp
plot(mean(abs(vst3brs2.eb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col3)
hold on
plot(mean(abs(vst3brs2sur1.eb_rho_diff),2),dep_grid,'-.o','LineWidth',2.5,'Color',col3,'MarkerSize',10)

plot(mean(abs(vst2brs2.eb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col2)
plot(mean(abs(vst2brs2sur1.eb_rho_diff),2),dep_grid,'-.o','LineWidth',2.5,'Color',col2,'MarkerSize',10)
plot(mean(abs(vst2brs2sur2.eb_rho_diff),2),dep_grid,'-.*','LineWidth',2.5,'Color',col2)

plot(mean(abs(vst1brs2.eb_rho_diff),2),dep_grid,'-','LineWidth',2.5,'Color',col1)
plot(mean(abs(vst1brs2sur1.eb_rho_diff),2),dep_grid,'-.o','LineWidth',2.5,'Color',col1,'MarkerSize',10)
plot(mean(abs(vst1brs2sur3.eb_rho_diff),2),dep_grid,'-..','LineWidth',2.5,'Color',col1)

axis ij
ylim([0 150])
title('EB1 mooring','FontWeight','normal')
xlabel('Mean absolute density error [kg m^-^3]')
set(gca,'FontSize',20)
text(0, 1.01, '(b)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

set(gcf,'Position',[1 76 1200 580])

%% - Overview of climatological gradients for McCarthy 2015 interpolation (m-vst_option = 1)

% --- loading ---
nmid = 'VIK';
load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_dTdp_dSdp.mat']);

load('romaO.mat');
idx = round(linspace(1, size(romaO,1), 12));
cmap12 = romaO(idx, :);

% --- plotting ---
% as function of pressure
figure()
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile
for i = 1:12
    plot(dTdp_wb(:,i),pres_ddp(1,:),'Color',cmap12(i,:),'LineWidth',2)
    hold on
end
axis ij
set(gca,'FontSize',15)
ylim([0 1200])
xlim([-0.3 0.01])
hold on
yline(200,'k')
xline(-0.1,'k')
title('WB4 mooring','FontWeight','normal')
ylabel('Pressure [dbar]')
xlabel('dT/dp [°C dbar^-^1]')
text(0, 1.01, '(a)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
for i = 1:12
    plot(dTdp_eb(:,i),pres_ddp(2,:),'Color',cmap12(i,:),'LineWidth',2)
    hold on
end
axis ij
set(gca,'FontSize',15)
ylim([0 1200])
xlim([-0.3 0.01])
hold on
yline(200,'k')
xline(-0.1,'k')
title('EB1 mooring','FontWeight','normal')
xlabel('dT/dp [°C dbar^-^1]')
text(0, 1.01, '(b)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
for i = 1:12
    plot(dSdp_wb(:,i),pres_ddp(1,:),'Color',cmap12(i,:),'LineWidth',2)
    hold on
end
axis ij
set(gca,'FontSize',15)
ylim([0 1200])
xlim([-0.03 0.22])
hold on
yline(200,'k')
xline(0.02,'k')
xline(-0.02,'k')
ylabel('Pressure [dbar]')
xlabel('dS/dp [psu dbar^-^1]')
text(0, 1.01, '(c)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
for i = 1:12
    plot(dSdp_eb(:,i),pres_ddp(2,:),'Color',cmap12(i,:),'LineWidth',2)
    hold on
end
axis ij
set(gca,'FontSize',15)
ylim([0 1200])
xlim([-0.03 0.22])
hold on
yline(200,'k')
xline(0.02,'k')
xline(-0.02,'k')
xlabel('dS/dp [psu dbar^-^1]')
text(0, 1.01, '(d)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

colormap(cmap12)
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Ticks = 0+1/24:1/12:1;
cb.TickLabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
cb.TickLength = 0;
set(gcf,'Position',[300 100 851 608])

%% BP Reconstruction

nmid = 'VIK';
ebb3 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc1_ebb3_ref1_vst1_stp1_' nmid '.mat']);
prc4 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc4_ebb0_ref1_vst1_stp1_' nmid '.mat']);
prc5 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc5_ebb0_ref1_vst1_stp1_' nmid '.mat']);
prc6 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc6_ebb0_ref1_vst1_stp1_' nmid '.mat']);
time_mly = prc5.time_mly;

% low-pass filter for WB
Fs = 1/(30*24*3600); % Sampling frequency in Hz (monthly)
Fc_lp = 1/(5*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(2, Fc_lp/(Fs/2), 'low'); % 2nd-order low-pass filter
bpa_wb300_orig_lp = filtfilt(b, a, prc5.bpa_wb300_mly_orig);
bpa_wb500_orig_lp = filtfilt(b, a, prc5.bpa_wb500_mly_orig);
prc4.bpa_wb300_lp = filtfilt(b, a, prc4.bpa_wb300_mly);
prc5.bpa_wb300_lp = filtfilt(b, a, prc5.bpa_wb300_mly);
prc6.bpa_wb300_lp = filtfilt(b, a, prc6.bpa_wb300_mly);
prc4.bpa_wb500_lp = filtfilt(b, a, prc4.bpa_wb500_mly);
prc5.bpa_wb500_lp = filtfilt(b, a, prc5.bpa_wb500_mly);
prc6.bpa_wb500_lp = filtfilt(b, a, prc6.bpa_wb500_mly);

% high-pass filter for EB
Fc_hp = 1/(2*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(1, Fc_hp/(Fs/2), 'high'); % 1st-order high-pass filter
bpa_eb300_orig_hp = filtfilt(b, a, ebb3.bpa_eb300_mly_orig);
bpa_eb500_orig_hp = filtfilt(b, a, ebb3.bpa_eb500_mly_orig);
ebb3.bpa_eb300_hp = filtfilt(b, a, ebb3.bpa_eb300_mly);
ebb3.bpa_eb500_hp = filtfilt(b, a, ebb3.bpa_eb500_mly);

% cutting half the filter length in the beginning and end
bpa_wb300_orig_lp = bpa_wb300_orig_lp(31:end-30);
bpa_wb500_orig_lp = bpa_wb500_orig_lp(31:end-30);
prc4.bpa_wb300_lp = prc4.bpa_wb300_lp(31:end-30);
prc5.bpa_wb300_lp = prc5.bpa_wb300_lp(31:end-30);
prc6.bpa_wb300_lp = prc6.bpa_wb300_lp(31:end-30);
prc4.bpa_wb500_lp = prc4.bpa_wb500_lp(31:end-30);
prc5.bpa_wb500_lp = prc5.bpa_wb500_lp(31:end-30);
prc6.bpa_wb500_lp = prc6.bpa_wb500_lp(31:end-30);
bpa_eb300_orig_hp = bpa_eb300_orig_hp(13:end-12);
bpa_eb500_orig_hp = bpa_eb500_orig_hp(13:end-12);
ebb3.bpa_eb300_hp = ebb3.bpa_eb300_hp(13:end-12);
ebb3.bpa_eb500_hp = ebb3.bpa_eb500_hp(13:end-12);

% plotting
figure()
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile
h(1) = plot(time_mly(31:end-30),bpa_wb300_orig_lp,'k','LineWidth',5);
hold on
h(4) = plot(time_mly(31:end-30),prc6.bpa_wb300_lp,'LineWidth',3);
h(3) = plot(time_mly(31:end-30),prc5.bpa_wb300_lp,'LineWidth',3);
h(2) = plot(time_mly(31:end-30),prc4.bpa_wb300_lp,'LineWidth',3);
ylim([-0.0151 0.011])
xticks(datenum(1985:5:2020,1,1))
datetick('x','keepticks')
xlim([time_mly(31) time_mly(end-30)])
set(gca,'FontSize',20)
title('WBb5 (low-pass filtered)')
legend(h,'Model truth','b-prc #4','b-prc #5a','b-prc #5b','Location','southwest')
legend box off
text(0, 1.01, '(a)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
plot(time_mly(13:end-12),bpa_eb300_orig_hp,'k','LineWidth',5)
hold on
plot(time_mly(13:end-12),ebb3.bpa_eb300_hp,'Color',[121 205 205]./255,'LineWidth',3) % alternative colour: 0 191 255 or 0 178 238
ylim([-0.018 0.024])
xticks(datenum(2014:2:2024,1,1))
datetick('x','keepticks')
xlim([datenum(2013,1,1) time_mly(end-12)])
title('EBb2 (high-pass filtered)')
legend('Model truth','b-ebb #3')
legend box off
set(gca,'FontSize',20)
text(0, 1.01, '(b)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
plot(time_mly(31:end-30),bpa_wb500_orig_lp,'k','LineWidth',5)
hold on
plot(time_mly(31:end-30),prc6.bpa_wb500_lp,'LineWidth',3)
plot(time_mly(31:end-30),prc5.bpa_wb500_lp,'LineWidth',3)
plot(time_mly(31:end-30),prc4.bpa_wb500_lp,'LineWidth',3)
ylim([-0.0151 0.011])
xticks(datenum(1985:5:2020,1,1))
datetick('x','keepticks')
xlim([time_mly(31) time_mly(end-30)])
set(gca,'FontSize',20)
title('WBb6 (low-pass filtered)')
text(0, 1.01, '(c)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
plot(time_mly(13:end-12),bpa_eb500_orig_hp,'k','LineWidth',5)
hold on
plot(time_mly(13:end-12),ebb3.bpa_eb500_hp,'Color',[121 205 205]./255,'LineWidth',3)
ylim([-0.018 0.024])
xticks(datenum(2014:2:2024,1,1))
datetick('x','keepticks')
xlim([datenum(2013,1,1) time_mly(end-12)])
title('EBb3 (high-pass filtered)')
set(gca,'FontSize',20)
text(0, 1.01, '(d)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

ylabel(t,'Bottom pressure [dbar]','FontSize',25)
set(gcf,'Position',[2 50 1400 720])

% --- related statistics (for text) ---
% bottom pressure trends for prc 0 ebr 0 (-> original data)
bp_orig = load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/bpr_pos22_prc0_ebb0_ref1_vst1_stp1_VIK.mat');
bp_orig.mdl_eb300 = fitlm(bp_orig.time_mly(bp_orig.time_mly>=datenum(1994,1,1)), bp_orig.bpa_eb300_mly(bp_orig.time_mly>=datenum(1994,1,1)));
bp_orig.mdl_eb500 = fitlm(bp_orig.time_mly(bp_orig.time_mly>=datenum(1994,1,1)), bp_orig.bpa_eb500_mly(bp_orig.time_mly>=datenum(1994,1,1)));
bp_orig.mdl_wb300 = fitlm(bp_orig.time_mly(bp_orig.time_mly>=datenum(1994,1,1)), bp_orig.bpa_wb300_mly(bp_orig.time_mly>=datenum(1994,1,1)));
bp_orig.mdl_wb500 = fitlm(bp_orig.time_mly(bp_orig.time_mly>=datenum(1994,1,1)), bp_orig.bpa_wb500_mly(bp_orig.time_mly>=datenum(1994,1,1)));

bp_orig.slope_per_year_eb300 = bp_orig.mdl_eb300.Coefficients.Estimate(2) * 365.25; % in dbar/year
bp_orig.slope_per_year_eb500 = bp_orig.mdl_eb500.Coefficients.Estimate(2) * 365.25;
bp_orig.slope_per_year_wb300 = bp_orig.mdl_wb300.Coefficients.Estimate(2) * 365.25;
bp_orig.slope_per_year_wb500 = bp_orig.mdl_wb500.Coefficients.Estimate(2) * 365.25;
bp_orig.slope_delta_300 = abs(bp_orig.slope_per_year_eb300 - bp_orig.slope_per_year_wb300);
bp_orig.slope_delta_500 = abs(bp_orig.slope_per_year_eb500 - bp_orig.slope_per_year_wb500);
