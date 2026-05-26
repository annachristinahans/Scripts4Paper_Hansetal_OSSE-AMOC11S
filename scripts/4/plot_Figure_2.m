%% Figure 2 - Overview AMOC in VIKING20X

% --- loading and preparing data ---
vik_v = load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/v_1100S_VIK.mat');
vik_vg = load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/vg_1100S_VIK.mat');
vik_ek = load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/TEk_1100S_VIK.mat');

% averaging to monthly values
[yy,mm,~] = datevec(vik_v.time);
time_mly = unique(datenum(yy,mm,15));
for i = 1:length(time_mly)
    vik_v.psi_mly(:,i) = mean(vik_v.psi(:,datenum(yy,mm,15)==time_mly(i)),2);
    vik_vg.psig_mly(:,i) = mean(vik_vg.psig(:,datenum(yy,mm,15)==time_mly(i)),2);
    vik_vg.Tg_mly(:,i) = mean(vik_vg.Tg_basin(:,datenum(yy,mm,15)==time_mly(i)),2);
    vik_ek.T_Ek_mly(:,i) = mean(vik_ek.T_Ek(datenum(yy,mm,15)==time_mly(i)));
end
[vik_v.AMOC_mly,vik_v.zmax_mly] = max(vik_v.psi_mly,[],1,'omitnan');
[vik_vg.AMOCg_mly,vik_vg.zmax_mly] = max(vik_vg.psig_mly,[],1,'omitnan');

vik_vg.mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), vik_vg.AMOCg_mly(time_mly>=datenum(1994,1,1)).*1e-6);
vik_vg.slope_per_decade = vik_vg.mdl.Coefficients.Estimate(2) * 365.25 * 10;

vik_ek.mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), vik_ek.T_Ek_mly(time_mly>=datenum(1994,1,1)).*1e-6);
vik_ek.slope_per_decade = vik_ek.mdl.Coefficients.Estimate(2) * 365.25 * 10;

% --- plotting ---
figure() % AMOCg and Ekman
tiledlayout(3,3,'TileSpacing','compact','Padding','compact')

nexttile(1,[2,3])
colororder([0, 0, 0; 130, 175, 40]./255)
yyaxis left
plot(time_mly,vik_vg.AMOCg_mly.*1e-6,'k','LineWidth',2)
hold on
plot(time_mly(time_mly>=datenum(1994,1,1)),predict(vik_vg.mdl,time_mly(time_mly>=datenum(1994,1,1))'),'w-','LineWidth',5)
plot(time_mly(time_mly>=datenum(1994,1,1)),predict(vik_vg.mdl,time_mly(time_mly>=datenum(1994,1,1))'),'k--','LineWidth',3)
ylim([13 33])
ylabel('AMOC_g [Sv]')
ax = gca;
ax.YTick = 15:5:30;
yyaxis right
plot(time_mly,vik_ek.T_Ek_mly.*1e-6,'LineWidth',2,'Color',[130, 175, 40]./255)
hold on
plot(time_mly(time_mly>=datenum(1994,1,1)),predict(vik_ek.mdl,time_mly(time_mly>=datenum(1994,1,1))'),'w-','LineWidth',5)
plot(time_mly(time_mly>=datenum(1994,1,1)),predict(vik_ek.mdl,time_mly(time_mly>=datenum(1994,1,1))'),'--','Color',[130, 175, 40]./255,'LineWidth',3)
ylim([-15 5])
ylabel('Ekman [Sv]')
ax = gca;
ax.XTick = datenum(1980,1,15):365.25*5:datenum(2023,1,1);
ax.XTickLabel = {'1980','1985','1990','1995','2000','2005','2010','2015','2020'};
xlim([time_mly(1) time_mly(end)])
set(gca,'FontSize',15)
text(0.005, 0.99, '(a)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','top')

nexttile(7,[1,1]) % missing ageostrophic processes
plot(time_mly,(vik_v.AMOC_mly - vik_vg.AMOCg_mly - vik_ek.T_Ek_mly).*1e-6,'LineWidth',2,'Color',[.5 .5 .5])
hold on
yline(0,'k')
datetick
ylabel({'AMOC - AMOC_g','- Ekman [Sv]'})
xlim([time_mly(1) time_mly(end)])
ylim([-2 1.5])
ax = gca;
ax.XTick = datenum(1980,1,15):365.25*10:datenum(2023,1,1);
ax.XTickLabel = {'1980','1990', '2000', '2010', '2020'};
set(gca,'FontSize',15)
text(0.015, 0.99, '(b)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','top')

nexttile(8,[1,1]) % impact fixed overturning depth
dep_fixed = 1100;
plot(time_mly,(vik_vg.AMOCg_mly - vik_vg.psig_mly(vik_vg.dep_grid == dep_fixed,:)).*1e-6,'LineWidth',2,'Color',[.5 .5 .5])
hold on
yline(0,'k')
datetick
ylabel({'AMOC_g', '- \Psi_g(1100 m) [Sv]'})
xlim([time_mly(1) time_mly(end)])
%ylim([-2 1.5])
ax = gca;
ax.XTick = datenum(1980,1,15):365.25*10:datenum(2023,1,1);
ax.XTickLabel = {'1980','1990', '2000', '2010', '2020'};
set(gca,'FontSize',15)
text(0.015, 0.99, '(c)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','top')

nexttile(9,[1,1]) % psig
plot(mean(vik_vg.psig,2).*1e-6,vik_vg.dep_grid,'k','LineWidth',2)
hold on
yline(mean(vik_vg.dep_grid(vik_vg.zmax_mly)),'k')
axis ij
xlabel('\Psi_g [Sv]')
ylabel('Depth [m]')
ax = gca;
ax.YTick = [0 1000 1100 2000 3000 4000 5000 6000];
ax.YTickLabel = {'0', '', '1100', '', '3000', '', '', '6000'};
set(gca,'FontSize',15)
text(0.015, 0.99, '(d)', 'Units','normalized', 'FontSize',15, 'HorizontalAlignment','left', 'VerticalAlignment','top')

set(gcf,'Position',[100 100 1000 600])

%% --- related statistics (for text) ---
% trend since 1994 for all components
mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), vik_v.AMOC_mly(time_mly>=datenum(1994,1,1)).*1e-6);
slope_per_decade_AMOC = mdl.Coefficients.Estimate(2) * 365.25 * 10;

mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), vik_vg.AMOCg_mly(time_mly>=datenum(1994,1,1)).*1e-6);
slope_per_decade_AMOCg = mdl.Coefficients.Estimate(2) * 365.25 * 10;

mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), vik_ek.T_Ek_mly(time_mly>=datenum(1994,1,1)).*1e-6);
slope_per_decade_Ek = mdl.Coefficients.Estimate(2) * 365.25 * 10;

ageo = vik_v.AMOC_mly - vik_vg.AMOCg_mly - vik_ek.T_Ek_mly;
mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), ageo(time_mly>=datenum(1994,1,1)).*1e-6);
slope_per_decade_ageo = mdl.Coefficients.Estimate(2) * 365.25 * 10;

% statistics for neglected ageostrophic component
offset_mean1 = mean(vik_v.AMOC_mly - vik_vg.AMOCg_mly - vik_ek.T_Ek_mly).*1e-6;

% statistics for FixedDepth AMOCg
offset_mean2 = mean(vik_vg.psig_mly(vik_vg.dep_grid == dep_fixed,:) - vik_vg.AMOCg_mly).*1e-6;
[r,p] = corrcoef(vik_vg.psig_mly(vik_vg.dep_grid == dep_fixed,:),vik_vg.AMOCg_mly);
r_all = r(1,2);
%rmse_all = sqrt(sum((vik_vg.psig_mly(vik_vg.dep_grid == dep_fixed,:) - vik_vg.AMOCg_mly).^2) ./length(vik_vg.AMOCg_mly)).*1e-6;
%rmse_meanrem = sqrt(sum((vik_vg.psig_mly(vik_vg.dep_grid == dep_fixed,:) - mean(vik_vg.psig_mly(vik_vg.dep_grid == dep_fixed,:)) - vik_vg.AMOCg_mly + mean(vik_vg.AMOCg_mly)).^2) ./length(vik_vg.AMOCg_mly)).*1e-6;
