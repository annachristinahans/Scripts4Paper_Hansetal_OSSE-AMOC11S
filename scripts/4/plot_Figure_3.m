%% Figure 3 - Comparing AMOCg for BPR method

% --- loading and preparing data ---
nmid = 'VIK';
file_name = ['pos22_prc0_ebb0_ref1_vst1_stp1_' nmid];
test1 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_' file_name '.mat']);
lgd1 = 'Idealised approach (BPR Case 1)';
file_name = ['pos22_prc1_ebb1_ref1_vst1_stp1_' nmid];
test2 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_' file_name '.mat']);
lgd2 = 'Current approach (BPR Case 9)';
time_mly = test1.time_mly;

test1.AMOCg = test1.AMOCg - mean(test1.AMOCg);
test2.AMOCg = test2.AMOCg - mean(test2.AMOCg);

% loading ideal transport for comparison
ideal = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/vg_1100S_' nmid '.mat']);
% calculate their monthly means
[yy,mm,~] = datevec(ideal.time);
for i = 1:length(time_mly)
    ideal.Tg_monthly(:,i) = mean(ideal.Tg_basin(:,datenum(yy,mm,15)==time_mly(i)),2);
    ideal.psig_monthly(:,i) = mean(ideal.psig(:,datenum(yy,mm,15)==time_mly(i)),2);
end
[ideal.AMOCg,~] = max(ideal.psig_monthly.*1e-6,[],1,'omitnan');
ideal.AMOCg = ideal.AMOCg - mean(ideal.AMOCg);

% --- filtering and calcuation of statistics for AMOCg

% filtering for high and low frequent variability seperately
Fs = 1/(30*24*3600); % Sampling frequency in Hz (monthly)
Fc_lp = 1/(5*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(2, Fc_lp/(Fs/2), 'low'); % 1st-order low-pass filter
test1.AMOCg_lp = filtfilt(b, a, test1.AMOCg);
test2.AMOCg_lp = filtfilt(b, a, test2.AMOCg);
ideal.AMOCg_lp = filtfilt(b, a, ideal.AMOCg);
Fc_hp = 1/(2*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(1, Fc_hp/(Fs/2), 'high'); % 1st-order low-pass filter
test1.AMOCg_hp = filtfilt(b, a, test1.AMOCg);
test2.AMOCg_hp = filtfilt(b, a, test2.AMOCg);
ideal.AMOCg_hp = filtfilt(b, a, ideal.AMOCg);

test1.AMOCg_lp = test1.AMOCg_lp(31:end-30);
test2.AMOCg_lp = test2.AMOCg_lp(31:end-30);
ideal.AMOCg_lp = ideal.AMOCg_lp(31:end-30);
test1.AMOCg_hp = test1.AMOCg_hp(13:end-12);
test2.AMOCg_hp = test2.AMOCg_hp(13:end-12);
ideal.AMOCg_hp = ideal.AMOCg_hp(13:end-12);

% rmse and R
[stats_test1, stats_ideal, ~] = f_stats_calc(time_mly,test1.AMOCg,ideal.AMOCg);
[stats_test2, ~, ~] = f_stats_calc(time_mly,test2.AMOCg,ideal.AMOCg);

% x/y limits for the following plots
min_AMOCg = min([test1.AMOCg test2.AMOCg ideal.AMOCg]);
max_AMOCg = max([test1.AMOCg test2.AMOCg ideal.AMOCg]);
min_AMOCg_lp = min([test1.AMOCg_lp test2.AMOCg_lp ideal.AMOCg_lp]);
max_AMOCg_lp = max([test1.AMOCg_lp test2.AMOCg_lp ideal.AMOCg_lp]);
min_AMOCg_hp = min([test1.AMOCg_hp test2.AMOCg_hp ideal.AMOCg_hp]);
max_AMOCg_hp = max([test1.AMOCg_hp test2.AMOCg_hp ideal.AMOCg_hp]);

% monthly climatology
[~,mm,~] = datevec(time_mly);
for i = 1:12
    clim_ideal(i) = mean(ideal.AMOCg(mm == i));
    clim_test1(i) = mean(test1.AMOCg(mm == i));
    clim_test2(i) = mean(test2.AMOCg(mm == i));
end

% --- plotting ---
figure()
rgb1 = [0.686, 0.496, 0.781]; % for idealised approach
rgb2 = [0 96 143]./255; % for current approach
t = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

nexttile(1,[2,1])
plot(ideal.AMOCg,test1.AMOCg,'*','MarkerSize',10,'Color',rgb1)
hold on
plot(ideal.AMOCg,test2.AMOCg,'*','MarkerSize',10,'Color',rgb2)
plot([min_AMOCg max_AMOCg],[min_AMOCg max_AMOCg],'k')
xlim([min_AMOCg max_AMOCg])
ylim([min_AMOCg max_AMOCg])
legend(['rmse = ' num2str(stats_test1.rmse_fs,'%.1f') ' Sv, R = ' num2str(stats_test1.r_fs,'%.2f')],['rmse = ' num2str(stats_test2.rmse_fs,'%.1f') ' Sv, R = ' num2str(stats_test2.r_fs,'%.2f')], 'Location','northwest')
legend('boxoff')
title('Full time series','FontWeight','normal')
ylabel('Reconstructed AMOC_g'' [Sv]','FontSize',18)
set(gca,'FontSize',18)
text(0, 1.01, '(a)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(2,[2,1]) % highpass filtered
plot(ideal.AMOCg_hp,test1.AMOCg_hp,'*','MarkerSize',10,'Color',rgb1)
hold on
plot(ideal.AMOCg_hp,test2.AMOCg_hp,'*','MarkerSize',10,'Color',rgb2)
plot([min_AMOCg_hp max_AMOCg_hp],[min_AMOCg_hp max_AMOCg_hp],'k')
xlim([min_AMOCg_hp max_AMOCg_hp])
ylim([min_AMOCg_hp max_AMOCg_hp])
yticks(-5:5:5)
legend(['rmse = ' num2str(stats_test1.rmse_hp,'%.1f') ' Sv, R = ' num2str(stats_test1.r_hp,'%.2f')],['rmse = ' num2str(stats_test2.rmse_hp,'%.1f') ' Sv, R = ' num2str(stats_test2.r_hp,'%.2f')], 'Location','northwest')
legend('boxoff')
title('2-year high-pass filtered','FontWeight','normal')
xlabel('Model-true AMOC_g'' [Sv]','FontSize',18)
set(gca,'FontSize',18)
text(0, 1.01, '(b)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(3,[2,1]) % lowpass filtered
plot(ideal.AMOCg_lp,test1.AMOCg_lp,'*','MarkerSize',10,'Color',rgb1)
hold on
plot(ideal.AMOCg_lp,test2.AMOCg_lp,'*','MarkerSize',10,'Color',rgb2)
plot([min_AMOCg_lp max_AMOCg_lp],[min_AMOCg_lp max_AMOCg_lp],'k')
xlim([min_AMOCg_lp max_AMOCg_lp])
ylim([min_AMOCg_lp max_AMOCg_lp])
yticks(-2:1)
legend(['rmse = ' num2str(stats_test1.rmse_lp,'%.1f') ' Sv, R = ' num2str(stats_test1.r_lp,'%.2f')],['rmse = ' num2str(stats_test2.rmse_lp,'%.1f') ' Sv, R = ' num2str(stats_test2.r_lp,'%.2f')], 'Location','northwest')
legend('boxoff')
text(0.22,0.15,lgd1,'Units','normalized','Color',rgb1,'FontSize',18)
text(0.22,0.08,lgd2,'Units','normalized','Color',rgb2,'FontSize',18)
title('5-year low-pass filtered','FontWeight','normal')
set(gca,'FontSize',18)
text(0, 1.01, '(c)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(7,[1,2])
plot(time_mly(13:end-12), test1.AMOCg_hp,'Color',rgb1,'LineWidth',2)
hold on
plot(time_mly(13:end-12), test2.AMOCg_hp,'Color',rgb2,'LineWidth',2)
plot(time_mly(13:end-12), ideal.AMOCg_hp,'k','LineWidth',2)
datetick
xlim([datenum(1982,06,1) datenum(2021,06,31)])
title('2-year high-pass filtered','FontWeight','normal')
ylabel('AMOC_g'' [Sv]')
xticks(datenum(2014:2:2024,1,1))
datetick('x','keepticks')
xlim([datenum(2013,1,1) time_mly(end-12)])
set(gca,'FontSize',18)
text(0, 1.01, '(d)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(9)
plot(1:13,[clim_ideal clim_ideal(1)],'k','LineWidth',3)
hold on
plot(1:13,[clim_test1 clim_test1(1)],'LineWidth',3,'Color',rgb1)
plot(1:13,[clim_test2 clim_test2(1)],'LineWidth',3,'Color',rgb2)
ax = gca;
ax.XTick = 1:13;
ax.XTickLabel = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A','S','O','N','D','J'};
xlim([1 13])
title('Monthly climatology','Fontweight','normal')
set(gca,'FontSize',18)
text(0, 1.01, '(e)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(10,[1,3])
plot(time_mly(31:end-30), test1.AMOCg_lp,'Color',rgb1,'LineWidth',2)
hold on
plot(time_mly(31:end-30), test2.AMOCg_lp,'Color',rgb2,'LineWidth',2)
plot(time_mly(31:end-30), ideal.AMOCg_lp,'k','LineWidth',2)
datetick
xlim([datenum(1982,06,1) datenum(2021,06,31)])
title('5-year low-pass filtered','FontWeight','normal')
ylabel('AMOC_g'' [Sv]')
ax = gca;
ax.XTick = datenum(1985,1,1):365.25*5:datenum(2023,1,1);
ax.XTickLabel = {'1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'};
set(gca,'FontSize',18)
text(0, 1.01, '(f)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

set(gcf,'Position',[2 250 1250 720])

%% --- related statistics (for text) ---

% monthly climatologies for observational period 2013-2018 (5 yrs)
[yy,mm,~] = datevec(time_mly);
mm(yy<2013) = 0;
mm(yy>=2018) = 0;
for i = 1:12
    clim_ideal_obs(i) = mean(ideal.AMOCg(mm == i));
    clim_test2_obs(i) = mean(test2.AMOCg(mm == i));
end

disp('Peak-to-peak amplitudes of seasonal cycles for 2013-2018')
disp('Model truth:')
disp(max(clim_ideal_obs) - min(clim_ideal_obs))
disp('Current approach:')
disp(max(clim_test2_obs) - min(clim_test2_obs))

% resulting over/underestimations compared to Herrford et al. (2021) 12.2 Sv
disp('Model underestimates Herrford by factor')
disp(12.2/(max(clim_ideal_obs) - min(clim_ideal_obs))) % value for section 2.1
disp(((12.2 - (max(clim_ideal_obs) - min(clim_ideal_obs)))/12.2)*100) % same in %
disp('Subsampled observational strategy still underestimates Herrford by factor')
disp(12.2/(max(clim_test2_obs) - min(clim_test2_obs))) % value for Discussion
disp(((12.2 - (max(clim_test2_obs) - min(clim_test2_obs)))/12.2)*100) % same in %
