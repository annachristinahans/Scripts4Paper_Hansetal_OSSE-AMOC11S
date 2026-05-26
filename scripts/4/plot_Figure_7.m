%% Figure 7 - Summary

nmid = 'VIK';
col1 = [0 96 143]./255; % for current approach
col2 = [175 195 0]./255; % for improvement

% --- loading and preparing ---
% --- LEFT SIDE - Time series ---

file_name = ['pos22_prc1_ebb1_ref1_vst1_stp1_' nmid ''];
bpr_l9 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_' file_name '.mat']);
time_mly = bpr_l9.time_mly;

file_name = ['pos22_prc0_ebb2_ref1_vst1_stp1_' nmid '']; % eventually change prc to 1
bpr_l21 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_' file_name '.mat']);

file_name = ['wbr3ebr2_vst1eof0sur1_idr2_' nmid ''];
moor_l8 = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/moor_' file_name '.mat']);

% loading ideal transport for comparison
ideal = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/vg_1100S_' nmid '.mat']);
% calculate their monthly means
[yy,mm,~] = datevec(ideal.time);
for i = 1:length(time_mly)
    ideal.psig_mly(:,i) = mean(ideal.psig(:,datenum(yy,mm,15)==time_mly(i)),2);
    ideal.Tg_wbeb_mly(:,i) = mean(ideal.Tg_wb(:,datenum(yy,mm,15)==time_mly(i)) + ideal.Tg_shelf_eb(:,datenum(yy,mm,15)==time_mly(i)),2);
end
[ideal.AMOCg_mly,~] = max(ideal.psig_mly.*1e-6,[],1,'omitnan');

% adding boundary transport to array estimate for mooring
ideal.AMOCg_wbeb_mly = sum(ideal.Tg_wbeb_mly(ideal.dep_grid<=max(moor_l8.dep_grid),:).*10,1).*1e-6;
moor_l8.AMOCg = moor_l8.AMOCg + ideal.AMOCg_wbeb_mly;

% compare anomalies
ideal.AMOCg_mly = ideal.AMOCg_mly - mean(ideal.AMOCg_mly);
ideal.AMOCg_wbeb_mly = ideal.AMOCg_wbeb_mly - mean(ideal.AMOCg_wbeb_mly);
moor_l8.AMOCg = moor_l8.AMOCg - mean(moor_l8.AMOCg);

% filtering
Fs = 1/(30*24*3600); % Sampling frequency in Hz (monthly)
Fc_lp = 1/(5*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(2, Fc_lp/(Fs/2), 'low'); % 1st-order low-pass filter
moor_l8.AMOCg_lp = filtfilt(b, a, moor_l8.AMOCg);
ideal.AMOCg_wbeb_lp = filtfilt(b, a, ideal.AMOCg_wbeb_mly);
ideal.AMOCg_lp = filtfilt(b, a, ideal.AMOCg_mly);

Fc_hp = 1/(2*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(1, Fc_hp/(Fs/2), 'high'); % 1st-order high-pass filter
bpr_l9.AMOCg_hp = filtfilt(b, a, bpr_l9.AMOCg);
bpr_l21.AMOCg_hp = filtfilt(b, a, bpr_l21.AMOCg);
ideal.AMOCg_hp = filtfilt(b, a, ideal.AMOCg_mly);

% cutting half the filter length in the beginning and end
moor_l8.AMOCg_lp = moor_l8.AMOCg_lp(31:end-30);
ideal.AMOCg_wbeb_lp = ideal.AMOCg_wbeb_lp(31:end-30);
ideal.AMOCg_lp = ideal.AMOCg_lp(31:end-30);

bpr_l9.AMOCg_hp = bpr_l9.AMOCg_hp(13:end-12);
bpr_l21.AMOCg_hp = bpr_l21.AMOCg_hp(13:end-12);
ideal.AMOCg_hp = ideal.AMOCg_hp(13:end-12);

% --- RIGHT SIDE - Sketches ---
topo = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/v_topo_tracos_' nmid '.mat']);
y_max = 2000;

% define stretch of xaxis by cutting in the middle
xmin = -36;
x1 = -34.5;
x2 = 12;
xmax = 13.5;
offset_east = 45;
mc_size = 30;

% define instrument positions
dep_VIK = 1e+3 .* [0.0030    0.0095    0.0164    0.0239    0.0322    0.0415    0.0519    0.0639    0.0776    0.0936    0.1123    0.1343   0.1603    0.1911    0.2276    0.2709    0.3220    0.3821    0.4524    0.5340    0.6279    0.7347    0.8551    0.9892    1.1369    1.2977    1.4709    1.6555    1.8504    2.0544    2.2665    2.4854    2.7101    2.9398    3.1736    3.4108    3.6507    3.8929    4.1370    4.3827    4.6295    4.8773    5.1259    5.3752    5.6250    5.8751];

moor_lat = [-10.27, -10.38, -10.61, -10.94, -10.83]; % WB1, WB2, WB3, WB4, EB1
moor_lon = [-35.86+0.016, -35.68, -35.39, -34.99, 13.00-1/20];
moor_dep = [900, 2400, 3500, 4090, 1250];

mc_wb4 = dep_VIK([1 7 10 14 16 18 20 21 23 25 29]); % WB4ext+: surface, 50m, 100m, 200m, 300m, 400m, 500m, 650m, 850m, 1100m, 1900m
mc_eb1 = dep_VIK([1 13 17 20 22 24 26]); % EB1ext: surface, 150m, 300m, 500m, 700m, 950m, 1200m, in observations: 300 500 700 950 1200
mc_eb1(end) = 1234; % adjust as bottom cell is actually thinner

bpr_lon = [-35.87-0.018, -35.86-0.016, -35.86+0.028, 13.23, 13.19+0.01, 13+0.03-1/20]; % WBb5, WBb6, WBb7, EBb2, EBb3, EBb4
bpr_dep = [272, 441, 1105, 295, 439, 1252-9];
bpr_Tg_dep = [284 440 1179];
ref_level = 1100;


% ----------------
% --- plotting ---
fig1 = figure();
tiledlayout(2,40,'TileSpacing','compact','Padding','compact')

% - Time series
nexttile(1,[1,24]) % hp version
colororder({'k','k'})
yyaxis left
yline(0,'Color',[.5 .5 .5])
hold on
h(1) = plot(time_mly(13:end-12),ideal.AMOCg_hp,'k-','LineWidth',3);
h(2) = plot(time_mly(13:end-12),bpr_l9.AMOCg_hp,'-','Color',col1,'LineWidth',3);
h(3) = plot(time_mly(13:end-12),bpr_l21.AMOCg_hp,'-','Color',col2,'LineWidth',3);
ylim([-12 5.5])
yticks([-5 0 5])
ylabel('Reconstructed AMOCg'' [Sv]')
yyaxis right
bar(time_mly(13:end-12),abs(ideal.AMOCg_hp - bpr_l9.AMOCg_hp),'FaceColor',col1,'FaceAlpha',0.7)
hold on
bar(time_mly(13:end-12),abs(ideal.AMOCg_hp - bpr_l21.AMOCg_hp),'FaceColor',col2,'FaceAlpha',0.7)
ylim([0 12])
yticks([0 2 4])
ylabel('Diff. to model truth [Sv]')

lgd = legend(h,'Model truth','Current approach (BPR Case 9)','Improvement (BPR Case 21)','Location','northoutside','Orientation','horizontal');
title(lgd,'2-year high-pass filtered','FontWeight','normal')
legend boxoff
set(gca,'FontSize',20)
xticks(datenum(2014:2:2024,1,1))
datetick('x','keepticks')
xlim([datenum(2013,1,1) time_mly(end-12)])
text(0, 1.01, '(a)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(41,[1,24]) % lp version (without detrending)
colororder({'k','k'})
yyaxis left
yline(0,'Color',[.5 .5 .5])
hold on
h(1) = plot(time_mly(31:end-30),ideal.AMOCg_lp,'k-','LineWidth',3);
h(2) = plot(time_mly(31:end-30),ideal.AMOCg_wbeb_lp,'-','Color',col1,'LineWidth',3);
h(3) = plot(time_mly(31:end-30),moor_l8.AMOCg_lp,'-','Color',col2,'LineWidth',3);
ylim([-4 3])
yticks([-3 0 3])
ylabel('Reconstructed AMOCg'' [Sv]')
yyaxis right
bar(time_mly(31:end-30),abs(ideal.AMOCg_lp - ideal.AMOCg_wbeb_lp),'FaceColor',col1,'FaceAlpha',0.7)
hold on
bar(time_mly(31:end-30),abs(ideal.AMOCg_lp - moor_l8.AMOCg_lp),'FaceColor',col2,'FaceAlpha',0.7)
ylim([0 8])
yticks([0 1 2])
ylabel('Diff. to model truth [Sv]')

lgd = legend(h,'Model truth','Current approach (MTS Case 1)','Improvement (MTS Case 8)','Location','northoutside','Orientation','horizontal');
title(lgd,'5-year low-pass filtered','FontWeight','normal')
legend boxoff
set(gca,'FontSize',20)
xticks(datenum(1985:5:2020,1,1))
datetick('x','keepticks')
xlim([time_mly(31) time_mly(end-30)])
text(0, 1.01, '(c)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')


% - Sketches
nexttile(26,[1,15])  % hp version
fill([topo.lon_tracos fliplr(topo.lon_tracos)],[topo.waterdep_tracos ones(1,length(topo.waterdep_tracos)).*y_max],[.7 .7 .7],'EdgeColor','none')
hold on
fill([topo.lon_tracos-offset_east fliplr(topo.lon_tracos-offset_east)],[topo.waterdep_tracos ones(1,length(topo.waterdep_tracos)).*y_max],[.7 .7 .7],'EdgeColor','none')
fill([x1 x2-offset_east x2-offset_east x1],[0 0 y_max y_max],[.95 .95 .95],'EdgeColor','none')
plot([-35.82 moor_lon(5)-offset_east],[ref_level ref_level],'--','color',[.2 .2 .2],'LineWidth',1)
plot([-35.895 moor_lon(5)-offset_east],[bpr_Tg_dep(1) bpr_Tg_dep(1)],'--','color',[.2 .2 .2],'LineWidth',1)
plot([-35.887 moor_lon(5)-offset_east],[bpr_Tg_dep(2) bpr_Tg_dep(2)],'--','color',[.2 .2 .2],'LineWidth',1)

plot(bpr_lon(1:2),bpr_dep(1:2),'s','Color',col1,'LineWidth',3,'Markersize',15,'MarkerFaceColor',col1) % WB
plot([moor_lon(5)-offset_east moor_lon(5)-offset_east],[mc_eb1(3) moor_dep(5)],'k','LineWidth',2) % EB
plot(moor_lon(5)-offset_east,mc_eb1(end),'s','Color',col2,'LineWidth',3,'Markersize',18,'MarkerFaceColor',col2)
plot(moor_lon(5)-offset_east,mc_eb1(end),'w.','Markersize',mc_size+10)
plot([moor_lon(5) moor_lon(5) moor_lon(5) moor_lon(5) moor_lon(5)]-offset_east,mc_eb1(3:end),'.','Color',col2,'Markersize',mc_size)

h1(1) = plot(bpr_lon(1),3000,'ks','LineWidth',3,'Markersize',15,'MarkerFaceColor','k');
h1(2) = plot(moor_lon(5),3000,'k.','Markersize',mc_size);
legend(h1,'BPR','MicroCAT','Location','southeast')
legend boxoff

ylim([0 y_max]); axis ij
xlim([xmin xmax-offset_east])
yticks(0:500:2000)
yticklabels({'0 m', '500 m', '1000 m', '1500 m', '2000 m'})

xticks(-36:0.5:-31.5)
xticklabels({'36°W', '', '35°W', '', '', '', '12°E', '', '13°E', ''})
xtickangle(0)
set(gca, 'Layer', 'top')
set(gca,'fontsize',20)
text(0, 1.01, '(b)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')


nexttile(66,[1,15])  % lp version
fill([xmin moor_lon(4) moor_lon(4) xmin],[0 0 ref_level ref_level],col1,'EdgeColor','none')
hold on
fill([moor_lon(5) xmax xmax moor_lon(5)]-offset_east,[0 0 ref_level ref_level],col1,'EdgeColor','none')
fill([topo.lon_tracos fliplr(topo.lon_tracos)],[topo.waterdep_tracos ones(1,length(topo.waterdep_tracos)).*y_max],[.7 .7 .7],'EdgeColor','none')
fill([topo.lon_tracos-offset_east fliplr(topo.lon_tracos-offset_east)],[topo.waterdep_tracos ones(1,length(topo.waterdep_tracos)).*y_max],[.7 .7 .7],'EdgeColor','none')
fill([x1 x2-offset_east x2-offset_east x1],[0 0 y_max y_max],[.95 .95 .95],'EdgeColor','none')
plot([moor_lon(4) moor_lon(5)-offset_east],[ref_level ref_level],'--','color',[.2 .2 .2],'LineWidth',1)

plot([moor_lon(4) moor_lon(4)],[mc_wb4(2) moor_dep(4)],'k','LineWidth',2) % WB
plot(repmat(moor_lon(4),1,length(mc_wb4)),mc_wb4,'.','Color',col2,'MarkerSize',mc_size,'LineWidth',2)

plot([moor_lon(5)-offset_east moor_lon(5)-offset_east],[mc_eb1(2) moor_dep(5)],'k','LineWidth',2) % EB
plot(repmat(moor_lon(5),1,length(mc_eb1))-offset_east,mc_eb1,'.','Color',col2,'Markersize',mc_size)

h2 = plot(moor_lon(5),3000,'k.','Markersize',mc_size);
legend(h2,'MicroCAT','Location','southeast','NumColumns',2)
legend boxoff

ylim([0 y_max]); axis ij
xlim([xmin xmax-offset_east])
yticks(0:500:2000)
yticklabels({'0 m', '500 m', '1000 m', '1500 m', '2000 m'})

xticks(-36:0.5:-31.5)
xticklabels({'36°W', '', '35°W', '', '', '', '12°E', '', '13°E', ''})
xtickangle(0)
set(gca, 'Layer', 'top')
set(gca,'fontsize',20)
text(0, 1.01, '(d)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

set(gcf,'Position',[1 50 1440 720])
