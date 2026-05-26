%% Figure 4 - Comparing Tg for BPR method

% --- loading and preparing data ---
nmid = 'VIK';

% ideal setup for comparison
ideal = load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/vg_1100S_' nmid '.mat']);
[yy,mm,~] = datevec(ideal.time);
ideal.time_mly = unique(datenum(yy,mm,15));
for i = 1:length(ideal.time_mly)
    ideal.Tg_basin_mly(:,i) = mean(ideal.Tg_basin(:,datenum(yy,mm,15)==ideal.time_mly(i)),2);
end
ideal.Tg_basin_mly = ideal.Tg_basin_mly - mean(ideal.Tg_basin_mly,2);

myColors = [
    224 201 0;  % yellow
    250 150 50; % orange
    15 160 235; % soft blue, original: 0 189 223
    250 150 50; % orange
    15 160 235; % soft blue
]./255;

% --- plotting ---
figure();
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for i = 1:5
    if i == 1
        load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc0_ebb0_ref1_vst1_stp1_' nmid '.mat'])
        %nid1 = 'eof_{no}, 2 BPR pairs';
        nid1 = '2 BPR pairs, EOF_{no}';
        LineStyleUse = '-';
    elseif i == 2
        load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc0_ebb0_ref1_vst2_stp1_' nmid '.mat'])
        %nid2 = 'eof_{known}, 2 BPR pairs';
        nid2 = '2 BPR pairs, EOF_{known}';
        LineStyleUse = '-';
    elseif i == 3
        load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc0_ebb0_ref1_vst3_stp1_' nmid '.mat'])
        %nid3 = 'linear, 2 BPR pairs';
        nid3 = '2 BPR pairs, Linear';
        LineStyleUse = '-';
    elseif i == 4
        load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc0_ebb0_ref1_vst2_stp4_' nmid '.mat'])
        %nid4 = 'eof_{known}, 3 BPR pairs';
        nid4 = '3 BPR pairs, EOF_{known}';
        LineStyleUse = 'o';
    elseif i == 5
        load(['/Users/ahans/Documents/PhD/DATA/Model/' nmid '20/interim_' nmid '_v4/bpr_pos22_prc0_ebb0_ref1_vst3_stp4_' nmid '.mat'])
        %nid5 = 'linear, 3 BPR pairs';
        nid5 = '3 BPR pairs, Linear';
        LineStyleUse = 'o';
    end
    
    nexttile(1) % corr (monthly or lp time series?????)
    for iz = 1:length(dep_grid)
        r = corrcoef(ideal.Tg_basin_mly(iz,:),Tg_basin(iz,:));
        r_mly(iz) = r(1,2);
    end
    hold on
    if i > 3
        plot(r_mly(1:10:end),dep_grid(1:10:end),LineStyleUse,'LineWidth',3,'Color',myColors(i,:))
    else
        plot(r_mly,dep_grid,LineStyleUse,'LineWidth',3,'Color',myColors(i,:))
    end
    
    nexttile(2) % rmse
    rmse_mly = (sqrt(sum((ideal.Tg_basin_mly(ideal.dep_grid<=max(dep_grid),:)- Tg_basin).^2,2) ./length(time_mly))).*1e-3;
    hold on
    if i > 3
        plot(rmse_mly(1:10:end),dep_grid(1:10:end),LineStyleUse,'LineWidth',3,'Color',myColors(i,:))
    else
        plot(rmse_mly,dep_grid,LineStyleUse,'LineWidth',3,'Color',myColors(i,:))
    end

    nexttile(3) % harmonic amplitude
    for iz = 1:length(dep_grid)
        [~,a,~] = anharm_edt(time_mly, datenum(2013,1,1):datenum(2014,1,1), Tg_basin(iz,:),2);
        a_bpr(:,iz) = a;
    end
    hold on
    if i > 3
        plot(a_bpr(2,1:10:end).*1e-3,dep_grid(1:10:end),LineStyleUse,'LineWidth',3,'Color',myColors(i,:))
    else
        plot(a_bpr(2,:).*1e-3,dep_grid,LineStyleUse,'LineWidth',3,'Color',myColors(i,:))
    end

    nexttile(4) % harmonic phase
    % include NaNs when there is a phase shift
    phase_diff = diff(a_bpr(4,:));
    a_bpr(4,abs(phase_diff)>100) = NaN;
    a_bpr(4,find(a_bpr(2,:)==0)) = NaN;
    hold on
    if i > 3
        plot(a_bpr(4,1:10:end),dep_grid(1:10:end),LineStyleUse,'Color',myColors(i,:),'LineWidth',3)
        plot(a_bpr(4,1:10:end)-365,dep_grid(1:10:end),LineStyleUse,'Color',myColors(i,:),'LineWidth',3)
        plot(a_bpr(4,1:10:end)+365,dep_grid(1:10:end),LineStyleUse,'Color',myColors(i,:),'LineWidth',3)
    else
        plot(a_bpr(4,:),dep_grid,LineStyleUse,'Color',myColors(i,:),'LineWidth',3)
        plot(a_bpr(4,:)-365,dep_grid,LineStyleUse,'Color',myColors(i,:),'LineWidth',3)
        plot(a_bpr(4,:)+365,dep_grid,LineStyleUse,'Color',myColors(i,:),'LineWidth',3)
    end

end

nexttile(1)
yline(280)
yline(440)
axis ij
ylim([0 max(dep_grid)])
xlabel('T_g'' Pearson correlation coefficient')
ylabel('Depth [m]')
set(gca,'FontSize',18)
box on
text(0.015, 0.99, '(a)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','top')

nexttile(2)
yline(280)
yline(440)
axis ij
ylim([0 max(dep_grid)])
xlabel('T_g'' rmse [Sv km^-^1]')
set(gca,'FontSize',18)
box on
text(0.015, 0.99, '(b)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','top')

nexttile(3)
for iz = 1:length(dep_grid)
    [~,a,~] = anharm_edt(time_mly, datenum(2013,1,1):datenum(2014,1,1), ideal.Tg_basin_mly(iz,:),2);
    a_idl(:,iz) = a;
end
plot(a_idl(2,:).*1e-3,dep_grid,'k','LineWidth',3)
yline(280)
yline(440)
axis ij
legend(nid1,nid2,nid3,nid4,nid5,'Model truth','Location','southeast')
legend('boxoff')
ylim([0 max(dep_grid)])
xlabel('T_g'' Annual harmonic amplitude [Sv km^-^1]')
ylabel('Depth [m]')
set(gca,'FontSize',18)
box on
text(0.015, 0.99, '(c)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','top')

nexttile(4)
phase_diff = diff(a_idl(4,:));
a_idl(4,abs(phase_diff)>100) = NaN;
plot(a_idl(4,:),dep_grid,'k','LineWidth',3)
plot(a_idl(4,:)-365,dep_grid,'k','LineWidth',3)
plot(a_idl(4,:)+365,dep_grid,'k','LineWidth',3)
xline(0)
xline(365)
yline(280)
yline(440)
axis ij
xlim([-180 365+180])
ylim([0 max(dep_grid)])
xlabel('T_g'' Annual harmonic phase [days]')
set(gca,'FontSize',18)
box on
text(0.015, 0.99, '(d)', 'Units','normalized', 'FontSize',18, 'HorizontalAlignment','left', 'VerticalAlignment','top')

set(gcf,'Position',[50 50 1300 700])
