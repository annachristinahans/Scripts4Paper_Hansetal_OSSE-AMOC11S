%% Figure 6 - Fits for IES method
% overview of the relation between travel time and rho, bp and foff

% --- loading and preparing data ---
load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/ies_test_tauts_VIK.mat')
load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/ies_test_bp_VIK.mat')
bp_wb300 = (bp_dyn_wb300 + bp_ssh_wb300) .*1e-4;
bp_wb500 = (bp_dyn_wb500 + bp_ssh_wb500) .*1e-4;
bp_eb1200 = (bp_dyn_eb1200 + bp_ssh_eb1200) .*1e-4;
load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/ies_test_foff_VIK.mat')

% average to annual data to omit seasonal noise
[yy,~,~] = datevec(time);
yyears = unique(yy);
%clear tau_wb300_yly rho_wb300_yly tau_wb500_yly rho_wb500_yly tau_eb1200_yly rho_eb1200_yly
%clear bp_wb300_yly bp_wb500_yly bp_eb1200_yly fof3_wb300_yly fof3_wb500_yly fof3_eb1200_yly fof5_wb500_yly fof5_eb1200_yly
for i = 1:length(yyears)
    tau_wb300_yly(i) = mean(tau_wb300(yy==yyears(i)));
    tau_wb500_yly(i) = mean(tau_wb500(yy==yyears(i)));
    tau_eb1200_yly(i) = mean(tau_eb1200(yy==yyears(i)));
    rho_wb300_yly(:,i) = mean(rho_wb300(:,yy==yyears(i)),2);
    rho_wb500_yly(:,i) = mean(rho_wb500(:,yy==yyears(i)),2);
    rho_eb1200_yly(:,i) = mean(rho_eb1200(:,yy==yyears(i)),2);
    tau_wb300_yly(i) = mean(tau_wb300(yy==yyears(i)));
    bp_wb300_yly(i) = mean(bp_wb300(yy==yyears(i)));
    bp_wb500_yly(i) = mean(bp_wb500(yy==yyears(i)));
    bp_eb1200_yly(i) = mean(bp_eb1200(yy==yyears(i)));
    fof3_wb300_yly(i) = mean(fof3_wb300(yy==yyears(i)));
    fof3_wb500_yly(i) = mean(fof3_wb500(yy==yyears(i)));
    fof3_eb1200_yly(i) = mean(fof3_eb1200(yy==yyears(i)));
    fof5_wb500_yly(i) = mean(fof5_wb500(yy==yyears(i)));
    fof5_eb1200_yly(i) = mean(fof5_eb1200(yy==yyears(i)));
end

% depth to diyplay for density
idep_wb300 = [3 10 15];
idep_wb500 = [10 15 18];
idep_eb1200 = [10 15 26];

% regular travel time vector
tau_q_wb300 = min(tau_wb300_yly):0.0001:max(tau_wb300_yly)+0.0001;
tau_q_wb500 = min(tau_wb500_yly):0.0001:max(tau_wb500_yly)+0.0001;
tau_q_eb1200 = min(tau_eb1200_yly):0.0001:max(tau_eb1200_yly)+0.0001;

% --- plotting ---
figure();
t = tiledlayout(5,3,'TileSpacing','compact','Padding','compact');
fig_lab = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
fig_lab_count = 1;

% rho
for i = 1:3
    nexttile % WB300
    [p,S] = polyfit(tau_wb300_yly,rho_wb300_yly(idep_wb300(i),:),3); % cubic polynomial
    [rho_fit, ~] = polyval(p,tau_q_wb300,S);
    rmse_fit = rmse(rho_wb300_yly(idep_wb300(i),:),polyval(p,tau_wb300_yly)).*1e+2;

    plot(tau_wb300_yly,rho_wb300_yly(idep_wb300(i),:),'kx','LineWidth',2)
    hold on
    plot(tau_q_wb300,rho_fit,'r','LineWidth',2)
    title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^-^2 kg m^-^3'],'FontWeight', 'normal')
    xlim([min(tau_q_wb300) max(tau_q_wb300)])
    ylabel('\rho [kg m^-^3]')
    text(0.02, 0.9, [num2str(round(dep_wb300(idep_wb300(i)))) ' m'], 'Units', 'normalized','FontSize',18)
    set(gca,'FontSize',20)
    text(0, 1.01, fig_lab{fig_lab_count}, 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
    fig_lab_count = fig_lab_count + 1;

    nexttile % WB500
    [p,S] = polyfit(tau_wb500_yly,rho_wb500_yly(idep_wb500(i),:),3); % cubic polynomial
    [rho_fit, ~] = polyval(p,tau_q_wb500,S);
    rmse_fit = rmse(rho_wb500_yly(idep_wb500(i),:),polyval(p,tau_wb500_yly)).*1e+2;

    plot(tau_wb500_yly,rho_wb500_yly(idep_wb500(i),:),'kx','LineWidth',2)
    hold on
    plot(tau_q_wb500,rho_fit,'r','LineWidth',2)
    title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^-^2 kg m^-^3'],'FontWeight', 'normal')
    xlim([min(tau_q_wb500) max(tau_q_wb500)])
    text(0.02, 0.9, [num2str(round(dep_wb500(idep_wb500(i)))) ' m'], 'Units', 'normalized','FontSize',18)
    set(gca,'FontSize',20)
    text(0, 1.01, fig_lab{fig_lab_count}, 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
    fig_lab_count = fig_lab_count + 1;

    nexttile % EB1200
    [p,S] = polyfit(tau_eb1200_yly,rho_eb1200_yly(idep_eb1200(i),:),3); % cubic polynomial
    [rho_fit, ~] = polyval(p,tau_q_eb1200,S);
    rmse_fit = rmse(rho_eb1200_yly(idep_eb1200(i),:),polyval(p,tau_eb1200_yly)).*1e+2;

    plot(tau_eb1200_yly,rho_eb1200_yly(idep_eb1200(i),:),'kx','LineWidth',2)
    hold on
    plot(tau_q_eb1200,rho_fit,'r','LineWidth',2)
    title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^-^2 kg m^-^3'],'FontWeight', 'normal')
    xlim([min(tau_q_eb1200) max(tau_q_eb1200)])
    text(0.02, 0.9, [num2str(round(dep_eb1200(idep_eb1200(i)))) ' m'], 'Units', 'normalized','FontSize',18)
    set(gca,'FontSize',20)
    text(0, 1.01, fig_lab{fig_lab_count}, 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
    fig_lab_count = fig_lab_count + 1;
end

% bp
nexttile
[p,S] = polyfit(tau_wb300_yly,bp_wb300_yly,3); % cubic polynomial
[bp_fit, ~] = polyval(p,tau_q_wb300,S);
rmse_fit = rmse(bp_wb300_yly,polyval(p,tau_wb300_yly)).*1e+3;

tau_wb300_mdl = fitlm(datenum(yyears(yyears>=1994),1,1), tau_wb300_yly(yyears>=1994));
tau_wb300_trend = predict(tau_wb300_mdl,datenum([1994 2023],1,1)');
bp_wb300_mdl = fitlm(datenum(yyears(yyears>=1994),1,1), bp_wb300_yly(yyears>=1994));
bp_wb300_trend = predict(bp_wb300_mdl,datenum([1994 2023],1,1)');

plot(tau_wb300_yly,bp_wb300_yly,'kx','LineWidth',2)
hold on
plot(tau_wb300_trend,bp_wb300_trend,':','Color',[255 193 37]./255,'LineWidth',6)
plot(tau_q_wb300,bp_fit,'r','LineWidth',2)
title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^-^3 dbar'],'FontWeight', 'normal')
ylabel('p_{bottom} [dbar]')
xlim([min(tau_q_wb300) max(tau_q_wb300)])
set(gca,'FontSize',20)
text(0, 1.01, '(j)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
[p,S] = polyfit(tau_wb500_yly,bp_wb500_yly,3); % cubic polynomial
[bp_fit, ~] = polyval(p,tau_q_wb500,S);
rmse_fit = rmse(bp_wb500_yly,polyval(p,tau_wb500_yly)).*1e+3;

tau_wb500_mdl = fitlm(datenum(yyears(yyears>=1994),1,1), tau_wb500_yly(yyears>=1994));
tau_wb500_trend = predict(tau_wb500_mdl,datenum([1994 2023],1,1)');
bp_wb500_mdl = fitlm(datenum(yyears(yyears>=1994),1,1), bp_wb500_yly(yyears>=1994));
bp_wb500_trend = predict(bp_wb500_mdl,datenum([1994 2023],1,1)');

plot(tau_wb500_yly,bp_wb500_yly,'kx','LineWidth',2)
hold on
plot(tau_wb500_trend,bp_wb500_trend,':','Color',[255 193 37]./255,'LineWidth',6)
plot(tau_q_wb500,bp_fit,'r','LineWidth',2)
title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^-^3 dbar'],'FontWeight', 'normal')
xlim([min(tau_q_wb500) max(tau_q_wb500)])
set(gca,'FontSize',20)
text(0, 1.01, '(k)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
[p,S] = polyfit(tau_eb1200_yly,bp_eb1200_yly,3); % cubic polynomial
[bp_fit, ~] = polyval(p,tau_q_eb1200,S);
rmse_fit = rmse(bp_eb1200_yly,polyval(p,tau_eb1200_yly)).*1e+3;

tau_eb1200_mdl = fitlm(datenum(yyears(yyears>=1994),1,1), tau_eb1200_yly(yyears>=1994));
tau_eb1200_trend = predict(tau_eb1200_mdl,datenum([1994 2023],1,1)');
bp_eb1200_mdl = fitlm(datenum(yyears(yyears>=1994),1,1), bp_eb1200_yly(yyears>=1994));
bp_eb1200_trend = predict(bp_eb1200_mdl,datenum([1994 2023],1,1)');

plot(tau_eb1200_yly,bp_eb1200_yly,'kx','LineWidth',2)
hold on
plot(tau_eb1200_trend,bp_eb1200_trend,':','Color',[255 193 37]./255,'LineWidth',6)
plot(tau_q_eb1200,bp_fit,'r','LineWidth',2)
title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^-^3 dbar'],'FontWeight', 'normal')
xlim([min(tau_q_eb1200) max(tau_q_eb1200)])
set(gca,'FontSize',20)
text(0, 1.01, '(l)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

% fof
nexttile
[p,S] = polyfit(tau_wb300_yly,fof3_wb300_yly,3); % cubic polynomial
[fof_fit, ~] = polyval(p,tau_q_wb300,S);
rmse_fit = rmse(fof3_wb300_yly,polyval(p,tau_wb300_yly)).*1e-4;

plot(tau_wb300_yly,fof3_wb300_yly.*1e-8,'kx','LineWidth',2)
hold on
plot(tau_q_wb300,fof_fit.*1e-8,'r','LineWidth',2)
title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^4 J m^-^2'],'FontWeight', 'normal')
ylabel('\chi_{270m} [10^8 J m^-^2]')
xlim([min(tau_q_wb300) max(tau_q_wb300)])
xlabel('\tau_{WBb5} [s]')
set(gca,'FontSize',20)
text(0, 1.01, '(m)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
[p,S] = polyfit(tau_wb500_yly,fof3_wb500_yly,3); % cubic polynomial
[fof_fit, ~] = polyval(p,tau_q_wb500,S);
rmse_fit = rmse(fof3_wb500_yly,polyval(p,tau_wb500_yly)).*1e-4;

plot(tau_wb500_yly,fof3_wb500_yly.*1e-8,'kx','LineWidth',2)
hold on
plot(tau_q_wb500,fof_fit.*1e-8,'r','LineWidth',2)
title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^4 J m^-^2'],'FontWeight', 'normal')
xlim([min(tau_q_wb500) max(tau_q_wb500)])
xlabel('\tau_{WBb6} [s]')
set(gca,'FontSize',20)
text(0, 1.01, '(n)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile
[p,S] = polyfit(tau_eb1200_yly,fof3_eb1200_yly,3); % cubic polynomial
[fof_fit, ~] = polyval(p,tau_q_eb1200,S);
rmse_fit = rmse(fof3_eb1200_yly,polyval(p,tau_eb1200_yly)).*1e-4;

plot(tau_eb1200_yly,fof3_eb1200_yly.*1e-8,'kx','LineWidth',2)
hold on
plot(tau_q_eb1200,fof_fit.*1e-8,'r','LineWidth',2)
title(['rmse: ' num2str(rmse_fit,'%.0f') ' x10^4 J m^-^2'],'FontWeight', 'normal')
xlim([min(tau_q_eb1200) max(tau_q_eb1200)])
xlabel('\tau_{EBb4} [s]')
set(gca,'FontSize',20)
text(0, 1.01, '(o)', 'Units','normalized', 'FontSize',20, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

set(gcf,'Position',[1500 76 1350 1100])


%% Auxilary
% --- related statistics (for text) ---
% mean density error

for idep = 1:16
    [p,~] = polyfit(tau_wb300_yly,rho_wb300_yly(idep,:),3); % cubic polynomial
    rho_diff = rho_wb300_yly(idep,:) - polyval(p,tau_wb300_yly);
    mean_err_wb300(idep) = mean(abs(rho_diff));
    
    [p,~] = polyfit(tau_wb500_yly,rho_wb500_yly(idep,:),3); % cubic polynomial
    rho_diff = rho_wb500_yly(idep,:) - polyval(p,tau_wb500_yly);
    mean_err_wb500(idep) = mean(abs(rho_diff));
    
    [p,~] = polyfit(tau_eb1200_yly,rho_eb1200_yly(idep,:),3); % cubic polynomial
    rho_diff = rho_eb1200_yly(idep,:) - polyval(p,tau_eb1200_yly);
    mean_err_eb1200(idep) = mean(abs(rho_diff));
end

figure();
plot(mean_err_wb300,dep_wb300(1:16))
hold on
plot(mean_err_wb500,dep_wb500(1:16))
plot(mean_err_eb1200,dep_eb1200(1:16))
axis ij
