%% Figure 1 - Overview TRACOS

section_tracos = load('/Users/ahans/Documents/PhD/DATA/Model/VIK20/interim_VIK_v4/v_topo_tracos_VIK.mat');
section_tracos.v_mean_tracos_filled = fillmissing(section_tracos.v_mean_tracos,'previous',2);

% define instrument positions
moor_lat = [-10.27, -10.38, -10.61, -10.94, -10.83]; % WB1, WB2, WB3, WB4, EB1
moor_lon = [-35.86+0.016, -35.68, -35.39, -34.99, 13.00-1/20];
moor_dep = [900, 2400, 3500, 4090, 1250];

dep_VIK = 1e+3 .* [0.0030    0.0095    0.0164    0.0239    0.0322    0.0415    0.0519    0.0639    0.0776    0.0936    0.1123    0.1343   0.1603    0.1911    0.2276    0.2709    0.3220    0.3821    0.4524    0.5340    0.6279    0.7347    0.8551    0.9892    1.1369    1.2977    1.4709    1.6555    1.8504    2.0544    2.2665    2.4854    2.7101    2.9398    3.1736    3.4108    3.6507    3.8929    4.1370    4.3827    4.6295    4.8773    5.1259    5.3752    5.6250    5.8751];

mc_wb1 = dep_VIK([20 21 23]);
mc_wb2 = dep_VIK([20 21 26 27 29 31]);
mc_wb3 = dep_VIK([20 21 29 33 36]);
mc_wb4 = dep_VIK([10 14 16 18 20 21 29 36 39]); % WB4ext: surface, 100m, 200m, 300m, 400m, 500m, 650m, 1900m
mc_wb4(end) = 4040; % adjust as bottom cell is actually thinner
mc_eb1 = dep_VIK([13 17 20 22 24 26]); % EB1ext: surface, 150m, 300m, 500m, 700m, 950m, 1200m, in observations: 300 500 700 950 1200
mc_eb1(end) = 1234;

adcps = [500 500 500 500 322];

cm_wb1 = dep_VIK(21);
cm_wb2 = dep_VIK([21 23 27 29]);
cm_wb3 = dep_VIK([21 23 27 29 32 34]);
cm_wb4 = dep_VIK([21 23 27 29 32 34]);
cm_eb1 = dep_VIK([21 23 25]);

bpr_lat = [ -10.23, -10.23, -10.27, -10.68, -10.71, -10.83]; % WBb5, WBb6, WBb7, EBb2, EBb3, EBb4
bpr_lon = [-35.87-0.018, -35.86-0.016, -35.86+0.028, 13.23, 13.19+0.01, 13+0.03-1/20+0.01];
bpr_dep = [272, 441, 1105, 295, 439, 1252-9];


% set properties for the lon-dep section
y_max = 4300;
clim_use = 5;
clevels = 25;
mc_size = 25;

% actual figure
fig = figure();
set(fig,'Position',[2 250 1440 620])
tiledlayout(3,3,'TileSpacing','compact','Padding','compact')

nexttile(1,[1 3])
insetplot.latlim = [-13.5 -7.5];
insetplot.lonlim = [-39 20];
m_proj('mercator','lat',insetplot.latlim,'lon',insetplot.lonlim);
[~,~] = m_etopo2('contourf',-6000:500:0,'edgecolor','none');
cmocean('-deep',12)
c = colorbar('east','Ticks',[-6000, -4000, -2000, 0],'TickLabels',{'6000 m', '4000 m', '2000 m', '0 m'});
caxis([-6000 0])
c.FontSize = 17;
hold on
m_plot(section_tracos.lon_tracos,section_tracos.lat_tracos,'w','LineWidth',4)
m_plot([-34.8 -34.8],[-12 -10],'w:','LineWidth',4)
m_plot([12.5 12.5],[-12 -10],'w:','LineWidth',4)
% m_plot([-34.8 -34.8],[-11.7 -10.3],'Color','k-','LineWidth',2)
% m_plot([12.5 12.5],[-11.7 -10.3],'Color','k-','LineWidth',2)
% m_plot([-34.8 -34.8],[-11.5 -10.5],'Color','w--','LineWidth',2)
% m_plot([12.5 12.5],[-11.5 -10.5],'Color','w--','LineWidth',2)
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_text(-38.7,-8.2,'Brazil', 'FontSize', 22, 'Color', 'k')
m_text(13.6,-8.2, 'Angola', 'FontSize', 22, 'Color', 'k')
m_grid('Fontsize',22,'linestyle','none','ticklen',1e-3); % with ticks
%m_grid('Fontsize',22,'tickdir','in','linestyle','none'); % without ticks
m_text(-39, -6.8, '(a)','FontSize',22,'Color','k')

nexttile(4,[2 1]) % WB
contourf(repmat(section_tracos.lon_tracos,46,1),section_tracos.dep',section_tracos.v_mean_tracos_filled'.*10,linspace(-clim_use,clim_use,clevels),'edgecolor','none')
caxis([-clim_use clim_use])
cmocean('bal',clevels-1)
hold on
line(cat(1,moor_lon,moor_lon),cat(1,adcps,moor_dep),'LineWidth',2,'Color','k')
line(cat(1,moor_lon(4),moor_lon(4)),cat(1,mc_wb4(1),adcps(4)),'LineWidth',2,'Color',[.45 .45 .45])
fill([section_tracos.lon_tracos'; flipud(section_tracos.lon_tracos')],[section_tracos.waterdep_tracos'; ones(length(section_tracos.waterdep_tracos),1).*y_max],[.7 .7 .7],'edgecolor','none')
% MCs
plot(repmat(moor_lon(1),1,length(mc_wb1)),mc_wb1,'k.','MarkerSize',mc_size,'LineWidth',3)
plot(repmat(moor_lon(2),1,length(mc_wb2)),mc_wb2,'k.','MarkerSize',mc_size,'LineWidth',3)
plot(repmat(moor_lon(3),1,length(mc_wb3)),mc_wb3,'k.','MarkerSize',mc_size,'LineWidth',3)
plot(repmat(moor_lon(4),1,length(mc_wb4(5:end))),mc_wb4(5:end),'k.','MarkerSize',mc_size,'LineWidth',3)
plot(repmat(moor_lon(4),1,4),mc_wb4(1:4),'.','Color',[139 139 0]./255,'MarkerSize',mc_size,'LineWidth',3)
%Aqds
plot(repmat(moor_lon(1),1,length(cm_wb1)),cm_wb1,'kx','MarkerSize',10,'LineWidth',2)
plot(repmat(moor_lon(2),1,length(cm_wb2)),cm_wb2,'kx','MarkerSize',10,'LineWidth',2)
plot(repmat(moor_lon(3),1,length(cm_wb3)),cm_wb3,'kx','MarkerSize',10,'LineWidth',2)
plot(repmat(moor_lon(4),1,length(cm_wb4)),cm_wb4,'kx','MarkerSize',10,'LineWidth',2)
% Adcps
plot(moor_lon,adcps,'k*','MarkerSize',15,'LineWidth',2)

scatter(bpr_lon,bpr_dep,200,'ks','filled')
scatter(bpr_lon(3),bpr_dep(3),200,[139 139 0]./255,'s','filled')

text(moor_lon(1)-0.15,moor_dep(1)+35,'WB1', 'FontSize', 22, 'Color', 'k')
text(moor_lon(2)-0.15,moor_dep(2)+55,'WB2', 'FontSize', 22, 'Color', 'k')
text(moor_lon(3)-0.15,moor_dep(3)+75,'WB3', 'FontSize', 22, 'Color', 'k')
text(moor_lon(4)-0.15,moor_dep(4)+85,'WB4', 'FontSize', 22, 'Color', 'k')
text(bpr_lon(1)-0.2,bpr_dep(1)+20,'WBb5', 'FontSize', 22, 'Color', 'k')
text(bpr_lon(2)-0.2,bpr_dep(2)+120,'WBb6', 'FontSize', 22, 'Color', 'k')
text(bpr_lon(3)-0.18,bpr_dep(3)+180,'WBb7', 'FontSize', 22, 'Color', [100 100 0]./255)

ylim([0 y_max]); axis ij
xlim([-36 -34.8])
yticks(0:500:4000)
yticklabels({'0 m', '', '1000 m', '', '2000 m', '', '3000 m', '', '4000 m'})
xticks(-36:0.5:-35)
xticklabels({'36°W', '35.5°W', '35°W'})
set(gca,'fontsize',22)
text(0, 1.01, '(b)', 'Units','normalized', 'FontSize',22, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(5,[2 1]) % Interior
contourf(repmat(section_tracos.lon_tracos,46,1),section_tracos.dep',section_tracos.v_mean_tracos_filled'.*10,linspace(-clim_use,clim_use,clevels),'edgecolor','none')
caxis([-clim_use clim_use])
cmocean('bal',clevels-1)
hold on
fill([section_tracos.lon_tracos'; flipud(section_tracos.lon_tracos')],[section_tracos.waterdep_tracos'; ones(length(section_tracos.waterdep_tracos),1).*y_max],[.7 .7 .7],'edgecolor','none')

ylim([0 y_max]); axis ij
xlim([-34.8 12.5])
set(gca,'YTickLabel',[])
xticks(-25:20:10)
xticklabels({'25°W', '0°'})
set(gca,'fontsize',22)
text(0, 1.01, '(c)', 'Units','normalized', 'FontSize',22, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

nexttile(6,[2 1]) % EB
contourf(repmat(section_tracos.lon_tracos,46,1),section_tracos.dep',section_tracos.v_mean_tracos_filled'.*10,linspace(-clim_use,clim_use,clevels),'edgecolor','none')
caxis([-clim_use clim_use])
cmocean('bal',clevels-1)
hold on
line(cat(1,moor_lon,moor_lon),cat(1,adcps,moor_dep),'LineWidth',2,'Color','k')
line(cat(1,moor_lon(5),moor_lon(5)),cat(1,mc_eb1(1),adcps(5)),'LineWidth',2,'Color',[.45 .45 .45])
fill([section_tracos.lon_tracos'; flipud(section_tracos.lon_tracos')],[section_tracos.waterdep_tracos'; ones(length(section_tracos.waterdep_tracos),1).*y_max],[.7 .7 .7],'edgecolor','none')
h(1) = plot(0,0,'ks','Markersize',15,'MarkerFaceColor','k'); % line only for legend
h(2) = plot(repmat(moor_lon(5),1,length(mc_eb1(2:end))),mc_eb1(2:end),'k.','MarkerSize',mc_size,'LineWidth',2);
plot(moor_lon(5),mc_eb1(1),'.','Color',[139 139 0]./255,'MarkerSize',mc_size,'LineWidth',2)
h(3) = plot(repmat(moor_lon(5),1,length(cm_eb1)),cm_eb1,'kx','MarkerSize',10,'LineWidth',2);
h(4) = plot(moor_lon(5),adcps(5),'k*','MarkerSize',15,'LineWidth',2);
scatter(bpr_lon,bpr_dep,200,'ks','filled')
scatter(bpr_lon(6),bpr_dep(6),200,[139 139 0]./255,'s','filled')
text(moor_lon(5)-0.07,moor_dep(5)+180,'EB1', 'FontSize', 22, 'Color', 'k')
text(bpr_lon(4)+0.03,bpr_dep(4),'EBb2', 'FontSize', 22, 'Color', 'k')
text(bpr_lon(5)+0.03,bpr_dep(5)+100,'EBb3', 'FontSize', 22, 'Color', 'k')
text(bpr_lon(6)+0.03,bpr_dep(6)+20,'EBb4', 'FontSize', 22, 'Color', [100 100 0]./255)

ylim([0 y_max]); axis ij
xlim([12.5 13.7])
set(gca,'YTickLabel',[]);
xticks(12.5:0.5:13.5)
xticklabels({'12.5°E','13°E', '13.5°E'})
set(gca,'fontsize',22)
text(0, 1.01, '(d)', 'Units','normalized', 'FontSize',22, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
legend(h,'PIES/BPR','MicroCAT','Current Meter','Adcp','Location','southwest')
legend boxoff

cb = colorbar;
cb.Location = 'east';
cb.Label.String = 'Meridional velocity [10 cm s^-^1]';
cb.Label.FontSize = 22;

img = imread('/Users/ahans/Documents/PhD/PLOTS/Model/Overview/globe_11Sline_grey.png');
axes('Position',[0.004 0.06 0.2 0.2])
imshow(img)

%% globe to overlay

insetplot.latlim = [-11 -11];
insetplot.lonlim = [-37 13.7];

fig = figure();
m_proj('ortho','lat',-11','long',-15');
m_coast('patch',[.7 .7 .7]);
m_grid('linest','-','xticklabels',[],'yticklabels',[]);
hold on
m_plot(insetplot.lonlim,insetplot.latlim,'k','Linewidth',3)
set(gcf,'Color',[.7 .7 .7])
set(gcf,'Position',[400 500 100 100])
exportgraphics(fig, '/Users/ahans/Documents/PhD/PLOTS/Model/Overview/globe_11Sline_grey.png', 'BackgroundColor', 'current','Resolution',600);
