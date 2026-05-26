%% INALT/VIKING IES Test - cut
% to run on nesh

% add a suffix in the saved files in code testing cases
test_suffix = ''; % '' if nothing is tested

% choose model: INA (1) or VIK (2)
mod_option = 2;
if mod_option == 1
    num_yy = 40;
    nmid = 'INA';
elseif mod_option == 2
    num_yy = 44;
    nmid = 'VIK';
end

addpath('/gxfs_home/geomar/smomw628/matlab_toolbox/seawater_ver3_3.1/')

% --- loading model data
if mod_option == 1
        files1 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*T_11S_box.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*T_11S_box.ncbotpressure_detrssh.nc');
        mmask = '/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/1_INALT20.L46-KFS119_mesh_mask_11S_box.nc';
elseif mod_option == 2
        files1 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*T_11S_box.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*T_11S_box.ncbotpressure_detrssh.nc');
        mmask = '/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/1_VIKING20X.L46-KFS003_mesh_mask_11S_box.nc';
end

lat = double(ncread(mmask,'gphit'));
lat = lat(1,:);
lon = double(ncread(mmask,'glamt'));
lon = lon(:,1);
dep = double(ncread(mmask,'gdept_0')); % in m, including partial cells
dep_cell_size = double(ncread(mmask,'e3t_0'));
qc = double(ncread(mmask,'tmask')); % 0 land, 1 water
qc(qc==0) = NaN;

% choose positions of the WB IESs
[~,ilat_wb300] = min(abs(lat+10.23));
[~,ilon_wb300] = min(abs(lon+35.87));
[~,ilat_wb500] = min(abs(lat+10.23));
[~,ilon_wb500] = min(abs(lon+35.86));
ilon_wb300 = ilon_wb300 - 2; % modification to fit actual BPR depth better
ilat_wb300 = ilat_wb300 - 2;
ilon_wb500 = ilon_wb500 - 1;
ilat_wb500 = ilat_wb500 - 1;
dep_wb300 = squeeze(dep(ilon_wb300,ilat_wb300,1:16));
cellsize_wb300 = squeeze(dep_cell_size(ilon_wb300,ilat_wb300,1:16));
dep_wb500 = squeeze(dep(ilon_wb500,ilat_wb500,1:19));
cellsize_wb500 = squeeze(dep_cell_size(ilon_wb500,ilat_wb500,1:19));

% choose position of K4 mooring
[~,ilat_k4] = min(abs(lat+10.94));
[~,ilon_k4] = min(abs(lon+34.99));
dep_k4 = squeeze(dep(ilon_k4,ilat_k4,1:39));
cellsize_k4 = squeeze(dep_cell_size(ilon_k4,ilat_k4,1:39));

% choose position of Angola mooring
[~,ilat_eb1200] = min(abs(lat+10.83));
[~,ilon_eb1200] = min(abs(lon-13.00));
ilon_eb1200 = ilon_eb1200 - 1; % so that the mooring reaches a depth of 1200m
dep_eb1200 = squeeze(dep(ilon_eb1200,ilat_eb1200,1:26));
cellsize_eb1200 = squeeze(dep_cell_size(ilon_eb1200,ilat_eb1200,1:26));


% actual loading of data
countt = 1;
for i = 1:num_yy % loop over years

    tim = double(ncread([files1(i).folder '/' files1(i).name],'time_counter'));
    time(countt:countt-1+length(tim)) = datenum(1900,1,1,0,0,tim);
    
    ptemp = double(ncread([files1(i).folder '/' files1(i).name],'votemper')); % in °C, potential temperature
    salt = double(ncread([files1(i).folder '/' files1(i).name],'vosaline')); % in 0.001, practical salinity
    ssh = double(ncread([files2(i).folder '/' files2(i).name],'sossheig')); % in m
    
    % applying land-sea mask
    ptemp = ptemp .* qc;
    salt = salt .* qc;
    ssh = ssh .* qc(:,:,1);

    % choosing relevant positions
    ptemp_wb300(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_wb300,ilat_wb300,1:16,:));
    salt_wb300(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_wb300,ilat_wb300,1:16,:));
    ssh_wb300(countt:countt-1+length(tim)) = squeeze(ssh(ilon_wb300,ilat_wb300,:));
    ptemp_wb500(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_wb500,ilat_wb500,1:19,:));
    salt_wb500(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_wb500,ilat_wb500,1:19,:));
    ssh_wb500(countt:countt-1+length(tim)) = squeeze(ssh(ilon_wb500,ilat_wb500,:));
    ptemp_k4(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_k4,ilat_k4,1:39,:));
    salt_k4(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_k4,ilat_k4,1:39,:));
    ssh_k4(countt:countt-1+length(tim)) = squeeze(ssh(ilon_k4,ilat_k4,:));
    ptemp_eb1200(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_eb1200,ilat_eb1200,1:26,:));
    salt_eb1200(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_eb1200,ilat_eb1200,1:26,:));
    ssh_eb1200(countt:countt-1+length(tim)) = squeeze(ssh(ilon_eb1200,ilat_eb1200,:));

    countt = countt + length(tim);
end
clear i files1 files2 mmask tim qc ptemp salt ssh countt

% --- processing and simulating acoustic travel times

% defining global variables
g = sw_g(-10.5,0);

pres_wb300 = sw_pres(dep_wb300,lat(ilat_wb300));
pres_wb500 = sw_pres(dep_wb500,lat(ilat_wb500));
pres_k4 = sw_pres(dep_k4,lat(ilat_k4));
pres_eb1200 = sw_pres(dep_eb1200,lat(ilat_eb1200));

temp_wb300 = sw_temp(salt_wb300, ptemp_wb300, pres_wb300,0);
temp_wb500 = sw_temp(salt_wb500, ptemp_wb500, pres_wb500,0);
temp_k4 = sw_temp(salt_k4, ptemp_k4, pres_k4,0);
temp_eb1200 = sw_temp(salt_eb1200, ptemp_eb1200, pres_eb1200,0);

% computing sound velocity (in m/s) for model bins
svel_wb300 = sw_svel(salt_wb300, temp_wb300, pres_wb300);
svel_wb500 = sw_svel(salt_wb500, temp_wb500, pres_wb500);
svel_k4 = sw_svel(salt_k4, temp_k4, pres_k4);
svel_eb1200 = sw_svel(salt_eb1200, temp_eb1200, pres_eb1200);

% deriving round trip travel time
tau_wb300 = sum(cellsize_wb300 ./ svel_wb300,1).*2;
tau_wb500 = sum(cellsize_wb500 ./ svel_wb500,1).*2;
tau_k4 = sum(cellsize_k4 ./ svel_k4,1).*2;
tau_eb1200 = sum(cellsize_eb1200 ./ svel_eb1200,1).*2;

% computing specific volume anomaly
rho_wb300 = sw_dens(salt_wb300, temp_wb300, pres_wb300);
rho_wb500 = sw_dens(salt_wb500, temp_wb500, pres_wb500);
rho_k4 = sw_dens(salt_k4, temp_k4, pres_k4);
rho_eb1200 = sw_dens(salt_eb1200, temp_eb1200, pres_eb1200);

% computing density
svan_wb300 = sw_svan(salt_wb300, temp_wb300, pres_wb300);
svan_wb500 = sw_svan(salt_wb500, temp_wb500, pres_wb500);
svan_k4 = sw_svan(salt_k4, temp_k4, pres_k4);
svan_eb1200 = sw_svan(salt_eb1200, temp_eb1200, pres_eb1200);

% computing geopotential anomaly
gpan_wb300 = sw_gpan(salt_wb300, temp_wb300, pres_wb300);
gpan_wb500 = sw_gpan(salt_wb500, temp_wb500, pres_wb500);
gpan_k4 = sw_gpan(salt_k4, temp_k4, pres_k4);
gpan_eb1200 = sw_gpan(salt_eb1200, temp_eb1200, pres_eb1200);


% interpolating on 10m grid
dep_wb300_grid = 0:10:270;
dep_wb500_grid = 0:10:440;
dep_k4_grid = 0:10:4100;
dep_eb1200_grid = 0:10:1250;
pres_wb300_grid = sw_pres(dep_wb300_grid,lat(ilat_wb300));
pres_wb500_grid = sw_pres(dep_wb500_grid,lat(ilat_wb500));
pres_k4_grid = sw_pres(dep_k4_grid,lat(ilat_k4));
pres_eb1200_grid = sw_pres(dep_eb1200_grid,lat(ilat_eb1200));

% for temperature and salinity
temp_wb300_grid = interp1(dep_wb300,temp_wb300,dep_wb300_grid,'linear','extrap');
temp_wb300_grid(1,:) = temp_wb300_grid(2,:); % assuming a mixed layer in the upper 10m
salt_wb300_grid = interp1(dep_wb300,salt_wb300,dep_wb300_grid,'linear','extrap');
salt_wb300_grid(1,:) = salt_wb300_grid(2,:);

temp_wb500_grid = interp1(dep_wb500,temp_wb500,dep_wb500_grid,'linear','extrap');
temp_wb500_grid(1,:) = temp_wb500_grid(2,:);
salt_wb500_grid = interp1(dep_wb500,salt_wb500,dep_wb500_grid,'linear','extrap');
salt_wb500_grid(1,:) = salt_wb500_grid(2,:);

temp_k4_grid = interp1(dep_k4,temp_k4,dep_k4_grid,'linear','extrap');
temp_k4_grid(1,:) = temp_k4_grid(2,:);
salt_k4_grid = interp1(dep_k4,salt_k4,dep_k4_grid,'linear','extrap');
salt_k4_grid(1,:) = salt_k4_grid(2,:);

temp_eb1200_grid = interp1(dep_eb1200,temp_eb1200,dep_eb1200_grid,'linear','extrap');
temp_eb1200_grid(1,:) = temp_eb1200_grid(2,:);
salt_eb1200_grid = interp1(dep_eb1200,salt_eb1200,dep_eb1200_grid,'linear','extrap');
salt_eb1200_grid(1,:) = salt_eb1200_grid(2,:);

% computing specific volume anomaly from interpolated data
svan_wb300_grid = sw_svan(salt_wb300_grid, temp_wb300_grid, pres_wb300_grid');
svan_wb500_grid = sw_svan(salt_wb500_grid, temp_wb500_grid, pres_wb500_grid');
svan_k4_grid = sw_svan(salt_k4_grid, temp_k4_grid, pres_k4_grid');
svan_eb1200_grid = sw_svan(salt_eb1200_grid, temp_eb1200_grid, pres_eb1200_grid');

% computing geopotential anomaly from interpolated data
gpan_wb300_grid = sw_gpan(salt_wb300_grid, temp_wb300_grid, pres_wb300_grid');
gpan_wb500_grid = sw_gpan(salt_wb500_grid, temp_wb500_grid, pres_wb500_grid');
gpan_k4_grid = sw_gpan(salt_k4_grid, temp_k4_grid, pres_k4_grid');
gpan_eb1200_grid = sw_gpan(salt_eb1200_grid, temp_eb1200_grid, pres_eb1200_grid');

% computing bottom pressure from interpolated data
rho_wb300_grid = sw_dens(salt_wb300_grid,temp_wb300_grid, pres_wb300_grid');
bp_dyn_wb300 = g .* trapz(dep_wb300_grid,rho_wb300_grid);
bp_ssh_wb300 = rho_wb300_grid(1,:) .* g .* ssh_wb300;
rho_wb500_grid = sw_dens(salt_wb500_grid,temp_wb500_grid, pres_wb500_grid');
bp_dyn_wb500 = g .* trapz(dep_wb500_grid,rho_wb500_grid);
bp_ssh_wb500 = rho_wb500_grid(1,:) .* g .* ssh_wb500;
rho_k4_grid = sw_dens(salt_k4_grid,temp_k4_grid, pres_k4_grid');
bp_dyn_k4 = g .* trapz(dep_k4_grid,rho_k4_grid);
bp_ssh_k4 = rho_k4_grid(1,:) .* g .* ssh_k4;
rho_eb1200_grid = sw_dens(salt_eb1200_grid,temp_eb1200_grid, pres_eb1200_grid');
bp_dyn_eb1200 = g .* trapz(dep_eb1200_grid,rho_eb1200_grid);
bp_ssh_eb1200 = rho_eb1200_grid(1,:) .* g .* ssh_eb1200;

% computing Fofonoff potential from interpolated data (for WBb5 and EBb4)
% 1. compute pressure
bp3_wb300(1,:) = bp_ssh_wb300;
bp3_wb500(1,:) = bp_ssh_wb500;
bp3_eb1200(1,:) = bp_ssh_eb1200;
bp5_wb500(1,:) = bp_ssh_wb500;
bp5_eb1200(1,:) = bp_ssh_eb1200;
for iz = 2:length(dep_wb500_grid) % until 440m
    bp5_wb500(iz,:) = g .* trapz(dep_wb500_grid(1:iz),rho_wb500_grid(1:iz)) + bp_ssh_wb500;
    bp5_eb1200(iz,:) = g .* trapz(dep_eb1200_grid(1:iz),rho_eb1200_grid(1:iz)) + bp_ssh_eb1200;
end
for iz = 2:length(dep_wb300_grid) % until 270m
    bp3_wb300(iz,:) = g .* trapz(dep_wb300_grid(1:iz),rho_wb300_grid(1:iz)) + bp_ssh_wb300;
    bp3_wb500(iz,:) = g .* trapz(dep_wb500_grid(1:iz),rho_wb500_grid(1:iz)) + bp_ssh_wb500;
    bp3_eb1200(iz,:) = g .* trapz(dep_eb1200_grid(1:iz),rho_eb1200_grid(1:iz)) + bp_ssh_eb1200;
end
% 2. integrate pressure to Foffenoff potential
fof3_wb300 = trapz(dep_wb300_grid,bp3_wb300);
fof3_wb500 = trapz(dep_wb300_grid,bp3_wb500);
fof3_eb1200 = trapz(dep_wb300_grid,bp3_eb1200);
fof5_wb500 = trapz(dep_wb500_grid,bp5_wb500);
fof5_eb1200 = trapz(dep_wb500_grid,bp5_eb1200);

% cubic fit for every depth level and location --> generation of GEM table
% and application to reconstruct variables
% -- WB 300m
tau_fit_wb300 = min(tau_wb300):0.0001:max(tau_wb300)+0.0001;
temp_wb300_fit = ones(length(dep_wb300_grid),length(tau_fit_wb300)).*NaN;
salt_wb300_fit = ones(length(dep_wb300_grid),length(tau_fit_wb300)).*NaN;
svan_wb300_fit = ones(length(dep_wb300_grid),length(tau_fit_wb300)).*NaN;
temp_wb300_reconstr = ones(length(dep_wb300_grid),length(tau_wb300)).*NaN;
salt_wb300_reconstr = ones(length(dep_wb300_grid),length(tau_wb300)).*NaN;
svan_wb300_reconstr = ones(length(dep_wb300_grid),length(tau_wb300)).*NaN;
for iz = 1:length(dep_wb300_grid)
    [p,S] = polyfit(tau_wb300,temp_wb300_grid(iz,:),3); % cubic polynomial
    [temp_fit, temp_fit_error] = polyval(p,tau_fit_wb300,S);
    temp_wb300_fit(iz,:) = temp_fit;
    temp_wb300_reconstr(iz,:) = polyval(p,tau_wb300);
    [p,S] = polyfit(tau_wb300,salt_wb300_grid(iz,:),3);
    [salt_fit, salt_fit_error] = polyval(p,tau_fit_wb300,S);
    salt_wb300_fit(iz,:) = salt_fit;
    salt_wb300_reconstr(iz,:) = polyval(p,tau_wb300);
    [p,S] = polyfit(tau_wb300,svan_wb300_grid(iz,:),3);
    [svan_fit, svan_fit_error] = polyval(p,tau_fit_wb300,S);
    svan_wb300_fit(iz,:) = svan_fit;
    svan_wb300_reconstr(iz,:) = polyval(p,tau_wb300);
end
[p,S] = polyfit(tau_wb300,gpan_wb300(end,:),3);
[gpan_wb300_fit, gpan_wb300_fit_error] = polyval(p,tau_fit_wb300,S);
gpan_wb300_reconstr = polyval(p,tau_wb300);

% -- WB 500m
tau_fit_wb500 = min(tau_wb500):0.0001:max(tau_wb500)+0.0001;
temp_wb500_fit = ones(length(dep_wb500_grid),length(tau_fit_wb500)).*NaN;
salt_wb500_fit = ones(length(dep_wb500_grid),length(tau_fit_wb500)).*NaN;
svan_wb500_fit = ones(length(dep_wb500_grid),length(tau_fit_wb500)).*NaN;
temp_wb500_reconstr = ones(length(dep_wb500_grid),length(tau_wb500)).*NaN;
salt_wb500_reconstr = ones(length(dep_wb500_grid),length(tau_wb500)).*NaN;
svan_wb500_reconstr = ones(length(dep_wb500_grid),length(tau_wb500)).*NaN;
for iz = 1:length(dep_wb500_grid)
    [p,S] = polyfit(tau_wb500,temp_wb500_grid(iz,:),3);
    [temp_fit, temp_fit_error] = polyval(p,tau_fit_wb500,S);
    temp_wb500_fit(iz,:) = temp_fit;
    temp_wb500_reconstr(iz,:) = polyval(p,tau_wb500);
    [p,S] = polyfit(tau_wb500,salt_wb500_grid(iz,:),3);
    [salt_fit, salt_fit_error] = polyval(p,tau_fit_wb500,S);
    salt_wb500_fit(iz,:) = salt_fit;
    salt_wb500_reconstr(iz,:) = polyval(p,tau_wb500);
    [p,S] = polyfit(tau_wb500,svan_wb500_grid(iz,:),3);
    [svan_fit, svan_fit_error] = polyval(p,tau_fit_wb500,S);
    svan_wb500_fit(iz,:) = svan_fit;
    svan_wb500_reconstr(iz,:) = polyval(p,tau_wb500);
end
[p,S] = polyfit(tau_wb500,gpan_wb500(end,:),3);
[gpan_wb500_fit, gpan_wb500_fit_error] = polyval(p,tau_fit_wb500,S);
gpan_wb500_reconstr = polyval(p,tau_wb500);

% -- K4
tau_fit_k4 = min(tau_k4):0.0001:max(tau_k4)+0.0001;
temp_k4_fit = ones(length(dep_k4_grid),length(tau_fit_k4)).*NaN;
salt_k4_fit = ones(length(dep_k4_grid),length(tau_fit_k4)).*NaN;
svan_k4_fit = ones(length(dep_k4_grid),length(tau_fit_k4)).*NaN;
temp_k4_reconstr = ones(length(dep_k4_grid),length(tau_k4)).*NaN;
salt_k4_reconstr = ones(length(dep_k4_grid),length(tau_k4)).*NaN;
svan_k4_reconstr = ones(length(dep_k4_grid),length(tau_k4)).*NaN;
for iz = 1:length(dep_k4_grid)
    [p,S] = polyfit(tau_k4,temp_k4_grid(iz,:),3);
    [temp_fit, temp_fit_error] = polyval(p,tau_fit_k4,S);
    temp_k4_fit(iz,:) = temp_fit;
    temp_k4_reconstr(iz,:) = polyval(p,tau_k4);
    [p,S] = polyfit(tau_k4,salt_k4_grid(iz,:),3);
    [salt_fit, salt_fit_error] = polyval(p,tau_fit_k4,S);
    salt_k4_fit(iz,:) = salt_fit;
    salt_k4_reconstr(iz,:) = polyval(p,tau_k4);
    [p,S] = polyfit(tau_k4,svan_k4_grid(iz,:),3);
    [svan_fit, svan_fit_error] = polyval(p,tau_fit_k4,S);
    svan_k4_fit(iz,:) = svan_fit;
    svan_k4_reconstr(iz,:) = polyval(p,tau_k4);
end
[p,S] = polyfit(tau_k4,gpan_k4(end,:),3);
[gpan_k4_fit, gpan_k4_fit_error] = polyval(p,tau_fit_k4,S);
gpan_k4_reconstr = polyval(p,tau_k4);

% -- EB 1200m
tau_fit_eb1200 = min(tau_eb1200):0.0001:max(tau_eb1200)+0.0001;
temp_eb1200_fit = ones(length(dep_eb1200_grid),length(tau_fit_eb1200)).*NaN;
salt_eb1200_fit = ones(length(dep_eb1200_grid),length(tau_fit_eb1200)).*NaN;
svan_eb1200_fit = ones(length(dep_eb1200_grid),length(tau_fit_eb1200)).*NaN;
temp_eb1200_reconstr = ones(length(dep_eb1200_grid),length(tau_eb1200)).*NaN;
salt_eb1200_reconstr = ones(length(dep_eb1200_grid),length(tau_eb1200)).*NaN;
svan_eb1200_reconstr = ones(length(dep_eb1200_grid),length(tau_eb1200)).*NaN;
for iz = 1:length(dep_eb1200_grid)
    [p,S] = polyfit(tau_eb1200,temp_eb1200_grid(iz,:),3);
    [temp_fit, temp_fit_error] = polyval(p,tau_fit_eb1200,S);
    temp_eb1200_fit(iz,:) = temp_fit;
    temp_eb1200_reconstr(iz,:) = polyval(p,tau_eb1200);
    [p,S] = polyfit(tau_eb1200,salt_eb1200_grid(iz,:),3);
    [salt_fit, salt_fit_error] = polyval(p,tau_fit_eb1200,S);
    salt_eb1200_fit(iz,:) = salt_fit;
    salt_eb1200_reconstr(iz,:) = polyval(p,tau_eb1200);
    [p,S] = polyfit(tau_eb1200,svan_eb1200_grid(iz,:),3);
    [svan_fit, svan_fit_error] = polyval(p,tau_fit_eb1200,S);
    svan_eb1200_fit(iz,:) = svan_fit;
    svan_eb1200_reconstr(iz,:) = polyval(p,tau_eb1200);
end
[p,S] = polyfit(tau_eb1200,gpan_eb1200(end,:),3);
[gpan_eb1200_fit, gpan_eb1200_fit_error] = polyval(p,tau_fit_eb1200,S);
gpan_eb1200_reconstr = polyval(p,tau_eb1200);

% --- save variables in interim folder, then download to plot (then script plot_MooringSetup) ---
save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/ies_test_tauts_' nmid test_suffix '.mat'],'time', ...
    'dep_eb1200','tau_eb1200','ptemp_eb1200','salt_eb1200', 'svan_eb1200', 'rho_eb1200',...
    'dep_wb300','tau_wb300','ptemp_wb300','salt_wb300', 'svan_wb300', 'rho_wb300', ...
    'dep_wb500','tau_wb500','ptemp_wb500','salt_wb500', 'svan_wb500', 'rho_wb500', ...
    'dep_k4','tau_k4','ptemp_k4','salt_k4', 'svan_k4', 'rho_k4','-v7.3')

save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/ies_test_reconstr_' nmid test_suffix '.mat'],'time', ...
    'dep_wb300_grid','temp_wb300_grid','temp_wb300_reconstr','salt_wb300_grid','salt_wb300_reconstr','svan_wb300_grid','svan_wb300_reconstr', ...
    'dep_wb500_grid','temp_wb500_grid','temp_wb500_reconstr','salt_wb500_grid','salt_wb500_reconstr','svan_wb500_grid','svan_wb500_reconstr', ...
    'dep_k4_grid','temp_k4_grid','temp_k4_reconstr','salt_k4_grid','salt_k4_reconstr','svan_k4_grid','svan_k4_reconstr', ...
    'dep_eb1200_grid','temp_eb1200_grid','temp_eb1200_reconstr','salt_eb1200_grid','salt_eb1200_reconstr','svan_eb1200_grid','svan_eb1200_reconstr','-v7.3')

save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/ies_test_gpan_' nmid test_suffix '.mat'],'time', ...
    'tau_wb300','tau_fit_wb300','gpan_wb300','gpan_wb300_fit','gpan_wb300_reconstr', ...
    'tau_wb500','tau_fit_wb500','gpan_wb500','gpan_wb500_fit','gpan_wb500_reconstr', ...
    'tau_k4','tau_fit_k4','gpan_k4','gpan_k4_fit','gpan_k4_reconstr', ...
    'tau_eb1200','tau_fit_eb1200','gpan_eb1200','gpan_eb1200_fit','gpan_eb1200_reconstr','-v7.3')

save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/ies_test_bp_' nmid test_suffix '.mat'],'time', ...
    'tau_wb300','bp_dyn_wb300','bp_ssh_wb300', ...
    'tau_wb500','bp_dyn_wb500','bp_ssh_wb500', ...
    'tau_k4','bp_dyn_k4','bp_ssh_k4', ...
    'tau_eb1200','bp_dyn_eb1200','bp_ssh_eb1200','-v7.3')

save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/ies_test_foff_' nmid test_suffix '.mat'],'time', ...
        'tau_wb300','fof3_wb300', ...
        'tau_wb500','fof3_wb500', 'fof5_wb500', ...
        'tau_eb1200','fof3_eb1200','fof5_eb1200','-v7.3')
