%% INALT/VIKING model true AMOC transport - cut
% to run on nesh
% gives total transport, geostrophic transport and Ekman transport

% choose INA (1) or VIK (2)
mod_option = 2;
if mod_option == 1
    num_yy = 40;
    nmid = 'INA';
    mmask = '/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/1_INALT20.L46-KFS119_mesh_mask_11S_box.nc';
elseif mod_option == 2
    num_yy = 44;
    nmid = 'VIK';
    mmask = '/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/1_VIKING20X.L46-KFS003_mesh_mask_11S_box.nc';
end

addpath('/gxfs_home/geomar/smomw628/matlab_toolbox/seawater_ver3_3.1/')

% choose latitude to compute for
lat_option = 1;
if lat_option == 1
    lat_use = -11;
    nlid = '1100';
elseif lat_option == 2
    lat_use = -10.75;
    nlid = '1075';
elseif lat_option == 3
    lat_use = -10.5;
    nlid = '1050';
end

%% total velocities and geostrophic velocities

if mod_option == 1
	files1 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*V_11S_box.nc');
    files2 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*V_11S_box.vgeostrophy.nc');
elseif mod_option == 2
	files1 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*V_11S_box.nc');
    files2 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*V_11S_box.vgeostrophy.nc');
end

lat = double(ncread(mmask,'gphiv'));
lat = lat(1,:);
lonV = double(ncread(mmask,'glamv'));
lonV = lonV(:,1);
x_spacing = double(ncread(mmask,'e1v')); % in m
z_spacing = double(ncread(mmask,'e3v_0')); % in m
dep = double(ncread(mmask,'gdepv')); % in m, including partial cells
qc = double(ncread(mmask,'vmask')); % 0 land, 1 water
qc(qc==0) = NaN;

% choose latitude (3 options to see whether there are differences)
ilat = find(lat(1,:) > lat_use,1,'first');
   
dep = squeeze(dep(:,ilat,:));
x_spacing = x_spacing(1,ilat);
z_spacing = squeeze(z_spacing(:,ilat,:));

countt = 1;
for i = 1:num_yy % loop over years

    % load all meridional velocities at one specified latitude (make 3 options)
    % and merge them to one file
    tim = double(ncread([files1(i).folder '/' files1(i).name],'time_counter'));
    time(countt:countt-1+length(tim)) = datenum(1900,1,1,0,0,tim);
    
    v_use = double(ncread([files1(i).folder '/' files1(i).name],'vomecrty')); % in m/s
    vg_use = double(ncread([files2(i).folder '/' files2(i).name],'vomecrty')); % in m/s

    % applying land-sea mask
    v_use = v_use .* qc;
    vg_use = vg_use .* qc;

    v(:,:,countt:countt-1+length(tim)) = squeeze(v_use(:,ilat,:,:));
    vg(:,:,countt:countt-1+length(tim)) = squeeze(vg_use(:,ilat,:,:));
    countt = countt + length(tim);
end

% --- deriving quantities to save and plot

% interpolation on higher resolution vertical grid
dep_grid = (0:10:dep(end))';
v_grid = ones(length(lonV),length(dep_grid),length(time)).*NaN;
vg_grid = ones(length(lonV),length(dep_grid),length(time)).*NaN;
for i = 1:length(lonV)
    if sum(~isnan(v(i,:,1)))>=2 % interpolation only possible if more than two depth contain data
        idep = ~isnan(v(i,:,1)); % determine maximum interpolation depth
        idep_maxnew = sum(idep .* z_spacing(i,:));
        v_grid(i,dep_grid<idep_maxnew,:) = interp1(dep(i,idep),squeeze(v(i,idep,:)),dep_grid(dep_grid<idep_maxnew),'spline');
    end
    if sum(~isnan(vg(i,:,1)))>=2
        idep = ~isnan(vg(i,:,1));
        idep_maxnew = sum(idep .* z_spacing(i,:));
        vg_grid(i,dep_grid<idep_maxnew,:) = interp1(dep(i,idep),squeeze(vg(i,idep,:)),dep_grid(dep_grid<idep_maxnew),'spline');
    end
end

% one mean over time, depth vs lon (-- lon, dep, v_mean)
v_mean = mean(v_grid,3);
vg_mean = mean(vg_grid,3);
% one inegral over lon, time vs depth (-- time, dep, vdx)
T_basin = squeeze(sum(v_grid,1,'omitnan')) .* x_spacing;
Tg_basin = squeeze(sum(vg_grid,1,'omitnan')) .* x_spacing;
% Tg also for different parts of the basin (longitudes choosen like in MooringSetup)
[~,ilon_k4] = min(abs(lonV+34.99));
[~,ilon_a1] = min(abs(lonV-13.00)); ilon_a1 = ilon_a1 - 1;
Tg_wb = squeeze(sum(vg_grid(1:ilon_k4,:,:),1,'omitnan')) .* x_spacing;
Tg_shelf_eb = squeeze(sum(vg_grid(ilon_a1:end,:,:),1,'omitnan')) .* x_spacing;
% overturning streamfunction
psi = cumtrapz(dep_grid,T_basin,1);
psig = cumtrapz(dep_grid,Tg_basin,1);

% --- save variables in interim folder, then download to plot (then script plot_IdealSetup)
save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/v_' nlid 'S_' nmid '.mat'],'time','lonV','dep_grid','v_mean','T_basin','psi','-v7.3')
save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/vg_' nlid 'S_' nmid '.mat'],'time','lonV','dep_grid','vg_mean','Tg_basin','Tg_wb','Tg_shelf_eb','psig','-v7.3')

%% Ekman transport

if mod_option == 1
        files1 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*U_11S_box.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*T_11S_box.nc');
elseif mod_option == 2
        files1 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*U_11S_box.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*T_11S_box.nc');
end

lat = double(ncread(mmask,'gphiu'));
lat = lat(1,:); % lat for grid U and T are the same
lonU = double(ncread(mmask,'glamu'));
lonU = lonU(:,1);
x_spacing = double(ncread(mmask,'e1u')); % in m
qcu = double(ncread(mmask,'umask')); % 0 land, 1 water
qcu(qcu==0) = NaN;
qct = double(ncread(mmask,'tmask')); % 0 land, 1 water
qct(qct==0) = NaN;

ilat = find(lat(1,:) > lat_use,1,'first');
x_spacing = x_spacing(1,ilat);

countt = 1;
for i = 1:num_yy % loop over years
    tim = double(ncread([files1(i).folder '/' files1(i).name],'time_centered'));
    time(countt:countt-1+length(tim)) = datenum(1900,1,1,0,0,tim);
    
    taux_use = double(ncread([files1(i).folder '/' files1(i).name],'sozotaux')); % in N/m2
    
    % density for determining rho0 as mean surface density
    ptemp_use = double(ncread([files2(i).folder '/' files2(i).name],'votemper')); % in °C, potential temperature
    salt_use = double(ncread([files2(i).folder '/' files2(i).name],'vosaline')); % in 0.001, practical salinity

    % applying land-sea mask
    taux_use = taux_use .* qcu(:,:,1);
    ptemp_use = ptemp_use .* qct;
    salt_use = salt_use .* qct;

    taux(:,countt:countt-1+length(tim)) = squeeze(taux_use(:,ilat,:));
    ptemp(:,countt:countt-1+length(tim)) = squeeze(ptemp_use(:,ilat,1,:));
    salt(:,countt:countt-1+length(tim)) = squeeze(salt_use(:,ilat,1,:));
    countt = countt + length(tim);
end


temp = sw_temp(salt,ptemp,0,0);
rho = sw_dens(salt,temp,0);

% --- deriving quantities to save and plot
f = sw_f(lat_use);
rho0 = mean(rho,'all','omitnan'); % --> 1024.2 kg/m^3
T_Ek = - 1/(rho0*f) .* sum(taux,1,'omitnan') .* x_spacing;

% --- save variables in interim folder, then download to plot (then script plot_IdealSetup)
save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/TEk_' nlid 'S_' nmid '.mat'],'time','T_Ek','-v7.3')
