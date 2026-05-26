function f_cut_MooringSetup(mod_option,wbr_option,ebr_option,vst_option,eof_option,sur_option,idr_option,test_suffix)

% function to cut INALT/VIKING Mooring AMOCg transport and to be run on nesh

% input:
% mod_option: choose model INALT20 (1) or VIKING20X (2)
% wbr_option: choose WB4 setup as ideal (0), initial (1), current (2), or enhanced (3)
% ebr_option: choose EB1 setup as
%       ideal (0), initial (1), current (2), or enhanced (3), or ideal averaged to climatological year (4) 
% vst_option: choose vertical structure as determined by
%       stepwise integration (1, Johns 2005), T'/S' interpolation (2, Williams 2015), EOF regression (3), or (0) if wbr and ebr ideal
% eof_option: choose number EOFs as whether the first 1 (1) or 2 (2) or 3 (3) EOFs are used, choose (0) if vst~=3
% sur_option: choose SSS like
%       normal measurement (0), climatological year (1), anomaly equals uppermost instrument - Williams 2015 (2),
%       seasonal extrapolation - McCarthy 2015 (3), climatological value - McCarthy 2015 (4)
% idr_option: choose integration approach
%       from surface (1) or from level of nomo (2), EB1 from EBb4 (BP1200) - highpass filtered (3),
%       EB1 from EBb4, WB4 from WBb6 - highpass filtered (4), same as 3 but not filtered (5)
% test_suffix: adds a suffix in the saved file name in code testing cases, string, default ''



% --- no modifications needed below ---

% choose model
if mod_option == 1
    num_yy = 40;
    nmid = 'INA';
elseif mod_option == 2
    num_yy = 44;
    nmid = 'VIK';
end

addpath('/gxfs_home/geomar/smomw628/matlab_toolbox/seawater_ver3_3.1/')

% --- Subsampling mooring data ---
if mod_option == 1
        files1 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*T_11S_box.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*T_11S_box.ncbotpressure_detrssh.nc');
        files3 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*V_11S_box.nc');
        mmask = '/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/1_INALT20.L46-KFS119_mesh_mask_11S_box.nc';
elseif mod_option == 2
        files1 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*T_11S_box.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*T_11S_box.ncbotpressure_detrssh.nc');
        files3 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*V_11S_box.nc');
        mmask = '/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/1_VIKING20X.L46-KFS003_mesh_mask_11S_box.nc';
end

lat = double(ncread(mmask,'gphit'));
lat = lat(1,:);
lon = double(ncread(mmask,'glamt'));
lon = lon(:,1);
dep = double(ncread(mmask,'gdept_0')); % in m, including partial cells
qc = double(ncread(mmask,'tmask')); % 0 land, 1 water
qc(qc==0) = NaN;

% choose positions close to moorings
%[~,ilat_k1] = min(abs(lat+10.27));
%[~,ilon_k1] = min(abs(lon+35.86));
%[~,ilat_k2] = min(abs(lat+10.38));
%[~,ilon_k2] = min(abs(lon+35.68));
%[~,ilat_k3] = min(abs(lat+10.61));
%[~,ilon_k3] = min(abs(lon+35.39));
[~,ilat_k4] = min(abs(lat+10.94));
[~,ilon_k4] = min(abs(lon+34.99));
[~,ilat_a1] = min(abs(lat+10.83));
[~,ilon_a1] = min(abs(lon-13.00));
ilon_a1 = ilon_a1 - 1; % so that the mooring reaches a depth of 1200m

if idr_option == 4 % find ilat_k4 seperately, ilon_k4 is the same on T and V grid
    latV = double(ncread(mmask,'gphiv'));
    latV = latV(1,:);
    [~,ilat_k4_V] = min(abs(latV+10.94));
    qcv = double(ncread(mmask,'vmask')); % 0 land, 1 water
    qcv(qcv==0) = NaN;
end

% get latitudes and depths for the used moorings
moor_lat = [lat(ilat_k4) lat(ilat_a1)];
lat_mean = mean(moor_lat);
moor_dep(1,:) = dep(ilon_k4,ilat_k4,1:29);
moor_dep(2,:) = dep(ilon_a1,ilat_a1,1:29);
moor_dep(:,1) = 0; % move uppermost value to surface (feasible as there should be no gradient within the mixed layer)

if idr_option == 4
    x_spacing = double(ncread(mmask,'e1v')); % in m
    x_spacing = x_spacing(1,ilat_k4);
end

% get temperature, salinity and SSH at mooring locations
countt = 1;
for i = 1:num_yy % loop over years
    tim = double(ncread([files1(i).folder '/' files1(i).name],'time_counter'));
    time(countt:countt-1+length(tim)) = datenum(1900,1,1,0,0,tim);
    
    ptemp_use = double(ncread([files1(i).folder '/' files1(i).name],'votemper')); % in °C, potential temperature
    salt_use = double(ncread([files1(i).folder '/' files1(i).name],'vosaline')); % in 0.001, practical salinity
    ssh_use = double(ncread([files2(i).folder '/' files2(i).name],'sossheig')); % in m

    % applying land-sea mask
    ptemp_use = ptemp_use .* qc;
    salt_use = salt_use .* qc;
    ssh_use = ssh_use .* qc(:,:,1);
    
    % temperature and salilnity up to 1850 m (bin 29) and SSH
    ptemp_dly(1,:,countt:countt-1+length(tim)) = squeeze(ptemp_use(ilon_k4,ilat_k4,1:29,:));
    salt_dly(1,:,countt:countt-1+length(tim)) = squeeze(salt_use(ilon_k4,ilat_k4,1:29,:));
    ssh_dly(1,countt:countt-1+length(tim)) = squeeze(ssh_use(ilon_k4,ilat_k4,:));

    ptemp_dly(2,:,countt:countt-1+length(tim)) = squeeze(ptemp_use(ilon_a1,ilat_a1,1:29,:));
    salt_dly(2,:,countt:countt-1+length(tim)) = squeeze(salt_use(ilon_a1,ilat_a1,1:29,:));
    ssh_dly(2,countt:countt-1+length(tim)) = squeeze(ssh_use(ilon_a1,ilat_a1,:));

    if idr_option == 4
        v_use = double(ncread([files3(i).folder '/' files3(i).name],'vomecrty')); % in m/s
        v_use = v_use .* qcv;
        v_WB(:,countt:countt-1+length(tim)) = squeeze(v_use(1:ilon_k4,ilat_k4_V,19,:)); % bin 19 is at 452 m depth and thus roughly matches BP500
    end

    countt = countt + length(tim);
end

clear ptemp_use salt_use ssh_use qc lon lat dep tim countt i mmask files1 files2
clear ilat_k4 ilon_k4 ilat_a1 ilon_a1


% --- deriving quantities to save and plot ---

% averaging to monthly data
[yy,mm,~] = datevec(time);
time_mly = unique(datenum(yy,mm,15));
[~,time_mly_mm,~] = datevec(time_mly);
ptemp = ones(2,size(moor_dep,2),length(time_mly)).*NaN;
salt = ones(2,size(moor_dep,2),length(time_mly)).*NaN;
ssh = ones(2,length(time_mly)).*NaN;
ssh_mly_eb = ones(1,length(time_mly)).*NaN;
if idr_option == 4
    v_WB_mly = ones(size(v_WB,1),length(time_mly)).*NaN;
end
for i = 1:length(time_mly)
    % averaging to monthly data for western boundary
    ptemp(1,:,i) = mean(ptemp_dly(1,:,datenum(yy,mm,15)==time_mly(i)),3);
    salt(1,:,i) = mean(salt_dly(1,:,datenum(yy,mm,15)==time_mly(i)),3);
    ssh(1,i) = mean(ssh_dly(1,datenum(yy,mm,15)==time_mly(i)),2);
    % ssh_mly_eb(i) = mean(ssh_dly(2,datenum(yy,mm,15)==time_mly(i)),2);
    if ebr_option < 4
        % averaging to monthly data for eastern boundary
        ptemp(2,:,i) = mean(ptemp_dly(2,:,datenum(yy,mm,15)==time_mly(i)),3);
        salt(2,:,i) = mean(salt_dly(2,:,datenum(yy,mm,15)==time_mly(i)),3);
        ssh(2,i) = mean(ssh_dly(2,datenum(yy,mm,15)==time_mly(i)),2);
    elseif ebr_option == 4
        % climatological year for eastern boundary
        ptemp(2,:,i) = mean(ptemp_dly(2,:,mm==time_mly_mm(i)),3);
        salt(2,:,i) = mean(salt_dly(2,:,mm==time_mly_mm(i)),3);
        ssh(2,i) = mean(ssh_dly(2,mm==time_mly_mm(i)),2);
        ptemp_eb_orig(:,i) = mean(ptemp_dly(2,:,datenum(yy,mm,15)==time_mly(i)),3);
        salt_eb_orig(:,i) = mean(salt_dly(2,:,datenum(yy,mm,15)==time_mly(i)),3);
    end
    if idr_option == 4
        v_WB_mly(:,i) = mean(v_WB(:,datenum(yy,mm,15)==time_mly(i)),2);
    end
end
if ebr_option == 4 % adding SSH signal
    % version 0: also SSH seasonal cycle (nothing added)

    % version 1: add linear trend (esp. sea level rise)
%     ssh_trend_eb = polyfit(time_mly,ssh_mly_eb,1);
%     ssh(2,:) = ssh(2,:) + polyval(ssh_trend_eb,time_mly) - mean(ssh(2,:),'all');
    
    % version 2: add annual (running) mean SSH
%     ssh_ann_eb = movmean(ssh_mly_eb,12);
%     ssh(2,:) = ssh(2,:) + ssh_ann_eb - mean(ssh(2,:),'all');
% 
    % version 3: add lowpass filtered SSH - beginning and end should still be capped
%     Fs = 1/(30*24*3600); % Sampling frequency in Hz (monthly)
%     Fc_lp = 1/(5*365*24*3600); % Cutoff frequency in Hz
%     [b, a] = butter(1, Fc_lp/(Fs/2), 'low'); % 1st-order low-pass filter
%     ssh_lp_eb = filtfilt(b,a,ssh_mly_eb);
%     ssh(2,:) = ssh(2,:) + ssh_lp_eb - mean(ssh(2,:),'all');
end
clear i mm yy ptemp_dly salt_dly ssh_dly ssh_mly_eb

% computing density from T,S following EOS80
moor_pres = sw_pres(moor_dep,moor_lat');
rho = ones(2,size(moor_dep,2),length(time_mly)).*NaN;
for i = 1:2
    temp = sw_temp(squeeze(salt(i,:,:)),squeeze(ptemp(i,:,:)),moor_pres(i,:)',0);
    rho(i,:,:) = sw_dens(squeeze(salt(i,:,:)),temp,moor_pres(i,:)');
end

% subsampling according to mooring setup
if wbr_option == 1
    wb_bins = [1 20 21 29]; % WB4: surface, 500m, 650m, 1900m
elseif wbr_option == 2
    wb_bins = [1 10 14 16 18 20 21 29]; % WB4ext: surface, 100m, 200m, 300m, 400m, 500m, 650m, 1900m
elseif wbr_option == 3
    wb_bins = [1 7 10 14 16 18 20 21 23 25 29]; % WB4ext+: surface, 50m, 100m, 200m, 300m, 400m, 500m, 650m, 850m, 1100m, 1900m
end

if ebr_option == 1
    eb_bins = [1 17 20 22 24 26]; % EB1: surface, 300m, 500m, 700m, 950m, 1200m
elseif ebr_option == 2
    eb_bins = [1 13 17 20 22 24 26]; % EB1ext: surface, 150m, 300m, 500m, 700m, 950m, 1200m
elseif ebr_option == 3
    eb_bins = [1 4 10 14 16 18 20 22 24 26]; % EB1ext+: surface, 20m, 100m, 200m, 300m, 400m, 500m, 700m, 950m, 1200m
    %eb_bins = [1 4 10 14 16 18 20 21 22 23 24 25 26]; % Angola ext+2: surface, 20m, 100m, 200m, 300m, 400m, 500m, 600m, 700m, 850m, 950m, 1100m, 1200m
    %eb_bins = [1 13 17 20 22 24 26]; % Angola plan 1: surface, 150m, 300m, 500m, 700m, 950m, 1200m
    %eb_bins = [1 17 18 20 22 24 26]; % Angola plan 2: surface, 300m, 400m, 500m, 700m, 950m, 1200m
end

% --- computing and applying the vertical structure ---
z_spacing = 10;
dep_grid = 0:z_spacing:1250; % 1250 as this is the model EB1 waterdepth
pres_grid = sw_pres(dep_grid,lat_mean);

% vertical interpolation for reference
wb_rho_grid_orig = interp1(moor_dep(1,1:26),squeeze(rho(1,1:26,:)),dep_grid,'spline');
if ebr_option ~= 4
    eb_rho_grid_orig = interp1(moor_dep(2,1:26),squeeze(rho(2,1:26,:)),dep_grid,'spline');
elseif ebr_option == 4
    temp_eb_orig = sw_temp(salt_eb_orig,ptemp_eb_orig,moor_pres(2,:)',0);
    rho_eb_orig = sw_dens(salt_eb_orig,temp_eb_orig,moor_pres(2,:)');
    eb_rho_grid_orig = interp1(moor_dep(2,1:26),rho_eb_orig(1:26,:),dep_grid,'spline');
end

% vertical interpolation for ideal moorings
if wbr_option == 0
    wb_rho_grid = wb_rho_grid_orig;
end
if ebr_option == 0
    eb_rho_grid = eb_rho_grid_orig;
elseif ebr_option == 4
    eb_rho_grid = interp1(moor_dep(2,1:26),squeeze(rho(2,1:26,:)),dep_grid,'spline');
end

salt_orig = salt;
if sur_option == 1 % compute climatological year for first salt bin
    for i = 1:length(time_mly)
        salt(:,1,i) = squeeze(mean(salt(:,1,time_mly_mm==time_mly_mm(i)),3));
    end
    for i = 1:2
        temp = sw_temp(squeeze(salt(i,:,:)),squeeze(ptemp(i,:,:)),moor_pres(i,:)',0);
        rho(i,:,:) = sw_dens(squeeze(salt(i,:,:)),temp,moor_pres(i,:)');
    end
elseif sur_option == 4 % compute mean for first salt bin
    salt(1,1,:) = mean(salt(1,1,:),3);
    salt(2,1,:) = mean(salt(2,1,:),3);
end

if vst_option == 3
    % computing the first 3 EOFs for density profiles
    EOFs_all = ones(2,size(moor_dep,2),3).*NaN; % moorings x depth_grid x EOFs
    for i = 1:2
        if i == 1
            rho4EOF = squeeze(rho(i,:,:))' - mean(squeeze(rho(i,:,:)),2)';
        else
            rho4EOF = squeeze(rho(i,1:26,:))' - mean(squeeze(rho(i,1:26,:)),2)';
        end
         
        % 1.0 calculate covariance matrix
        C = cov(rho4EOF);
        
        % 2.0 do eigenanalysis
        % columns of EOFs are eigenvectors, D diagonal matrix of eigenvalues
        [EOFs,D] = eigs(C);
        
        % 3.0 identify the first modes and calculate eyplained variance
        lambda = diag(D);
        exp_var = lambda ./ sum(lambda,'omitnan');
        disp('Explained variance [%]')
        disp(exp_var.*100)
        
        if i == 1
            EOFs_all(i,:,:) = EOFs(:,1:3);
        else
            EOFs_all(i,1:26,:) = EOFs(:,1:3);
        end
    end

    % regressing EOFs on subsampled data
    nb_eof_use = 1:eof_option;

    % WB
    if wbr_option > 0
        EOF_use = squeeze(EOFs_all(1,wb_bins,nb_eof_use));
        if eof_option == 1
            EOF_use = EOF_use';
        end
        obs_use = squeeze(rho(1,wb_bins,:)) - mean(squeeze(rho(1,wb_bins,:)),2); % note: if obs would be a different data set than EOF, then remove mean that was removed for EOF
        alpha = (EOF_use' * EOF_use)\(EOF_use' * obs_use);
        if eof_option == 1
            wb_rho_recon = squeeze(EOFs_all(1,1:26,nb_eof_use))' * alpha + mean(rho(1,1:26,:),3)';
        elseif eof_option >= 2
            wb_rho_recon = squeeze(EOFs_all(1,1:26,nb_eof_use)) * alpha + mean(rho(1,1:26,:),3)';
        end
        wb_rho_grid = interp1(moor_dep(1,1:26),wb_rho_recon,dep_grid,'spline');
    end

    % EB
    if ebr_option ~= 0 && ebr_option ~= 4
        EOF_use = squeeze(EOFs_all(2,eb_bins, nb_eof_use));
        if eof_option == 1
            EOF_use = EOF_use';
        end
        obs_use = squeeze(rho(2,eb_bins,:)) - mean(squeeze(rho(2,eb_bins,:)),2);
        alpha = (EOF_use' * EOF_use)\(EOF_use' * obs_use);
        if eof_option == 1
            eb_rho_recon = squeeze(EOFs_all(2,1:26,nb_eof_use))' * alpha + mean(rho(2,1:26,:),3)';
        elseif eof_option >= 2
            eb_rho_recon = squeeze(EOFs_all(2,1:26,nb_eof_use)) * alpha + mean(rho(2,1:26,:),3)';
        end
        eb_rho_grid = interp1(moor_dep(2,1:26),eb_rho_recon,dep_grid,'spline');
    end

elseif vst_option == 1
% interpolating between instrument levels by using climatological gradients of temperature and salinity (following Johns et al, 2005)
    
    % calculation of climatological gradients
    % 1. T,S profils for each month of the year
    [~,time_mly_mm,~] = datevec(time_mly);
    ptemp_ann = ones(2,size(moor_dep,2),length(time_mly)).*NaN;
    salt_ann = ones(2,size(moor_dep,2),length(time_mly)).*NaN;
    for i = 1:length(time_mly)
        ptemp_ann(:,:,i) = mean(ptemp(:,:,time_mly_mm==time_mly_mm(i)),3);
        salt_ann(:,:,i) = mean(salt_orig(:,:,time_mly_mm==time_mly_mm(i)),3);
    end

    % 2. vertical gradients
    pres_ddp = ones(2,size(moor_dep,2)-1).*NaN;
    ptemp_ddp = ones(2,size(moor_dep,2)-1,length(time_mly)).*NaN;
    for iz = 1:size(moor_dep,2)-1
        pres_ddp(:,iz) = mean(moor_pres(:,iz:iz+1),2);
        ptemp_ddp(:,iz,:) = mean(ptemp_ann(:,iz:iz+1,:),2); % would be needed if gridded on temp instead of pres
    end

    dTdp_wb = diff(squeeze(ptemp_ann(1,:,:)))./diff(moor_pres(1,:))';
    dTdp_eb = diff(squeeze(ptemp_ann(2,:,:)))./diff(moor_pres(2,:))';
    dSdp_wb = diff(squeeze(salt_ann(1,:,:)))./diff(moor_pres(1,:))';
    dSdp_eb = diff(squeeze(salt_ann(2,:,:)))./diff(moor_pres(2,:))';

    save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/moor_dTdp_dSdp.mat'], 'time_mly','ptemp_ddp','pres_ddp','dTdp_wb','dTdp_eb','dSdp_wb','dSdp_eb','-v7.3')

    % 3. vertical gradients as a function of pressure (not temperature as in McCarthy 2015)
    pres_gridhr = 0:1870;
    dTdp_gridhr_wb = ones(length(pres_gridhr),length(time_mly)).*NaN;
    dTdp_gridhr_eb = ones(length(pres_gridhr),length(time_mly)).*NaN;
    dSdp_gridhr_wb = ones(length(pres_gridhr),length(time_mly)).*NaN;
    dSdp_gridhr_eb = ones(length(pres_gridhr),length(time_mly)).*NaN;
    for i = 1:length(time_mly)
        dTdp_gridhr_wb(:,i) = interp1(squeeze(pres_ddp(1,:)),dTdp_wb(:,i),pres_gridhr,'linear',0);
        dTdp_gridhr_eb(:,i) = interp1(squeeze(pres_ddp(2,1:25)),dTdp_eb(1:25,i),pres_gridhr,'linear','extrap');
        dTdp_gridhr_eb(1:10,i) = interp1(squeeze(pres_ddp(2,1:25)),dTdp_eb(1:25,i),pres_gridhr(1:10),'linear',0);
        dSdp_gridhr_wb(:,i) = interp1(squeeze(pres_ddp(1,:)),dSdp_wb(:,i),pres_gridhr,'linear',0);
        dSdp_gridhr_eb(:,i) = interp1(squeeze(pres_ddp(2,1:25)),dSdp_eb(1:25,i),pres_gridhr,'linear','extrap');
        dSdp_gridhr_eb(1:10,i) = interp1(squeeze(pres_ddp(2,1:25)),dSdp_eb(1:25,i),pres_gridhr(1:10),'linear',0);
    end

    % calculation of T&S on desired grid
    if wbr_option > 0
        bb_1 = 1;
    else
        bb_1 = 2;
    end
    if ebr_option ~= 0 && ebr_option ~= 4
        bb_2 = 2;
    else
        bb_2 = 1;
    end

    for bb = bb_1:bb_2 % for both boundaries
        if bb == 1
            bb_use = 'wb';
        elseif bb == 2
            bb_use = 'eb';
        end
        eval(['bins_use = ' bb_use '_bins;'])
        eval(['dTdp_gridhr_use = dTdp_gridhr_' bb_use ';'])
        eval(['dSdp_gridhr_use = dSdp_gridhr_' bb_use ';'])

        ptemp_grid_use = ones(length(dep_grid),length(time_mly)).*NaN;
        salt_grid_use = ones(length(dep_grid),length(time_mly)).*NaN;

        for iz = 1:length(dep_grid)
            % finding adjacent instruments
            instr_upper = bins_use(find(moor_dep(bb,bins_use)<=dep_grid(iz),1,'last'));
            instr_lower = bins_use(find(moor_dep(bb,bins_use)>dep_grid(iz),1,'first'));

            if ~isempty(instr_lower)
                % actual interpolation / stepwise integration
                weight_1 = 1 - (abs(pres_grid(iz)-moor_pres(bb,instr_upper))) / (moor_pres(bb,instr_lower)-moor_pres(bb,instr_upper));
                weight_2 = 1 - (abs(pres_grid(iz)-moor_pres(bb,instr_lower))) / (moor_pres(bb,instr_lower)-moor_pres(bb,instr_upper));
                
                pp_1 = find(pres_gridhr <= round(pres_grid(iz)) & pres_gridhr >= round(moor_pres(bb,instr_upper)));
                ptemp_step_1 = trapz(pres_gridhr(pp_1),dTdp_gridhr_use(pp_1,:),1);
                salt_step_1 = trapz(pres_gridhr(pp_1),dSdp_gridhr_use(pp_1,:),1);
                
                pp_2 = find(pres_gridhr >= round(pres_grid(iz)) & pres_gridhr <= round(moor_pres(bb,instr_lower)));
                ptemp_step_2 = -trapz(pres_gridhr(pp_2),dTdp_gridhr_use(pp_2,:),1);
                salt_step_2 = -trapz(pres_gridhr(pp_2),dSdp_gridhr_use(pp_2,:),1);
    
                ptemp_grid_use(iz,:) = weight_1.* (squeeze(ptemp(bb,instr_upper,:))' + ptemp_step_1) + weight_2.* (squeeze(ptemp(bb,instr_lower,:))' + ptemp_step_2);
                salt_grid_use(iz,:) = weight_1.* (squeeze(salt(bb,instr_upper,:))' + salt_step_1) + weight_2.* (squeeze(salt(bb,instr_lower,:))' + salt_step_2);
            else
                pp_1 = find(pres_gridhr <= round(pres_grid(iz)) & pres_gridhr >= round(moor_pres(bb,instr_upper)));
                ptemp_grid_use(iz,:) = squeeze(ptemp(bb,instr_upper,:))' + trapz(pres_gridhr(pp_1),dTdp_gridhr_use(pp_1,:),1);
                salt_grid_use(iz,:) = squeeze(salt(bb,instr_upper,:))' + trapz(pres_gridhr(pp_1),dSdp_gridhr_use(pp_1,:),1);
            end
        end
        
        % include surface option
        if sur_option == 3
            % 1. determine reference depth as depth of first instrument
            i_zr = round(moor_dep(bb,bins_use(2)),-1)/z_spacing+1;
            % 2. linear extrapolation of S
            salt_lin_extrap = interp1(moor_dep(bb,bins_use(2:3))',squeeze(salt(bb,bins_use(2:3),:)),dep_grid(1:i_zr)','linear','extrap');
            % 3. remove linear extrapolation from climatological salt year
            salt_lin_extrap_annomaly = interp1(moor_dep(bb,:)',squeeze(salt_ann(bb,:,:)),dep_grid(1:i_zr)','spline') - salt_lin_extrap;
            % 4. fit of quadratic and cubic term
            model = @(a, z) a(1) * (z - dep_grid(i_zr)).^2 + a(2) * (z - dep_grid(i_zr)).^3;
            a0 = [1, 1]; % initial guess
            salt_residual_fit = ones(i_zr,length(time_mly)) .* NaN;
            for i = 1:length(time_mly)
                a_estimate = lsqcurvefit(model, a0, dep_grid(1:i_zr)', salt_lin_extrap_annomaly(:,i)); % % Fit using lsqcurvefit
                salt_residual_fit(:,i) = model(a_estimate, dep_grid(1:i_zr)');
            end
            % 5. add linear and residual fit
            salt_final_extrap = salt_lin_extrap + salt_residual_fit;
            salt_grid_use(1:i_zr,:) = salt_final_extrap;
        end

        % calculating density for interpolated profiles
        temp = sw_temp(salt_grid_use,ptemp_grid_use,pres_grid',0);
        rho_grid_use = sw_dens(salt_grid_use,temp,pres_grid');
        eval([bb_use '_rho_grid = rho_grid_use;'])
    end


elseif vst_option == 2
% interpolating between instrument levels by using anomalies relative
% to climatological profiles (following Williams et al, 2015)
    
    % calculation of reference profiles for the anomalies
    ptemp_ref = mean(ptemp,3);
    salt_ref = mean(salt_orig,3);
    wb_ptemp_ref_F = griddedInterpolant(moor_pres(1,:),ptemp_ref(1,:),'spline');
    eb_ptemp_ref_F = griddedInterpolant(moor_pres(2,1:26),ptemp_ref(2,1:26),'spline');
    wb_salt_ref_F = griddedInterpolant(moor_pres(1,:),salt_ref(1,:),'spline');
    eb_salt_ref_F = griddedInterpolant(moor_pres(2,1:26),salt_ref(2,1:26),'spline');
    
    if wbr_option > 0 % WB
        % calculation of anomalies for subsampled data
        wb_ptemp_anm = squeeze(ptemp(1,wb_bins,:) - wb_ptemp_ref_F(moor_pres(1,wb_bins)));
        wb_salt_anm = squeeze(salt(1,wb_bins,:) - wb_salt_ref_F(moor_pres(1,wb_bins)));
        if sur_option == 2
            wb_salt_anm(1,:) = wb_salt_anm(2,:);
        end
        % interpolating on regular pressure levels
        % --linear
        wb_ptemp_grid = interp1(moor_pres(1,wb_bins),wb_ptemp_anm,pres_grid,'linear','extrap') + wb_ptemp_ref_F(pres_grid)';
        wb_salt_grid = interp1(moor_pres(1,wb_bins),wb_salt_anm,pres_grid,'linear','extrap') + wb_salt_ref_F(pres_grid)';
        % calculating density for interpolated profiles
        temp = sw_temp(wb_salt_grid,wb_ptemp_grid,pres_grid',0);
        wb_rho_grid = sw_dens(wb_salt_grid,temp,pres_grid');
    end

    if ebr_option ~= 0 && ebr_option ~= 4 % same for EB
        eb_ptemp_anm = squeeze(ptemp(2,eb_bins,:) - eb_ptemp_ref_F(moor_pres(2,eb_bins)));
        eb_salt_anm = squeeze(salt(2,eb_bins,:) - eb_salt_ref_F(moor_pres(2,eb_bins)));
        if sur_option == 2
            eb_salt_anm(1,:) = eb_salt_anm(2,:);
        end
        eb_ptemp_grid = interp1(moor_pres(2,eb_bins),eb_ptemp_anm,pres_grid,'linear','extrap') + eb_ptemp_ref_F(pres_grid)';
        eb_salt_grid = interp1(moor_pres(2,eb_bins),eb_salt_anm,pres_grid,'linear','extrap') + eb_salt_ref_F(pres_grid)';
        temp = sw_temp(eb_salt_grid,eb_ptemp_grid,pres_grid',0);
        eb_rho_grid = sw_dens(eb_salt_grid,temp,pres_grid');
    end

end

% --- computing Tg(z,t) and AMOCg(t) ---
g = sw_g(lat_mean,0); % [m/s2] Earth gravitational acceleration
f = sw_f(lat_mean);
rho0 = 1025; % [kg/m3] reference density (value from NEMO, e.g. cdfgeostrophy)
irefdep = 1100/10 +1; % depth until which transport is calculated

Tg = ones(irefdep,length(time_mly)).*NaN;
if idr_option == 1 % integrate density from surface
    for iz = 1:irefdep
        Tg(iz,:) = g/(rho0*f) .* (trapz(dep_grid(1:iz),eb_rho_grid(1:iz,:)) - trapz(dep_grid(1:iz),wb_rho_grid(1:iz,:)) + rho0.*(ssh(2,:)-ssh(1,:)));
    end
elseif idr_option == 2 % integrate density from level of nomo
    for iz = 1:irefdep-1
        Tg(iz,:) = - g/(rho0*f) .* trapz(dep_grid(iz:irefdep),eb_rho_grid(iz:irefdep,:) - wb_rho_grid(iz:irefdep,:));
    end
    Tg(end,:) = 0;
elseif idr_option == 5 % integrate WB4 from surface and EB1 from BP1200
    bp_eb_1200 = g .*  trapz(dep_grid,eb_rho_grid_orig) + g.*rho0 .* ssh(2,:);

    for iz = 1:irefdep
        Tg(iz,:) = g/(rho0*f) .* (bp_eb_1200./g - trapz(dep_grid(iz:end),eb_rho_grid(iz:end,:)) - trapz(dep_grid(1:iz),wb_rho_grid(1:iz,:)) - rho0.*ssh(1,:));
    end
elseif idr_option == 3 % same as 5 but with bp highpass filtered and then anomalies
    bprs = load(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/bpr_pos22_prc1_ebb0_ref1_vst1_stp4_' nmid '.mat']); % load de-drifted BPR values, caution: in dbar

    % highpass filter (for density contribution to match bp)
    Fs = 1/(30*24*3600); % Sampling frequency in Hz
    Fc = 1/(2*365*24*3600); % Cutoff frequency in Hz
    [b, a] = butter(3, Fc/(Fs/2), 'high'); % 3rd-order low-pass filter

    for iz = 1:irefdep
        Tg(iz,:) = g/(rho0*f) .* (bprs.bpa_eb1200_mly .* 1e4 ./g - filtfilt(b,a,trapz(dep_grid(iz:end),eb_rho_grid(iz:end,:))) - trapz(dep_grid(1:iz),wb_rho_grid(1:iz,:)) - rho0.*ssh(1,:) + mean(filtfilt(b,a,trapz(dep_grid(iz:end),eb_rho_grid(iz:end,:))) + trapz(dep_grid(1:iz),wb_rho_grid(1:iz,:)) + rho0.*ssh(1,:)));
    end
elseif idr_option == 4 % integrate both WB4 and EB1 relative to BPRs
    bprs = load(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/bpr_pos22_prc1_ebb0_ref1_vst1_stp4_' nmid '.mat']); % load de-drifted BPR values, caution: in dbar
    iz_bp_wb500 = 440/10 +1;
    
    % highpass filter (for density contribution to match bp)
    Fs = 1/(30*24*3600); % Sampling frequency in Hz
    Fc = 1/(2*365*24*3600); % Cutoff frequency in Hz
    [b, a] = butter(3, Fc/(Fs/2), 'high'); % 3rd-order low-pass filter
    
    bp_wb500_vdx = f .* rho0 .* squeeze(sum(v_WB_mly,1,'omitnan')) .* x_spacing; 

    for iz = 1:irefdep
        pEB1_part = bprs.bpa_eb1200_mly .* 1e4 - g.* filtfilt(b,a,trapz(dep_grid(iz:end),eb_rho_grid(iz:end,:))) + mean(g.* filtfilt(b,a,trapz(dep_grid(iz:end),eb_rho_grid(iz:end,:))));
        if iz < iz_bp_wb500
            pWB4_part = bprs.bpa_wb500_mly .* 1e4 + filtfilt(b,a,bp_wb500_vdx) - g.* filtfilt(b,a,trapz(dep_grid(iz:iz_bp_wb500),eb_rho_grid(iz:iz_bp_wb500,:))) - mean(filtfilt(b,a,bp_wb500_vdx)) + mean(g.* filtfilt(b,a,trapz(dep_grid(iz:iz_bp_wb500),eb_rho_grid(iz:iz_bp_wb500,:))));
        elseif iz == iz_bp_wb500
            pWB4_part = bprs.bpa_wb500_mly .* 1e4 + filtfilt(b,a,bp_wb500_vdx) - mean(filtfilt(b,a,bp_wb500_vdx));
        elseif iz > iz_bp_wb500
            pWB4_part = bprs.bpa_wb500_mly .* 1e4 + filtfilt(b,a,bp_wb500_vdx) + g.* filtfilt(b,a,trapz(dep_grid(iz_bp_wb500:iz),wb_rho_grid(iz_bp_wb500:iz,:))) - mean(filtfilt(b,a,bp_wb500_vdx)) - mean(g.* filtfilt(b,a,trapz(dep_grid(iz_bp_wb500:iz),wb_rho_grid(iz_bp_wb500:iz,:))));
        end
        Tg(iz,:) = 1/(rho0*f) .* (pEB1_part - pWB4_part);
    end
end
clear iz

AMOCg = sum(Tg.*z_spacing,1).*1e-6;


% --- save variables in interim folder, then download to plot (then script plot_MooringSetup) ---
% transport
save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/moor_wbr' num2str(wbr_option) 'ebr' num2str(ebr_option) '_vst' num2str(vst_option) 'eof' num2str(eof_option) 'sur' num2str(sur_option) '_idr' num2str(idr_option) '_' nmid test_suffix '.mat'],'time_mly','dep_grid','Tg','AMOCg','-v7.3')
% density
save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/moor_density_wbr' num2str(wbr_option) 'ebr' num2str(ebr_option) '_vst' num2str(vst_option) 'eof' num2str(eof_option) 'sur' num2str(sur_option) '_idr' num2str(idr_option) '_' nmid test_suffix '.mat'],'time_mly','ssh','dep_grid','wb_rho_grid','eb_rho_grid','wb_rho_grid_orig','eb_rho_grid_orig','-v7.3')
