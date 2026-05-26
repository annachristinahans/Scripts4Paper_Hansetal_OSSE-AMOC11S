function f_cut_BPRSetup(mod_option,posbpr_option,possla_option,prc_option,ebb_option,ref_option,vst_option,stp_option,test_suffix)

% function to cut INALT/VIKING BPR AMOCg transport and to be run on nesh

% input:
% mod_option: choose model INALT20 (1) or VIKING20X (2)
% posbpr_option: choose which postition to use for BPR at WB (1 or 2)
% possla_option: choose SLA at BPR 300m/1200m position (1) or at the boundary (2)
% prc_option: choose processing of BPRs as
%       original timeseries (0), every 2.3yrs detrended (1), 2 yrs highpass filtered (2),
%       original values + WB low drift (AZA version) corresponding to Sonardyne maximum (31), RBR maximum (32), RBR median (33), fictive (34),
%       every 2.3yrs detrended + WB IES annual data (4),
%       EB 2.3yrs detrended and WB with velocity field and WB4 current setup density (5),
%       EB 2.3yrs detrended and WB with velocity field and WB4 enhanced setup density (6),
%       or like 2 and also SSH highpass filtered (7)
% ebb_option: choose EB gap filling as
%       using original data (0), annual + semiannual harmonics (1), reconstruction from 1200m with density (2),
%       same as 2 with with density contribution high pass filtered (3)
% ref_option: choose reference level at 1100m (1), 900m (2), or 1300m (3)
% vst_option: choose vertical structure as
%       EOF regression with no motion (1), EOF regression with known motion (2),
%       linear interpolation (3), spline interpolation (4),
%       EOF regression and linear interpolation combined - possible for stp2 or stp6 (5)
% stp_option: choose setup as reconstruction of
%       300m and 500m (1), only 500m (2), only 300m (3), 300m, 500m and 1200m (4), 500m and 1200m (5), only 1200m (6)
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

% --- preparing model data
if mod_option == 1
        files1 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*T_11S_box.ncbotpressure_detrssh.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*T_11S_box.nc');
        files3 = dir('/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/*V_11S_box.nc');
        mmask = '/gxfs_work/geomar/smomw628/INALT20.L46-KFS119_11S/1_INALT20.L46-KFS119_mesh_mask_11S_box.nc';
elseif mod_option == 2
        files1 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*T_11S_box.ncbotpressure_detrssh.nc');
        files2 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*T_11S_box.nc');
        files3 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*V_11S_box.nc');
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

% choose positions of the BPRs
% EB
[~,ilat_eb300] = min(abs(lat+10.68));
[~,ilon_eb300] = min(abs(lon-13.23));
waterdep_eb300 = sum(squeeze(dep_cell_size(ilon_eb300,ilat_eb300,1:16)));
[~,ilat_eb500] = min(abs(lat+10.71));
[~,ilon_eb500] = min(abs(lon-13.19));
waterdep_eb500 = sum(squeeze(dep_cell_size(ilon_eb500,ilat_eb500,1:19)));
[~,ilat_eb1200] = min(abs(lat+10.83));
[~,ilon_eb1200] = min(abs(lon-13.00));
ilon_eb1200 = ilon_eb1200 - 1; % so that the mooring reaches a depth of 1200m
dep_eb1200 = squeeze(dep(ilon_eb1200,ilat_eb1200,1:26));
waterdep_eb1200 = sum(squeeze(dep_cell_size(ilon_eb1200,ilat_eb1200,1:26)));

% WB
[~,ilat_wb300] = min(abs(lat+10.23));
[~,ilon_wb300] = min(abs(lon+35.87));
[~,ilat_wb500] = min(abs(lat+10.23));
[~,ilon_wb500] = min(abs(lon+35.86));
% modification of WB positions so that depth are fitting better to BPRs
if posbpr_option == 1
    ilon_wb300 = ilon_wb300 - 1;
    ilat_wb500 = ilat_wb500 + 1;
elseif posbpr_option == 2
    ilon_wb300 = ilon_wb300 - 2;
    ilat_wb300 = ilat_wb300 - 2;
    ilon_wb500 = ilon_wb500 - 1;
    ilat_wb500 = ilat_wb500 - 1;
end
dep_wb300 = squeeze(dep(ilon_wb300,ilat_wb300,1:16));
waterdep_wb300 = sum(squeeze(dep_cell_size(ilon_wb300,ilat_wb300,1:16)));
dep_wb500 = squeeze(dep(ilon_wb500,ilat_wb500,1:19));
waterdep_wb500 = sum(squeeze(dep_cell_size(ilon_wb500,ilat_wb500,1:19)));
[~,ilat_wb1200] = min(abs(lat+10.27)); % equals in model position of K1, (use ilat-1 for next depth step at 1450m)
[~,ilon_wb1200] = min(abs(lon+35.86));
dep_wb1200 = squeeze(dep(ilon_wb1200,ilat_wb1200,1:25));
waterdep_wb1200 = sum(squeeze(dep_cell_size(ilon_wb1200,ilat_wb1200,1:25)));

% WB4 mooring
latV = double(ncread(mmask,'gphiv'));
latV = latV(1,:);
[~,ilat_WB4] = min(abs(latV+10.94));
[~,ilon_WB4] = min(abs(lon+34.99));
depV = double(ncread(mmask,'gdepv')); % in m, including partial cells
dep_WB4 = squeeze(depV(ilon_WB4,ilat_WB4,1:29));
x_spacing = double(ncread(mmask,'e1v')); % in m
x_spacing = x_spacing(1,ilat_WB4);
qcv = double(ncread(mmask,'vmask')); % 0 land, 1 water
qcv(qcv==0) = NaN;

% --- loading of data
countt = 1;
for i = 1:num_yy % loop over years

    tim = double(ncread([files1(i).folder '/' files1(i).name],'time_counter'));
    time(countt:countt-1+length(tim)) = datenum(1900,1,1,0,0,tim);
    
    bp = double(ncread([files1(i).folder '/' files1(i).name],'sobotpres')) .* 1e-4; % in dbar (10-4 Pa)
    ssh = double(ncread([files1(i).folder '/' files1(i).name],'sossheig')); % in m
    ptemp = double(ncread([files2(i).folder '/' files2(i).name],'votemper')); % in °C, potential temperature
    salt = double(ncread([files2(i).folder '/' files2(i).name],'vosaline')); % in 0.001, practical salinity
    
    % applying land-sea mask
    bp = bp .* qc(:,:,1);
    ssh = ssh .* qc(:,:,1);
    ptemp = ptemp .* qc;
    salt = salt .* qc;


    % selecting bottom pressure
    bp_wb300(countt:countt-1+length(tim)) = squeeze(bp(ilon_wb300,ilat_wb300,:));
    bp_wb500(countt:countt-1+length(tim)) = squeeze(bp(ilon_wb500,ilat_wb500,:));
    bp_wb1200(countt:countt-1+length(tim)) = squeeze(bp(ilon_wb1200,ilat_wb1200,:));
    
    bp_eb300(countt:countt-1+length(tim)) = squeeze(bp(ilon_eb300,ilat_eb300,:));
    bp_eb500(countt:countt-1+length(tim)) = squeeze(bp(ilon_eb500,ilat_eb500,:));
    bp_eb1200(countt:countt-1+length(tim)) = squeeze(bp(ilon_eb1200,ilat_eb1200,:));

    % selecting T,S,SSH for IES and/or moored density reconstruction
    ptemp_eb1200(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_eb1200,ilat_eb1200,1:26,:));
    salt_eb1200(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_eb1200,ilat_eb1200,1:26,:));
    ssh_eb1200(countt:countt-1+length(tim)) = squeeze(ssh(ilon_eb1200,ilat_eb1200,:));

    ptemp_wb300(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_wb300,ilat_wb300,1:16,:));
    salt_wb300(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_wb300,ilat_wb300,1:16,:));
    ssh_wb300(countt:countt-1+length(tim)) = squeeze(ssh(ilon_wb300,ilat_wb300,:));
    ptemp_wb500(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_wb500,ilat_wb500,1:19,:));
    salt_wb500(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_wb500,ilat_wb500,1:19,:));
    ssh_wb500(countt:countt-1+length(tim)) = squeeze(ssh(ilon_wb500,ilat_wb500,:));
    ptemp_wb1200(:,countt:countt-1+length(tim)) = squeeze(ptemp(ilon_wb1200,ilat_wb1200,1:25,:));
    salt_wb1200(:,countt:countt-1+length(tim)) = squeeze(salt(ilon_wb1200,ilat_wb1200,1:25,:));
    ssh_wb1200(countt:countt-1+length(tim)) = squeeze(ssh(ilon_wb1200,ilat_wb1200,:));


    % choosing SSH position
    if possla_option == 1
        if ebb_option == 0 || ebb_option == 1
            ilat_ssh_eb = ilat_eb300;
            ilon_ssh_eb = ilon_eb300;
        else
            ilat_ssh_eb = ilat_eb1200;
            iont_ssh_eb = ilon_eb1200;
        end
        ilon_ssh_wb = ilon_wb300;
        ilat_ssh_wb = ilat_wb300;
    elseif possla_option == 2
        if ebb_option == 0 || ebb_option == 1
            ilat_ssh_eb = ilat_eb300;
        else
            ilat_ssh_eb = ilat_eb1200;
        end
        ilon_ssh_eb = find(~isnan(squeeze(ssh(:,ilat_ssh_eb,1))),1,'last');
        ilat_ssh_wb = ilat_wb300;
        ilon_ssh_wb = find(~isnan(squeeze(ssh(:,ilat_ssh_wb,1))),1,'first');
    end

    ssh_eb(countt:countt-1+length(tim)) = squeeze(ssh(ilon_ssh_eb,ilat_ssh_eb,:));
    ssh_wb(countt:countt-1+length(tim)) = squeeze(ssh(ilon_ssh_wb,ilat_ssh_wb,:));

    if prc_option == 5 || prc_option == 6
        ssh_WB4(countt:countt-1+length(tim)) = squeeze(ssh(ilon_WB4,ilat_WB4,:));

        v_use = double(ncread([files3(i).folder '/' files3(i).name],'vomecrty')); % in m/s
        v_use = v_use .* qcv;
        v_WB(:,:,countt:countt-1+length(tim)) = squeeze(v_use(1:ilon_WB4,ilat_WB4,1:29,:));
    end

    countt = countt + length(tim);
end
clear i files1 files2 mmask tim qc bp ptemp salt ssh countt

% --- processing and reconstrucing BPR data

f = sw_f(-10.5); % [rad/s] Coriolis parameter
g = sw_g(-10.5,0); % [m/s2] Earth gravitational acceleration
rho0 = 1025; % [kg/m3] reference density (value from NEMO, e.g. cdfgeostrophy)


% BPR processing
if prc_option == 0
    bp_eb300_dtr = bp_eb300;
    bp_eb500_dtr = bp_eb500;
    bp_eb1200_dtr = bp_eb1200;
    bp_wb300_dtr = bp_wb300;
    bp_wb500_dtr = bp_wb500;
    bp_wb1200_dtr = bp_wb1200;

elseif prc_option > 30
    % adding low drift (AZA version)
    if prc_option == 31 % Sonardyne: better than 2mm/year
        low_drift = 0.002/365.25 * (time-time(1)); % a
    elseif prc_option == 32 % RBR: better than 1mm/year
        low_drift = 0.001/365.25 * (time-time(1)); % b
    elseif prc_option == 33 % RBR: median of 7e-5m/year
        low_drift = 0.00007/365.25 * (time-time(1)); % c
    elseif prc_option == 34 % fictive - not realistic
        low_drift = 0.00001/365.25 * (time-time(1)); % d
    end
    % apply drift for each boundary same direction so that it does not cancel out
    bp_wb300_dtr = bp_wb300 + low_drift;
    bp_wb500_dtr = bp_wb500 + low_drift;
    bp_wb1200_dtr = bp_wb1200 + low_drift;
    if ebb_option == 1 % assuming no AZAs at EB
        bp_eb300_dtr = bp_eb300;
        bp_eb500_dtr = bp_eb500;
        bp_eb1200_dtr = bp_eb1200;
    else % assuming AZAs at EB
        bp_eb300_dtr = bp_eb300 - low_drift;
        bp_eb500_dtr = bp_eb500 - low_drift;
        bp_eb1200_dtr = bp_eb1200 - low_drift;
    end

elseif prc_option == 1 || prc_option == 4 || prc_option == 5 || prc_option == 6
    % detrending BPRs every 2.2/2.3 yrs (simulate deployment periods)
    if mod_option == 1
        depl_length = 811; % in days (2.2 years)
    elseif mod_option == 2
        depl_length = 845; % in days (2.3 years)
    end
    time_depl = time(1):depl_length:time(end);
    for i = 1:length(time_depl)-1
        ind = find(time>=time_depl(i) & time<time_depl(i+1));
        % compute annual harmonics
        yout = anharm_bp(time(ind), time(ind), bp_eb300(ind),2);
        % compute exp-lin trend
        [~,bp_exfit]  = exp_lin_fit_bp(time(ind),bp_eb300(ind) - yout,[0.5 0.05 0.005 0]);
        % remove exp-lin trend and mean
        bp_eb300_dtr(ind) = bp_eb300(ind) - bp_exfit - mean(bp_eb300(ind) - bp_exfit);
        
        % same thing for other bprs
        yout = anharm_bp(time(ind), time(ind), bp_eb500(ind),2);
        [~,bp_exfit]  = exp_lin_fit_bp(time(ind),bp_eb500(ind) - yout,[0.5 0.05 0.005 0]);
        bp_eb500_dtr(ind) = bp_eb500(ind) - bp_exfit - mean(bp_eb500(ind) - bp_exfit);

        yout = anharm_bp(time(ind), time(ind), bp_eb1200(ind),2);
        [~,bp_exfit]  = exp_lin_fit_bp(time(ind),bp_eb1200(ind) - yout,[0.5 0.05 0.005 0]);
        bp_eb1200_dtr(ind) = bp_eb1200(ind) - bp_exfit - mean(bp_eb1200(ind) - bp_exfit);
        
        yout = anharm_bp(time(ind), time(ind), bp_wb300(ind),2);
        [~,bp_exfit]  = exp_lin_fit_bp(time(ind),bp_wb300(ind) - yout,[0.5 0.05 0.005 0]);
        bp_wb300_dtr(ind) = bp_wb300(ind) - bp_exfit - mean(bp_wb300(ind) - bp_exfit);
        
        yout = anharm_bp(time(ind), time(ind), bp_wb500(ind),2);
        [~,bp_exfit]  = exp_lin_fit_bp(time(ind),bp_wb500(ind) - yout,[0.5 0.05 0.005 0]);
        bp_wb500_dtr(ind) = bp_wb500(ind) - bp_exfit - mean(bp_wb500(ind) - bp_exfit);
        
        yout = anharm_bp(time(ind), time(ind), bp_wb1200(ind),2);
        [~,bp_exfit]  = exp_lin_fit_bp(time(ind),bp_wb1200(ind) - yout,[0.5 0.05 0.005 0]);
        bp_wb1200_dtr(ind) = bp_wb1200(ind) - bp_exfit - mean(bp_wb1200(ind) - bp_exfit);
    end
    clear ind yout bp_exfit
    % adapt time vector to deployment periods also for other variables
    time = time(1:find(time==time_depl(i+1))-1);
    ssh_eb = ssh_eb(1:length(time));
    ssh_wb  = ssh_wb(1:length(time));
    ptemp_wb300 = ptemp_wb300(:,1:length(time));
    salt_wb300 = salt_wb300(:,1:length(time));
    ssh_wb300  = ssh_wb300(1:length(time));
    ptemp_wb500 = ptemp_wb500(:,1:length(time));
    salt_wb500 = salt_wb500(:,1:length(time));
    ssh_wb500  = ssh_wb500(1:length(time));
    ptemp_wb1200 = ptemp_wb1200(:,1:length(time));
    salt_wb1200 = salt_wb1200(:,1:length(time));
    ssh_wb1200  = ssh_wb1200(1:length(time));
    ptemp_eb1200 = ptemp_eb1200(:,1:length(time));
    salt_eb1200 = salt_eb1200(:,1:length(time));
    ssh_eb1200  = ssh_eb1200(1:length(time));

    if prc_option == 4 % additionally add IES reconstruction for the WB
        % preparing regular 10m grid for travel time calculation
        delta_z_ies = 10;
        dep_wb300_grid_ies = 0:delta_z_ies:270;
        dep_wb500_grid_ies = 0:delta_z_ies:440;
        dep_wb1200_grid_ies = 0:delta_z_ies:1100;
        
        % calculating reconstructed yearly bottom pressure
        [bp_wb300_reconstr_yly,bp_wb300_yly,yyears] = bp_from_ies_yly(dep_wb300,ptemp_wb300,salt_wb300,bp_wb300,dep_wb300_grid_ies,ilat_wb300, lat,time,delta_z_ies);
        [bp_wb500_reconstr_yly,bp_wb500_yly,~] = bp_from_ies_yly(dep_wb500,ptemp_wb500,salt_wb500,bp_wb500,dep_wb500_grid_ies,ilat_wb500, lat,time,delta_z_ies);
        [bp_wb1200_reconstr_yly,bp_wb1200_yly,~] = bp_from_ies_yly(dep_wb1200,ptemp_wb1200,salt_wb1200,bp_wb1200,dep_wb1200_grid_ies,ilat_wb1200, lat,time,delta_z_ies);
        
        % interpolating yearly to daily data while conserving the annual mean
        bp_wb300_reconstr_dly = bp_from_ies_yly2dly(bp_wb300_reconstr_yly,time,yyears);
        bp_wb500_reconstr_dly = bp_from_ies_yly2dly(bp_wb500_reconstr_yly,time,yyears);
        bp_wb1200_reconstr_dly = bp_from_ies_yly2dly(bp_wb1200_reconstr_yly,time,yyears);
        
        % adding reconstructed bottom pressure to de-drifted data
        bp_wb300_dtr = bp_wb300_dtr + bp_wb300_reconstr_dly;
        bp_wb500_dtr = bp_wb500_dtr + bp_wb500_reconstr_dly;
        bp_wb1200_dtr = bp_wb1200_dtr + bp_wb1200_reconstr_dly;

%         dep_eb1200_grid_ies = 0:delta_z_ies:1250;
%         [bp_eb1200_reconstr_yly,~,~] = bp_from_ies_yly(dep_eb1200,ptemp_eb1200,salt_eb1200,bp_eb1200,dep_eb1200_grid_ies,ilat_eb1200, lat,time,delta_z_ies);
%         bp_3b1200_reconstr_dly = bp_from_ies_yly2dly(bp_eb1200_reconstr_yly,time,yyears);
%         bp_eb1200_dtr = bp_eb1200_dtr + bp_eb1200_reconstr_dly;
    end
    clear delta_z_ies dep_wb300_grid_ies dep_wb500_grid_ies dep_wb1200_grid_ies dep_eb1200_grid_ies

elseif prc_option == 2 || prc_option == 7
    % 2 yrs high pass filter for BPRs
    Fs = 1/(24*3600); % Sampling frequency in Hz
    Fc = 1/(2*365*24*3600); % Cutoff frequency in Hz
    [b, a] = butter(3, Fc/(Fs/2), 'high'); % 3rd-order low-pass filter
    bp_eb300_dtr = filtfilt(b, a, bp_eb300);
    bp_eb500_dtr = filtfilt(b, a, bp_eb500);
    bp_eb1200_dtr = filtfilt(b, a, bp_eb1200);
    bp_wb300_dtr = filtfilt(b, a, bp_wb300);
    bp_wb500_dtr = filtfilt(b, a, bp_wb500);
    bp_wb1200_dtr = filtfilt(b, a, bp_wb1200);

    if prc_option == 7
        % addionally high pass filter ssh
        ssh_eb = filtfilt(b, a, ssh_eb);
        ssh_wb = filtfilt(b, a, ssh_wb);
    end
end

% computing anomalies as absolute pressure differences are not known
sla_eb = ssh_eb - mean(ssh_eb);
sla_wb = ssh_wb - mean(ssh_wb);
bpa_eb300 = bp_eb300_dtr - mean(bp_eb300_dtr);
bpa_eb500 = bp_eb500_dtr - mean(bp_eb500_dtr);
bpa_eb1200 = bp_eb1200_dtr - mean(bp_eb1200_dtr);
bpa_wb300 = bp_wb300_dtr - mean(bp_wb300_dtr);
bpa_wb500 = bp_wb500_dtr - mean(bp_wb500_dtr);
bpa_wb1200 = bp_wb1200_dtr - mean(bp_wb1200_dtr);

% approximating the measurement depth
depth_bp300 = mean([waterdep_wb300 waterdep_eb300]);
depth_bp500 = mean([waterdep_wb500 waterdep_eb500]);
depth_bp1200 = mean([waterdep_wb1200 waterdep_eb1200]);

% dealing with data gaps at the EB I
if ebb_option == 1 % harmonics for EB, like in observations, harmonics are computed from first 2 years
    yout_eb300 = anharm_bp(time(time<datenum(1982,1,1)),time,bpa_eb300(time<datenum(1982,1,1)),2);
    bpa_eb300 = yout_eb300;
    yout_eb500 = anharm_bp(time(time<datenum(1982,1,1)),time,bpa_eb500(time<datenum(1982,1,1)),2);
    bpa_eb500 = yout_eb500;
end

% averaging to monthly data (so that geostrophy is valid)
[yy,mm,~] = datevec(time);
time_mly = unique(datenum(yy,mm,15));
bpa_eb300_mly = ones(1,length(time_mly)).*NaN;
bpa_eb500_mly = ones(1,length(time_mly)).*NaN;
bpa_eb1200_mly = ones(1,length(time_mly)).*NaN;
bpa_wb300_mly = ones(1,length(time_mly)).*NaN;
bpa_wb500_mly = ones(1,length(time_mly)).*NaN;
bpa_wb1200_mly = ones(1,length(time_mly)).*NaN;
sla_eb_mly = ones(1,length(time_mly)).*NaN;
sla_wb_mly = ones(1,length(time_mly)).*NaN;
if prc_option == 5 || prc_option == 6
    ssh_WB4_mly = ones(1,length(time_mly)).*NaN;
    v_WB_mly = ones(size(v_WB,1),size(v_WB,2),length(time_mly)).*NaN;
end

for i = 1:length(time_mly)
    bpa_eb300_mly(1,i) = mean(bpa_eb300(datenum(yy,mm,15)==time_mly(i)));
    bpa_eb500_mly(1,i) = mean(bpa_eb500(datenum(yy,mm,15)==time_mly(i)));
    bpa_eb1200_mly(1,i) = mean(bpa_eb1200(datenum(yy,mm,15)==time_mly(i)));
    bpa_wb300_mly(1,i) = mean(bpa_wb300(datenum(yy,mm,15)==time_mly(i)));
    bpa_wb500_mly(1,i) = mean(bpa_wb500(datenum(yy,mm,15)==time_mly(i)));
    bpa_wb1200_mly(1,i) = mean(bpa_wb1200(datenum(yy,mm,15)==time_mly(i)));
    sla_eb_mly(1,i) = mean(sla_eb(datenum(yy,mm,15)==time_mly(i)));
    sla_wb_mly(1,i) = mean(sla_wb(datenum(yy,mm,15)==time_mly(i)));
    if prc_option == 5 || prc_option == 6
        ssh_WB4_mly(1,i) = mean(ssh_WB4(datenum(yy,mm,15)==time_mly(i)));
        v_WB_mly(:,:,i) = mean(v_WB(:,:,datenum(yy,mm,15)==time_mly(i)),3);
    end
end

% combination with mooring data (as the density reconstruction is available monthly)
if prc_option == 5 || prc_option == 6 % replace WB by WB4 reconstruction
    bpa_wb300_mly_orig = ones(1,length(time_mly)).*NaN;
    bpa_wb500_mly_orig = ones(1,length(time_mly)).*NaN;
    bpa_wb1200_mly_orig = ones(1,length(time_mly)).*NaN;
    for i = 1:length(time_mly)
        bpa_wb300_mly_orig(1,i) = mean(bp_wb300(datenum(yy,mm,15)==time_mly(i)));
        bpa_wb500_mly_orig(1,i) = mean(bp_wb500(datenum(yy,mm,15)==time_mly(i)));
        bpa_wb1200_mly_orig(1,i) = mean(bp_wb1200(datenum(yy,mm,15)==time_mly(i)));
    end
    bpa_wb300_mly_orig = bpa_wb300_mly_orig - mean(bpa_wb300_mly_orig);
    bpa_wb500_mly_orig = bpa_wb500_mly_orig - mean(bpa_wb500_mly_orig);
    bpa_wb1200_mly_orig = bpa_wb1200_mly_orig - mean(bpa_wb1200_mly_orig);


    bpa_wb300_mly = bp_from_WB4(waterdep_wb300,dep_WB4,ssh_WB4_mly,v_WB_mly,nmid,g,f,rho0,x_spacing,prc_option);
    bpa_wb500_mly = bp_from_WB4(waterdep_wb500,dep_WB4,ssh_WB4_mly,v_WB_mly,nmid,g,f,rho0,x_spacing,prc_option);
    bpa_wb1200_mly = bp_from_WB4(waterdep_wb1200,dep_WB4,ssh_WB4_mly,v_WB_mly,nmid,g,f,rho0,x_spacing,prc_option);
end
clear i

if ebb_option == 2 || ebb_option == 3 % dealing with data gaps at the EB II -> reconstructing the EB BPRs at 300m and 500m

    rho_EB1 = load(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/moor_density_wbr0ebr1_vst1eof0sur1_idr1_' nmid '.mat']); % for initial setup
    eb_rho_grid = rho_EB1.eb_rho_grid;

    bpa_eb300_mly_orig = ones(1,length(time_mly)).*NaN;
    bpa_eb500_mly_orig = ones(1,length(time_mly)).*NaN;
    for i = 1:length(time_mly)
        bpa_eb300_mly_orig(1,i) = mean(bp_eb300(datenum(yy,mm,15)==time_mly(i)));
        bpa_eb500_mly_orig(1,i) = mean(bp_eb500(datenum(yy,mm,15)==time_mly(i)));
    end
    bpa_eb300_mly_orig = bpa_eb300_mly_orig - mean(bpa_eb300_mly_orig);
    bpa_eb500_mly_orig = bpa_eb500_mly_orig - mean(bpa_eb500_mly_orig);

    % reconstruct pressure at 300m and 500m depth (currently in dbar, not Pa!)
    if ebb_option == 2
        bpa_eb300_mly = bpa_eb1200_mly - g .* trapz(rho_EB1.dep_grid(round((depth_bp300+40)/10)+1:end),eb_rho_grid(round((depth_bp300+40)/10)+1:end,:)) .* 1e-4 + mean(g .* trapz(rho_EB1.dep_grid(round((depth_bp300+40)/10)+1:end),eb_rho_grid(round((depth_bp300+40)/10)+1:end,:)) .* 1e-4); % bin + 40m to fit mc depth
        bpa_eb500_mly = bpa_eb1200_mly - g .* trapz(rho_EB1.dep_grid(round(depth_bp500/10)+1:end),eb_rho_grid(round(depth_bp500/10)+1:end,:)) .* 1e-4 + mean(g .* trapz(rho_EB1.dep_grid(round(depth_bp500/10)+1:end),eb_rho_grid(round(depth_bp500/10)+1:end,:)) .* 1e-4);
    elseif ebb_option == 3
        % with high pass filtering the density contribution
        bpa_eb300_denscontr = g .* trapz(rho_EB1.dep_grid(round((depth_bp300+40)/10)+1:end),eb_rho_grid(round((depth_bp300+40)/10)+1:end,:)) .* 1e-4;
        bpa_eb500_denscontr = g .* trapz(rho_EB1.dep_grid(round(depth_bp500/10)+1:end),eb_rho_grid(round(depth_bp500/10)+1:end,:)) .* 1e-4;
        % 2 yrs high pass filter for BPRs
        Fs = 1/(30*24*3600); % Sampling frequency in Hz (daily values)
        Fc = 1/(2*365*24*3600); % Cutoff frequency in Hz
        [b, a] = butter(3, Fc/(Fs/2), 'high'); % 3rd-order low-pass filter
        bpa_eb300_mly_orig = filtfilt(b,a,bpa_eb300_mly_orig);
        bpa_eb500_mly_orig = filtfilt(b,a,bpa_eb500_mly_orig);
        bpa_eb300_denscontr_hp = filtfilt(b, a, bpa_eb300_denscontr);
        bpa_eb500_denscontr_hp = filtfilt(b, a, bpa_eb500_denscontr);
        bpa_eb300_mly = bpa_eb1200_mly - bpa_eb300_denscontr_hp;
        bpa_eb500_mly = bpa_eb1200_mly - bpa_eb500_denscontr_hp;
    end
end
clear rho_EB1 yy mm


% computing anomalies again as monthly averaging added a small bias
sla_eb_mly = sla_eb_mly - mean(sla_eb_mly);
sla_wb_mly = sla_wb_mly - mean(sla_wb_mly);
bpa_eb300_mly = bpa_eb300_mly - mean(bpa_eb300_mly);
bpa_eb500_mly = bpa_eb500_mly - mean(bpa_eb500_mly);
bpa_eb1200_mly = bpa_eb1200_mly - mean(bpa_eb1200_mly);
bpa_wb300_mly = bpa_wb300_mly - mean(bpa_wb300_mly);
bpa_wb500_mly = bpa_wb500_mly - mean(bpa_wb500_mly);
bpa_wb1200_mly = bpa_wb1200_mly - mean(bpa_wb1200_mly);

% --- deriving quantities to save and plot

% calculating the geostrophic transport at the 'observation' depth
Tg_0m = (g/f) .* (sla_eb_mly - sla_wb_mly);
Tg_300m =  (bpa_eb300_mly - bpa_wb300_mly) .* 1e4 ./ (f*rho0);
Tg_500m =  (bpa_eb500_mly - bpa_wb500_mly) .* 1e4 ./ (f*rho0);
Tg_1200m = (bpa_eb1200_mly - bpa_wb1200_mly) .* 1e4 ./ (f*rho0);

% setting up the vertical structure
if ref_option == 1
    reference_level = 1100;
elseif ref_option == 2
    reference_level = 900;
elseif ref_option == 3
    reference_level = 1300;
end

if vst_option == 3 || vst_option == 4 % interpolate Tg between given points
    if vst_option == 3
        interp_method = 'linear';
    elseif vst_option == 4
        interp_method = 'spline';
    end

    delta_z = 10;
    dep_grid = 0:delta_z:reference_level;
    Tg_basin = ones(length(dep_grid),length(time_mly)).*NaN;
    if stp_option == 1 % BPRs at 300m, 500m
        for ii = 1:length(time_mly)
            Tg_basin(:,ii) = interp1([0 depth_bp300 depth_bp500 reference_level],[Tg_0m(ii) Tg_300m(ii) Tg_500m(ii) 0],dep_grid,interp_method,'extrap');
        end
    elseif stp_option == 2 % BPRs at 500m
        for ii = 1:length(time_mly)
            Tg_basin(:,ii) = interp1([0 depth_bp500 reference_level],[Tg_0m(ii) Tg_500m(ii) 0],dep_grid,interp_method,'extrap');
        end
    elseif stp_option == 3 % BPRs at 300m
        for ii = 1:length(time_mly)
            Tg_basin(:,ii) = interp1([0 depth_bp300 reference_level],[Tg_0m(ii) Tg_300m(ii) 0],dep_grid,interp_method,'extrap');
        end
    elseif stp_option == 4 % BPRs at 300m, 500m, 1200m
        for ii = 1:length(time_mly)
            Tg_basin(:,ii) = interp1([0 depth_bp300 depth_bp500 depth_bp1200],[Tg_0m(ii) Tg_300m(ii) Tg_500m(ii) Tg_1200m(ii)],dep_grid,interp_method,'extrap');
        end
    elseif stp_option == 5 % BPRs at 500m, 1200m
        for ii = 1:length(time_mly)
            Tg_basin(:,ii) = interp1([0 depth_bp500 depth_bp1200],[Tg_0m(ii) Tg_500m(ii) Tg_1200m(ii)],dep_grid,interp_method,'extrap');
        end
    elseif stp_option == 6 % BPRs at 1200m
        for ii = 1:length(time_mly)
            Tg_basin(:,ii) = interp1([0 depth_bp1200],[Tg_0m(ii) Tg_1200m(ii)],dep_grid,interp_method,'extrap');
        end
    end

elseif vst_option == 1 || vst_option == 2 % regress EOF vertical structure for given Tg points
    % computing EOF pattern
    vstruc_eof = load(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/vg_1100S_' nmid '.mat']);
    delta_z = vstruc_eof.dep_grid(2) - vstruc_eof.dep_grid(1);

    if vst_option == 1 % set level of no motion at reference_level
        vstruc_eof.Tg_basin = vstruc_eof.Tg_basin - vstruc_eof.Tg_basin(vstruc_eof.dep_grid==reference_level,:);
    end
    idep_use = find(vstruc_eof.dep_grid<1400); % define depth range over which EOF is computed (must be larger than reference_level)
    vstruc_eof.Tg4EOF = vstruc_eof.Tg_basin(idep_use,:)' - mean(vstruc_eof.Tg_basin(idep_use,:),2)';
    vstruc_eof.C = cov(vstruc_eof.Tg4EOF);
    [vstruc_eof.EOFs,vstruc_eof.D] = eigs(vstruc_eof.C);

    % explained variance by EOF pattern
    vstruc_eof.lambda = diag(vstruc_eof.D);
    vstruc_eof.exp_var = vstruc_eof.lambda ./ sum(vstruc_eof.lambda,'omitnan');
    disp('Explained variance of EOFs in %:')
    disp(vstruc_eof.exp_var.*100)
    
    % regressing EOF pattern
    nb_EOF_use = 2; % standard is 2, can be changed for tests
    if stp_option == 1 % BPRs at 300m, 500m
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp300,-1) | vstruc_eof.dep_grid == round(depth_bp500,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_300m,Tg_500m);
    elseif stp_option == 2 % BPRs at 500m
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp500,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_500m);
    elseif stp_option == 3 % BPRs at 300m
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp300,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_300m);
    elseif stp_option == 4 % BPRs at 300m, 500m, 1200m
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp300,-1) |  vstruc_eof.dep_grid == round(depth_bp500,-1) | vstruc_eof.dep_grid == round(depth_bp1200,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_300m,Tg_500m,Tg_1200m);
    elseif stp_option == 5 % BPRs at 500m, 1200m
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp500,-1) | vstruc_eof.dep_grid == round(depth_bp1200,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_500m,Tg_1200m);
    elseif stp_option == 6 % BPRs at 1200m
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp1200,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_1200m);
    end
    obs_use = obs_use - mean(obs_use,2,'omitnan');
    alpha = (EOF_use' * EOF_use)\(EOF_use' * obs_use);
    Tg_basin = vstruc_eof.EOFs(vstruc_eof.dep_grid <= reference_level,1:nb_EOF_use) * alpha;
    dep_grid = vstruc_eof.dep_grid(vstruc_eof.dep_grid <= reference_level);

elseif vst_option == 5 % combination of EOF regression and interpolation for Tg
    % computing EOF pattern
    vstruc_eof = load(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/vg_1100S_' nmid '.mat']);
    delta_z = vstruc_eof.dep_grid(2) - vstruc_eof.dep_grid(1);
    vstruc_eof.Tg_basin = vstruc_eof.Tg_basin - vstruc_eof.Tg_basin(vstruc_eof.dep_grid==reference_level,:); % referencing to level of no net motion
    idep_use = find(vstruc_eof.dep_grid<500);
    vstruc_eof.Tg4EOF = vstruc_eof.Tg_basin(idep_use,:)' - mean(vstruc_eof.Tg_basin(idep_use,:),2)';
    vstruc_eof.C = cov(vstruc_eof.Tg4EOF);
    [vstruc_eof.EOFs,~] = eigs(vstruc_eof.C);

    % refressing EOF pattern until 500m BPR
    if stp_option == 4 % BPRs at 300m, 500m, 1200m
        nb_EOF_use = 3;
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp300,-1) |  vstruc_eof.dep_grid == round(depth_bp500,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_300m,Tg_500m);
    elseif stp_option == 5 % BPRs at 500m, 1200m
        nb_EOF_use = 2;
        EOF_use = vstruc_eof.EOFs(vstruc_eof.dep_grid == 0 | vstruc_eof.dep_grid == round(depth_bp500,-1),1:nb_EOF_use);
        obs_use = cat(1,Tg_0m,Tg_500m);
    end
    obs_use = obs_use - mean(obs_use,2,'omitnan');
    alpha = (EOF_use' * EOF_use)\(EOF_use' * obs_use);
    Tg_basin_upper = vstruc_eof.EOFs(vstruc_eof.dep_grid <= round(depth_bp500,-1),1:nb_EOF_use) * alpha;
    dep_grid = vstruc_eof.dep_grid(vstruc_eof.dep_grid <= reference_level);

    % linear interpolation below 500m BPR
    Tg_basin_lower = ones(length(dep_grid(dep_grid > round(depth_bp500,-1))),length(time_mly)).*NaN;
    for ii = 1:length(time_mly)
        Tg_basin_lower(:,ii) = interp1([round(depth_bp500,-1) depth_bp1200],[Tg_basin_upper(end,ii) Tg_1200m(ii)],dep_grid(dep_grid > round(depth_bp500,-1)),'linear','extrap');
    end
    Tg_basin = cat(1,Tg_basin_upper,Tg_basin_lower);
end

% vertical integration for geostrophic AMOC transport
AMOCg = sum(Tg_basin.*delta_z,1).*1e-6;
    

% --- save variables in interim folder, then download to plot (then script plot_BPRSetup)
if ebb_option == 2 || ebb_option == 3
    vars2save = {'bpa_eb300_mly_orig','bpa_eb500_mly_orig'};
end
if prc_option == 4
    if exist('vars2save','var')
        vars2save = {vars2save{:},'yyears','bp_wb300_reconstr_yly','bp_wb500_reconstr_yly','bp_wb1200_reconstr_yly','bp_wb300_yly','bp_wb500_yly','bp_wb1200_yly'};
    else
        vars2save = {'yyears','bp_wb300_reconstr_yly','bp_wb500_reconstr_yly','bp_wb1200_reconstr_yly','bp_wb300_yly','bp_wb500_yly','bp_wb1200_yly'};
    end
elseif prc_option == 5 || prc_option == 6
    if exist('vars2save','var')
        vars2save = {vars2save{:},'bpa_wb300_mly_orig','bpa_wb500_mly_orig','bpa_wb1200_mly_orig'};
    else
        vars2save = {'bpa_wb300_mly_orig','bpa_wb500_mly_orig','bpa_wb1200_mly_orig'};
    end
end

if exist('vars2save','var')
    save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/bpr_pos' num2str(posbpr_option) num2str(possla_option) '_prc' num2str(prc_option) '_ebb' num2str(ebb_option) '_ref' num2str(ref_option) '_vst' num2str(vst_option) '_stp' num2str(stp_option) '_' nmid test_suffix '.mat'],'time_mly','dep_grid','Tg_basin','AMOCg','bpa_eb300_mly','bpa_eb500_mly','bpa_eb1200_mly','bpa_wb300_mly','bpa_wb500_mly','bpa_wb1200_mly',vars2save{:},'-v7.3')
else
    save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/bpr_pos' num2str(posbpr_option) num2str(possla_option) '_prc' num2str(prc_option) '_ebb' num2str(ebb_option) '_ref' num2str(ref_option) '_vst' num2str(vst_option) '_stp' num2str(stp_option) '_' nmid test_suffix '.mat'],'time_mly','dep_grid','Tg_basin','AMOCg','bpa_eb300_mly','bpa_eb500_mly','bpa_eb1200_mly','bpa_wb300_mly','bpa_wb500_mly','bpa_wb1200_mly','-v7.3')
end

end



% --- defining functions used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bp_xxx_reconstr_yly,bp_xxx_yly,yyears] = bp_from_ies_yly(dep_xxx,ptemp_xxx,salt_xxx,bp_xxx,dep_xxx_grid_ies,ilat_xxx, lat,time,delta_z_ies)
% to reconstruct yearly bottom pressure data from ies travel time
    
    % interpolating T&S on 10m grid 
    pres_xxx_grid_ies = sw_pres(dep_xxx_grid_ies,lat(ilat_xxx));

    ptemp_xxx_grid_ies = interp1(dep_xxx,ptemp_xxx,dep_xxx_grid_ies,'linear','extrap');
    ptemp_xxx_grid_ies(1,:) = ptemp_xxx_grid_ies(2,:); % assuming a mixed layer in the upper 10m
    salt_xxx_grid_ies = interp1(dep_xxx,salt_xxx,dep_xxx_grid_ies,'linear','extrap');
    salt_xxx_grid_ies(1,:) = salt_xxx_grid_ies(2,:);
        
    % computing sound velocity (in m/s) for model bins
    temp_xxx_grid_ies = sw_temp(salt_xxx_grid_ies,ptemp_xxx_grid_ies,pres_xxx_grid_ies',0);
    svel_xxx_grid_ies = sw_svel(salt_xxx_grid_ies, temp_xxx_grid_ies, pres_xxx_grid_ies');
      
    % deriving round trip travel time
    tau_xxx = 2.* (sum(delta_z_ies ./ svel_xxx_grid_ies(2:end-1,:),1) + sum(0.5*delta_z_ies ./ svel_xxx_grid_ies([1 end],:),1));

    % fit of travel time on annual bottom pressure data
    % to omit noise from seasonal cycle --> considering annual mean data
    [yy,~,~] = datevec(time);
    yyears = unique(yy);
    for i = 1:length(yyears)
        tau_xxx_yly(i) = mean(tau_xxx(yy==yyears(i)));
        bp_xxx_yly(i) = mean(bp_xxx(yy==yyears(i)));
    end
            
    [p,S] = polyfit(tau_xxx_yly,bp_xxx_yly,3); % cubic polynomial
    [bp_xxx_reconstr_yly, ~] = polyval(p,tau_xxx_yly,S);
    rmse(bp_xxx_yly,bp_xxx_reconstr_yly)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bp_xxx_reconstr_dly] = bp_from_ies_yly2dly(bp_xxx_reconstr_yly,time,yyears)
% to interpolate yearly to monthly data while conserving the annual mean
            bp_xxx_reconstr_dly_unscaled = interp1(datenum(yyears,7,1), bp_xxx_reconstr_yly, time, 'spline'); % extrapolates first and last 6 months
            bp_xxx_reconstr_dly = bp_xxx_reconstr_dly_unscaled;
            [yy,~,~] = datevec(time);
            for k = 1:length(yyears)
                target = bp_xxx_reconstr_yly(k);
                current_mean = mean(bp_xxx_reconstr_dly_unscaled(yy == yyears(k)));
                if current_mean ~= 0
                    factor = target / current_mean;
                    bp_xxx_reconstr_dly(yy == yyears(k)) = bp_xxx_reconstr_dly_unscaled(yy == yyears(k)) * factor;
                else
                    % fallback: set constant months equal to target (avoids divide by zero)
                    bp_xxx_reconstr_dly(yy == yyears(k)) = target;
                end
            end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bpa_xxx_mly] = bp_from_WB4(waterdep_xxx,dep_WB4,ssh_WB4_mly,v_WB_mly,nmid,g,f,rho0,x_spacing,prc_option)
% to reconstruct WB bottom pressure from the moorings

    % 1. import interpolated WB4 density profile
    if prc_option == 5
        rho_WB4 = load(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/moor_density_wbr2ebr0_vst1eof0sur1_idr1_' nmid '.mat']); % for current setup
    elseif prc_option == 6
        rho_WB4 = load(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/moor_density_wbr3ebr0_vst1eof0sur1_idr1_' nmid '.mat']); % for enhanced setup
    end
    
    % 2. compute WB4 pressure profile
    idep_xxx = round(waterdep_xxx/10)+1;
    bp_xxx_WB4 = (g .*  trapz(rho_WB4.dep_grid(1:idep_xxx),rho_WB4.wb_rho_grid(1:idep_xxx,:),1) + g.*rho0 .* ssh_WB4_mly) .* 1e-4; % in dbar
    
    % 3. interpolate velocities
    [~,idep_xxx] = min(abs(dep_WB4 - waterdep_xxx));
    bp_xxx_vdx = f .* rho0 .* squeeze(sum(v_WB_mly(:,idep_xxx,:),1,'omitnan')) .* x_spacing .* 1e-4; % in dbar
    
    % 4. add contributions from 2 and 3
    bpa_xxx_mly = bp_xxx_WB4 - bp_xxx_vdx';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yout = anharm_bp(x1,x2,y,n)
% to compute annual harmonics
% inputs:
%   x1 = time vector corresponding to y
%   x2 = time vector for desired yout
%   y = data vector to fit harmonics to
%   n = 1 annual harmonics only, n = 2 annual and semi-annual harmonics
% outputs:
%   yout = time series of fitted harmonics

    if nargin < 3
      n = 1;
    end
    
    year = 365.25;
    f = 2*pi*[1:n]/year;

    % --- fit harmonics to data
    x1 = x1(:);
    y = y(:);
    
    V = ones(length(x1),2*length(f)+1);
    for j = 1:length(f)
      V(:,2*j) = cos(f(j)*x1);
      V(:,2*j+1) = sin(f(j)*x1);
    end
    
    p = V\y;

    % --- evaluate harmonics
    d = size(x2);
    x2 = reshape(x2,prod(d),1);

    W = ones(length(x2),length(f));
    for j = 1:length(f)
      W(:,2*j) = cos(f(j)*x2);
      W(:,2*j+1) = sin(f(j)*x2);
    end
    
    yout = reshape(W*p,d);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,yfit] = exp_lin_fit_bp(x,y,a)
% for exponential linear fit:
% yfit = a1*(1-exp(-a2*(x-x(1)))) + a3*(x-x(1)) + a4 of some function y = y(x) using gaussian method
% returns coefficients of model a = [a(1) a(2) a(3) a(4)];
%
% input   y    -- function to which fit will be applied 
%         x    -- y = y(x);
%         a    -- a = [a1 a2 a3 a4] first guess 
%                 for coefficients of fit, default = [0 0 0 0]
%
% output  a    -- a = [a(1) a(2) a(3) a(4)]: coefficients of fit 
%         yfit -- yfit = a1*(1-exp(-a2*(x-x(1)))) + a3*(x-x(1)) + a4

    % --- Model ----------- 
    x0 = x(1);
    funca = @(a,x)(a(1)*(1-exp(-a(2)*(x-x0))) + a(3)*(x-x0) + a(4));
    
    % --- default input values --------
    if nargin < 3 || isempty(a)
        a1 = 0;
        a2 = 0; 
        a3 = 0; 
        a4 = 0; 
    else 
        a1 = a(1);
        a2 = a(2);
        a3 = a(3);
        a4 = a(3);
    end
    
    % ---  simplex - least square fit ----
    a = nlinfit_old(x,y,funca,[a1 a2 a3 a4]);
    yfit = a(1)*(1-exp(-a(2)*(x-x0))) + a(3)*(x-x0) + a(4);

end

function [beta,r,J] = nlinfit_old(X,y,model,beta0)
%NLINFIT Nonlinear least-squares data fitting by the Gauss-Newton method.
%   NLINFIT(X,Y,'MODEL',BETA0) finds the coefficients of the nonlinear 
%   function described in MODEL. MODEL is a user supplied function having 
%   the form y = f(beta,x). That is MODEL returns the predicted values of y
%   given initial parameter estimates, beta, and the independent variable, X.   
%   [BETA,R,J] = NLINFIT(X,Y,'MODEL',BETA0) returns the fitted coefficients
%   BETA the residuals, R, and the Jacobian, J, for use with NLINTOOL to
%   produce error estimates on predictions.

%   B.A. Jones 12-06-94.
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.12 $  $Date: 1998/09/09 19:39:39 $

    n = length(y);
    if min(size(y)) ~= 1
       error('Requires a vector second input argument.');
    end
    y = y(:);
    
    if size(X,1) == 1 % turn a row vector into a column vector.
       X = X(:);
    end
    
    p = length(beta0);
    beta0 = beta0(:);
    
    J = zeros(n,p);
    beta = beta0;
    betanew = beta + 1;
    maxiter = 100;
    iter = 0;
    betatol = 1.0E-4;
    rtol = 1.0E-4;
    sse = 1;
    sseold = sse;
    seps = sqrt(eps);
    zbeta = zeros(size(beta));
    s10 = sqrt(10);
    eyep = eye(p);
    zerosp = zeros(p,1);
    
    while (norm((betanew-beta)./(beta+seps)) > betatol | abs(sseold-sse)/(sse+seps) > rtol) & iter < maxiter
       if iter > 0 
          beta = betanew;
       end
    
       iter = iter + 1;
       yfit = feval(model,beta,X);
       r = y - yfit;
       sseold = r'*r;
    
       for k = 1:p
          delta = zbeta;
          delta(k) = seps*beta(k);
          yplus = feval(model,beta+delta,X);
          J(:,k) = (yplus - yfit)/delta(k);
       end
    
       Jplus = [J;(1.0E-2)*eyep];
       rplus = [r;zerosp];
    
       % Levenberg-Marquardt type adjustment 
       % Gauss-Newton step -> J\r
       % LM step -> inv(J'*J+constant*eye(p))*J'*r
       step = Jplus\rplus;
       
       betanew = beta + step;
       yfitnew = feval(model,betanew,X);
       rnew = y - yfitnew;
       sse = rnew'*rnew;
       iter1 = 0;
       while sse > sseold & iter1 < 12
          step = step/s10;
          betanew = beta + step;
          yfitnew = feval(model,betanew,X);
          rnew = y - yfitnew;
          sse = rnew'*rnew;
          iter1 = iter1 + 1;
       end
    end
    if iter == maxiter
       disp('NLINFIT did NOT converge. Returning results from last iteration.');
    end
end
