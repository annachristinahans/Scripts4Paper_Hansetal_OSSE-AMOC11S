%% INALT/VIKING mean TRACOS section - cut
% to run on nesh
% gives topography and time-mean velocity for the TRACOS section

% --- defining the TRACOS section
% define positions of CTD stations to create the section
ctd_wb_lat = [14.20 14.60 15.30 16.00 19.50 22.80 27.40 32.00 36.50 41.40 46.40 51.40 56.40 07.60+60]./60;
ctd_wb_lat = -10 - ctd_wb_lat;
ctd_wb_lon = [54.20 53.60 52.60 51.70 46.10 40.80 34.90 29.30 23.60 17.60 11.60 05.60 59.60-60 43.90-60]./60;
ctd_wb_lon = -35 - ctd_wb_lon;

ctd_eb_lat = [-11.0725 -10.4666];
ctd_eb_lon = [12.6400 13.5300];

% interpolate ctd stations on model lonV grid
lon_wb = -36.5:0.05:-34.7;
lat_wb = interp1(ctd_wb_lon,ctd_wb_lat,lon_wb,'linear','extrap');

lon_eb = 12.5:0.05:14;
lat_eb = interp1(ctd_eb_lon,ctd_eb_lat,lon_eb,'linear','extrap');

% merge to one section where the interior is at 11°S
tracos.lon_tracos = cat(2,lon_wb(1:33),lon_wb(34):0.05:lon_eb(5),lon_eb(6:end));
tracos.lat_tracos = cat(2,lat_wb(1:33),ones(1,length(lon_wb(34):0.05:lon_eb(5)))-12,lat_eb(6:end));
clearvars -except tracos

% --- loading model data
num_yy = 44;
nmid = 'VIK';
mmask = '/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/1_VIKING20X.L46-KFS003_mesh_mask_11S_box.nc';

files1 = dir('/gxfs_work/geomar/smomw628/VIKING20X.L46-KFS003_11S/*V_11S_box.nc');

latV = double(ncread(mmask,'gphiv'));
latV = latV(1,:);
lonV = double(ncread(mmask,'glamv'));
lonV = lonV(:,1)';
depV = double(ncread(mmask,'gdepv')); % in m, including partial cells
dep_cell_size = double(ncread(mmask,'e3v_0')); % in m
qc = double(ncread(mmask,'vmask')); % 0 land, 1 water
qc(qc==0) = NaN;

% cutting the TRACOS section
ilon = find(lonV <= max(tracos.lon_tracos) & lonV >= min(tracos.lon_tracos));
ilat = interp1(latV,1:numel(latV),tracos.lat_tracos,'nearest');
ilon = ilon(~isnan(ilat));
ilat = ilat(~isnan(ilat));
lon_tracos = lonV(ilon);
lat_tracos = latV(ilat);

% defining the topography
for i = 1:length(ilon)
    idep = find(squeeze(qc(ilon(i),ilat(i),:)) == 1,1,'last');
    waterdep_tracos(i) = sum(squeeze(dep_cell_size(ilon(i),ilat(i),1:idep)));
    dep(i,:) = depV(ilon(i),ilat(i),:);
end

% loading meridional velocities
countt = 1;
for i = 1:num_yy % loop over years

    % load all meridional velocities and merge them to one file
    tim = double(ncread([files1(i).folder '/' files1(i).name],'time_counter'));
    time(countt:countt-1+length(tim)) = datenum(1900,1,1,0,0,tim);
    
    v_use = double(ncread([files1(i).folder '/' files1(i).name],'vomecrty')); % in m/s

    % applying land-sea mask
    v_use = v_use .* qc;
    
    for j = 1:length(ilon)
        v(j,:,countt:countt-1+length(tim)) = squeeze(v_use(ilon(j),ilat(j),:,:));
    end
    countt = countt + length(tim);
end

v_mean_tracos = mean(v,3);

save(['/gxfs_work/geomar/smomw628/interim_' nmid '_v4/v_topo_tracos_' nmid '.mat'],'lat_tracos','lon_tracos','waterdep_tracos','dep','v_mean_tracos','-v7.3')
