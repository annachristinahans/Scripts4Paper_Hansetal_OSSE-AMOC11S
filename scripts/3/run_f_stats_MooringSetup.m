%% RUN function to compute statistics for INALT/VIKING Mooring AMOCg transport

% define combinations of options to be run
wbr_option = [3 0 0 2 3 2 3 3 2 2 3 3 2 3];
ebr_option = [0 3 4 4 4 2 2 3 1 2 2 3 1 3];
idr_option = [1 1 1 1 1 2 2 2 3 3 3 3 4 4];

% prepare structure to save lines
nb_stats_ideal = 7;
nb_stats_tracos = 10;
stats_tracos_all = ones(length(wbr_option),nb_stats_tracos) .* NaN;
colcoding_all = ones(length(wbr_option),nb_stats_tracos) .* NaN;
stats_ideal_all = ones(1,nb_stats_ideal) .* NaN;
stats_bcs_all = ones(2,nb_stats_tracos) .* NaN;

% stats for BPR Setup for all defined combinations
for i = 1:length(wbr_option)
    [stats_tracos, stats_ideal, stats_bcs, colcoding, colcoding_bcs] = f_stats_MooringSetup(2,wbr_option(i),ebr_option(i),1,0,1,idr_option(i),'');

    if i == 1
        fields_ideal = fieldnames(stats_ideal);
        for j = 1:nb_stats_ideal
            stats_ideal_all(1,j) = stats_ideal.(fields_ideal{j});
        end
        fields_bcs = fieldnames(stats_bcs);
        fields_colcoding_bcs = fieldnames(colcoding_bcs);
        for j = 1:nb_stats_tracos
            stats_bcs_all(1,j) = stats_bcs.(fields_bcs{j});
            stats_bcs_all(2,j) = colcoding_bcs.(fields_colcoding_bcs{j});
        end
    end

    fields_tracos = fieldnames(stats_tracos);
    for j = 1:nb_stats_tracos
        stats_tracos_all(i,j) = stats_tracos.(fields_tracos{j});
    end

    fields_colcoding = fieldnames(colcoding);
    for j = 1:nb_stats_tracos
        colcoding_all(i,j) = colcoding.(fields_colcoding{j});
    end

    disp([num2str(i) ' out of ' num2str(length(wbr_option)) ' done'])
end

% saving options
options = [1:length(wbr_option); wbr_option; ebr_option; idr_option]';
writematrix(options,'/Users/ahans/Documents/PhD/PLOTS/Model/MooringSetup/table/stats_MooringSetup_options_VIK_v4.txt','Delimiter','\t');

% saving stats tracos
fid = fopen('/Users/ahans/Documents/PhD/PLOTS/Model/MooringSetup/table/stats_MooringSetup_tracos_VIK_v4.txt','w');
% write header
fprintf(fid,'%s\t', fields_tracos{1:end-1});
fprintf(fid,'%s\n', fields_tracos{end});
% write data row by row
for i = 1:size(stats_tracos_all,1)
    fprintf(fid, '%.2f\t', stats_tracos_all(i,1:end-1));
    fprintf(fid, '%.2f\n', stats_tracos_all(i,end));
end
fclose(fid);

% saving stats ideal
fid = fopen('/Users/ahans/Documents/PhD/PLOTS/Model/MooringSetup/table/stats_MooringSetup_modeltruth_VIK_v4.txt','w');
% write header
fprintf(fid,'%s\t', fields_ideal{1:end-1});
fprintf(fid,'%s\n', fields_ideal{end});
% write data
fprintf(fid, '%.2f\t', stats_ideal_all);
fclose(fid);

% saving colcoding
fid = fopen('/Users/ahans/Documents/PhD/PLOTS/Model/MooringSetup/table/stats_MooringSetup_colcoding_VIK_v4.txt','w');
% write header
fprintf(fid,'%s\t', fields_colcoding{1:end-1});
fprintf(fid,'%s\n', fields_colcoding{end});
% write data row by row
for i = 1:size(stats_tracos_all,1)
    fprintf(fid, '%.0f\t', colcoding_all(i,1:end-1));
    fprintf(fid, '%.0f\n', colcoding_all(i,end));
end
fclose(fid);

% saving stats boundary currents only
fid = fopen('/Users/ahans/Documents/PhD/PLOTS/Model/MooringSetup/table/stats_MooringSetup_bcs_VIK_v4.txt','w');
% write header
fprintf(fid,'%s\t', fields_bcs{1:end-1});
fprintf(fid,'%s\n', fields_bcs{end});
% write data
for i = 1:size(stats_bcs_all,1)
    fprintf(fid, '%.2f\t', stats_bcs_all(i,1:end-1));
    fprintf(fid, '%.2f\n', stats_bcs_all(i,end));
end
fclose(fid);
