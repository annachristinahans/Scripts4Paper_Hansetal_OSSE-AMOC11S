%% RUN function to compute statistics for INALT/VIKING BPR AMOCg transport

% define combinations of options to be run
prc_option = [0 1 2 31 32 33 34 0 1 0 0 0 0 0 0 0 0 0 0 0 0 5 6];
ebb_option =     [0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 2 0 0];
ref_option =     [1 1 1 1 1 1 1 1 1 2 3 1 1 1 1 1 1 1 1 1 1 1 1];
vst_option =     [1 1 1 1 1 1 1 1 1 1 1 2 3 4 1 1 3 2 5 5 1 1 1];
stp_option =     [1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 4 4 4 5 1 1 1];

% prepare structure to save lines
nb_stats_ideal = 6;
nb_stats_tracos = 9;
stats_tracos_all = ones(length(prc_option),nb_stats_tracos) .* NaN;
colcoding_all = ones(length(prc_option),nb_stats_tracos) .* NaN;
stats_ideal_all = ones(1,nb_stats_ideal) .* NaN;

% stats for BPR Setup for all defined combinations
for i = 1:length(prc_option)
    [stats_tracos, stats_ideal, colcoding] = f_stats_BPRSetup(2,2,2,prc_option(i),ebb_option(i),ref_option(i),vst_option(i),stp_option(i),'');

    if i == 1
        fields_ideal = fieldnames(stats_ideal);
        for j = 1:nb_stats_ideal
            stats_ideal_all(1,j) = stats_ideal.(fields_ideal{j});
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

    disp([num2str(i) ' out of ' num2str(length(prc_option)) ' done'])
end

% saving options
options = [1:length(prc_option); prc_option; ebb_option; ref_option; vst_option; stp_option]';
writematrix(options,'/Users/ahans/Documents/PhD/PLOTS/Model/BPRSetup/table/stats_BPRSetup_options_VIK_v4.txt','Delimiter','\t');

% saving stats tracos
fid = fopen('/Users/ahans/Documents/PhD/PLOTS/Model/BPRSetup/table/stats_BPRSetup_tracos_VIK_v4.txt','w');
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
fid = fopen('/Users/ahans/Documents/PhD/PLOTS/Model/BPRSetup/table/stats_BPRSetup_modeltruth_VIK_v4.txt','w');
% write header
fprintf(fid,'%s\t', fields_ideal{1:end-1});
fprintf(fid,'%s\n', fields_ideal{end});
% write data
fprintf(fid, '%.2f\t', stats_ideal_all);
fclose(fid);

% saving colcoding
fid = fopen('/Users/ahans/Documents/PhD/PLOTS/Model/BPRSetup/table/stats_BPRSetup_colcoding_VIK_v4.txt','w');
% write header
fprintf(fid,'%s\t', fields_colcoding{1:end-1});
fprintf(fid,'%s\n', fields_colcoding{end});
% write data row by row
for i = 1:size(stats_tracos_all,1)
    fprintf(fid, '%.0f\t', colcoding_all(i,1:end-1));
    fprintf(fid, '%.0f\n', colcoding_all(i,end));
end
fclose(fid);
