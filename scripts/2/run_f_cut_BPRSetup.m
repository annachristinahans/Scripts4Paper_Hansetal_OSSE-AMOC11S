%% RUN function to cut INALT/VIKING BPR AMOCg transport on nesh

% define combinations of options to be run
prc_option = [0 1 2 31 32 33 34 0 1 0 0 0 0 0 0 0 0 0 0 0 0 5 6 1 4];
ebb_option =     [0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 2 0 0 3 0];
ref_option =     [1 1 1 1 1 1 1 1 1 2 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
vst_option =     [1 1 1 1 1 1 1 1 1 1 1 2 3 4 1 1 3 2 5 5 1 1 1 1 1];
stp_option =     [1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 4 4 4 5 1 1 1 1 1];

% cut BPR Setup for all defined combinations
for i = 1:length(prc_option)
    f_cut_BPRSetup(2,2,2,prc_option(i),ebb_option(i),ref_option(i),vst_option(i),stp_option(i),'')
    disp([num2str(i) ' out of ' num2str(length(prc_option)) ' done']) 
end
