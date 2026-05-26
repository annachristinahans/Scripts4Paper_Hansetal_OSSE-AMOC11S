%% RUN function to cut INALT/VIKING Mooring AMOCg transport on nesh

%% for table
% define combinations of options to be run
wbr_option = [3 0 0 2 3 2 3 3 2 2 3 3 2 3];
ebr_option = [0 3 4 4 4 2 2 3 1 2 2 3 1 3];
idr_option = [1 1 1 1 1 2 2 2 3 3 3 3 4 4];

% cut Mooring Setup for all defined combinations
for i = 1:length(wbr_option)
    f_cut_MooringSetup(2,wbr_option(i),ebr_option(i),1,0,1,idr_option(i),'')
    disp([num2str(i) ' out of ' num2str(length(wbr_option)) ' done (table)']) 
end

%% for plots
% define combinations of options to be run
wbr_option = [1 1 1 2 2 2 2 3 2 2 2 2 2 2];
ebr_option = [1 1 1 2 2 2 2 3 2 2 2 2 2 2];
vst_option = [1 2 3 1 2 3 3 1 1 1 2 2 3 3];
eof_option = [0 0 1 0 0 1 2 0 0 0 0 0 1 2];
sur_option = [0 0 0 0 0 0 0 0 1 3 1 2 1 1];

% cut Mooring Setup for all defined combinations
for i = 1:length(wbr_option)
    f_cut_MooringSetup(2,wbr_option(i),ebr_option(i),vst_option(i),eof_option(i),sur_option(i),1,'')
    disp([num2str(i) ' out of ' num2str(length(wbr_option)) ' done (plots)']) 
end
