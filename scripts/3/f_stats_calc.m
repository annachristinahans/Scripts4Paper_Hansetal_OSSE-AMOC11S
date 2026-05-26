function [stats_tracos, stats_ideal, colcoding] = f_stats_calc(time_mly,tracos_AMOCg,ideal_AMOCg)

% function to compute statistics to evaluate performance of AMOCg reconstruction

% input:
% time_mly: time vector at mly resolution
% tracos_AMOCg: recosntructed AMOCg at mly resolution
% ideal_AMOCg: model-true AMOCg at mly resolution

% output:
% stats_tracos : statistics for the performance of the tested AMOCg reconstruction
% stats_ideal: statistics of the model-truth for comparison
% colcoding: colour coding level used to categorize performance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. computation of different time scales
% detrending
tracos_AMOCg_dtr = detrend(tracos_AMOCg,'linear');
ideal_AMOCg_dtr = detrend(ideal_AMOCg,'linear');
% filtering
Fs = 1/(30*24*3600); % Sampling frequency in Hz (monthly)
Fc_lp = 1/(5*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(2, Fc_lp/(Fs/2), 'low'); % 2nd-order low-pass filter
tracos_AMOCg_lp = filtfilt(b, a, tracos_AMOCg_dtr);
ideal_AMOCg_lp = filtfilt(b, a, ideal_AMOCg_dtr);
Fc_hp = 1/(2*365*24*3600); % Cutoff frequency in Hz
[b, a] = butter(1, Fc_hp/(Fs/2), 'high'); % 1st-order high-pass filter
tracos_AMOCg_hp = filtfilt(b, a, tracos_AMOCg_dtr);
ideal_AMOCg_hp = filtfilt(b, a, ideal_AMOCg_dtr);
% cutting half the filter length in the beginning and end
tracos_AMOCg_lp = tracos_AMOCg_lp(31:end-30);
tracos_AMOCg_hp = tracos_AMOCg_hp(13:end-12);
ideal_AMOCg_lp = ideal_AMOCg_lp(31:end-30);
ideal_AMOCg_hp = ideal_AMOCg_hp(13:end-12);
% rmse
rmse_fs = sqrt(sum((tracos_AMOCg_dtr - ideal_AMOCg_dtr).^2) ./length(ideal_AMOCg_dtr));
rmse_lp = sqrt(sum((tracos_AMOCg_lp - ideal_AMOCg_lp).^2) ./length(ideal_AMOCg_lp));
rmse_hp = sqrt(sum((tracos_AMOCg_hp - ideal_AMOCg_hp).^2) ./length(ideal_AMOCg_hp));
% Pearson correlation coefficient
[r,~] = corrcoef(tracos_AMOCg_dtr,ideal_AMOCg_dtr);
r_fs = r(1,2);
[r,~] = corrcoef(tracos_AMOCg_lp,ideal_AMOCg_lp);
r_lp = r(1,2);
[r,~] = corrcoef(tracos_AMOCg_hp,ideal_AMOCg_hp);
r_hp = r(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Evaluating annual cycle
% by mean year (=climatology)
[~,mm,~] = datevec(time_mly);
for i = 1:12
    ideal_AMOCg_clim(i) = mean(ideal_AMOCg_dtr(mm == i));
    tracos_AMOCg_clim(i) = mean(tracos_AMOCg_dtr(mm == i));
end
ideal_climamp = (max(ideal_AMOCg_clim)-min(ideal_AMOCg_clim));
tracos_climamp = (max(tracos_AMOCg_clim)-min(tracos_AMOCg_clim));
[~, ideal_climpha] = max(ideal_AMOCg_clim);
[~, tracos_climpha] = max(tracos_AMOCg_clim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. linear slope since 1994 (as a measure for the trend)
ideal_mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), ideal_AMOCg(time_mly>=datenum(1994,1,1)));
ideal_slope_per_decade = ideal_mdl.Coefficients.Estimate(2) * 365.25 * 10;

tracos_mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), tracos_AMOCg(time_mly>=datenum(1994,1,1)));
tracos_slope_per_decade = tracos_mdl.Coefficients.Estimate(2) * 365.25 * 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put together the output variables
stats_tracos.rmse_fs = rmse_fs;
stats_tracos.r_fs = r_fs;
stats_tracos.rmse_hp = rmse_hp;
stats_tracos.r_hp = r_hp;
stats_tracos.climamp = tracos_climamp;
stats_tracos.climpha = tracos_climpha;
stats_tracos.rmse_lp = rmse_lp;
stats_tracos.r_lp = r_lp;
stats_tracos.slope_per_decade = tracos_slope_per_decade;

stats_ideal.std_fs = std(ideal_AMOCg_dtr);
stats_ideal.std_hp = std(ideal_AMOCg_hp);
stats_ideal.climamp = ideal_climamp;
stats_ideal.climpha = ideal_climpha;
stats_ideal.std_lp = std(ideal_AMOCg_lp);
stats_ideal.slope_per_decade = ideal_slope_per_decade;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Colour-coding

% choose stats to reference to for colour-coding
stats_ref = stats_ideal;

% evaluate stats compared to reference
error_climamp = abs(stats_ref.climamp - stats_tracos.climamp) / stats_ref.climamp;
error_climpha = abs(stats_ref.climpha - stats_tracos.climpha);

error_rmse_fs = stats_tracos.rmse_fs / stats_ref.std_fs;
error_rmse_hp = stats_tracos.rmse_hp / stats_ref.std_hp;
error_rmse_lp = stats_tracos.rmse_lp / stats_ref.std_lp;

% match colours
limits_r = [0.7 0.85]; % 0.7 -> about 50% of the variance explained, 0.85 -> about 75% of variance explained
limits_rmse = [1 0.5];
limits_climpha = [2 1];
limits_climamp = [0.4 0.2];

% for rmse_fs
if error_rmse_fs > limits_rmse(1)
    colcoding.rmse_fs = 0;
elseif error_rmse_fs > limits_rmse(2)
    colcoding.rmse_fs = 1;
else
    colcoding.rmse_fs = 2;
end

% for r_fs
if stats_tracos.r_fs < limits_r(1)
    colcoding.r_fs = 0;
elseif stats_tracos.r_fs < limits_r(2)
    colcoding.r_fs = 1;
else
    colcoding.r_fs = 2;
end

% for rmse_hp
if error_rmse_hp > limits_rmse(1)
    colcoding.rmse_hp = 0;
elseif error_rmse_hp > limits_rmse(2)
    colcoding.rmse_hp = 1;
else
    colcoding.rmse_hp = 2;
end

% for r_hp
if stats_tracos.r_hp < limits_r(1)
    colcoding.r_hp = 0;
elseif stats_tracos.r_hp < limits_r(2)
    colcoding.r_hp = 1;
else
    colcoding.r_hp = 2;
end

% for climamp
if error_climamp >= limits_climamp(1)
    colcoding.climamp = 0;
elseif error_climamp >= limits_climamp(2)
    colcoding.climamp = 1;
else
    colcoding.climamp = 2;
end

% for climpha
if error_climpha >= limits_climpha(1)
    colcoding.climpha = 0;
elseif error_climpha >= limits_climpha(2)
    colcoding.climpha = 1;
else
    colcoding.climpha = 2;
end

% for rmse_lp
if error_rmse_lp > limits_rmse(1)
    colcoding.rmse_lp = 0;
elseif error_rmse_lp > limits_rmse(2)
    colcoding.rmse_lp = 1;
else
    colcoding.rmse_lp = 2;
end

% for r_lp
if stats_tracos.r_lp < limits_r(1)
    colcoding.r_lp = 0;
elseif stats_tracos.r_lp < limits_r(2)
    colcoding.r_lp = 1;
else
    colcoding.r_lp = 2;
end

% for slope_per_decade
% whether the trend found is not significantly different from trend of model-truth

% approach 1:
mdl = fitlm(time_mly(time_mly>=datenum(1994,1,1)), tracos_AMOCg(time_mly>=datenum(1994,1,1))-ideal_AMOCg(time_mly>=datenum(1994,1,1)));
if mdl.Coefficients.pValue(2) < 0.05 % -> slopes are significantly different
    colcoding.slope_per_decade = 0;
else % -> no evidence slopes differ
    colcoding.slope_per_decade = 1;
end

% approach 2:
% x = [time_mly(time_mly>=datenum(1994,1,1))'; time_mly(time_mly>=datenum(1994,1,1))'];
% y = [ideal_AMOCg(time_mly>=datenum(1994,1,1))'; tracos_AMOCg(time_mly>=datenum(1994,1,1))'];
% group = [zeros(length(time_mly(time_mly>=datenum(1994,1,1))),1); ones(length(time_mly(time_mly>=datenum(1994,1,1))),1)];
% tbl = table(x, group, y);
% mdl = fitlm(tbl, 'y ~ x + group + x:group'); % model: Predict y using x, group, and their interaction; mathematically: y=β0​+β1​x+β2​group+β3​(x⋅group)+ε and slopes equal if β3​=0
% if mdl.Coefficients.pValue('x:group') < 0.05 % -> slopes are significantly different
%     colcoding.slope_per_decade = 0;
% else % -> no evidence slopes differ
%     colcoding.slope_per_decade = 1;
% end
