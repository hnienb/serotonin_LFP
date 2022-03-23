function plot_c2s
%%
% plot results from spike prediction from LFPs
%

% Use the cleaned version with 68 sessions
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-2), '/');
data_path = [dir_path, '/resources/Data/'];
addpath(data_path)

if nargin < 1
    load([data_path, '/c2s/met_cv10.mat'], 'met');
end


metrics = {'correlations', 'info'};
metricnames = {'correlation', 'MI'};
lenm = length(metrics);
drugs = {'NaCl', '5HT'};

figure;
for k = 1:lenm % metric
    for d = 1:2 % NaCl or 5HT
        subplot(lenm, 3, d+3*(k-1))
        x = met(met(:, 2)==d-1, 3+2*(k-1));
        y = met(met(:, 2)==d-1, 4+2*(k-1));
        z = met(met(:, 2)==d-1, 1);
        unity_scatter(x, y, z)
        hold on;
        if k==1
            title([drugs{d}])
        end
        if d==1
            ylabel({metricnames{k}, 'drug'})
        else
            xlabel('baseline')
        end
    end
    
    % FR control
    subplot(lenm, 3, 3+3*(k-1))
    x = met(:, 7+2*(k-1));
    y = met(:, 8+2*(k-1));
    z = met(:, 1);
    unity_scatter(x, y, z)
    hold on;
    if k==1
        title('FR control')
    else
        xlabel('high FR')
    end
    ylabel('low FR')
end
    

% Collect the data for correlations and information for each monkey and
% condition
kaki_nacl_corr = met(met(:,1)==0 & met(:,2)==0, 3:4);
kaki_5HT_corr = met(met(:,1)==0 & met(:,2)==1, 3:4);
kaki_nacl_info = met(met(:,1)==0 & met(:,2)==0, 5:6);
kaki_5HT_info = met(met(:,1)==0 & met(:,2)==1, 5:6);

mango_nacl_corr = met(met(:,1)==1 & met(:,2)==0, 3:4);
mango_5HT_corr = met(met(:,1)==1 & met(:,2)==1, 3:4);
mango_nacl_info = met(met(:,1)==1 & met(:,2)==0, 5:6);
mango_5HT_info = met(met(:,1)==1 & met(:,2)==1, 5:6);


% Perform statistical tests for Correlation data for each monkey (and
% together)
disp('Stats for Correlation Data')
stats_kaki_corr = pair_tests(kaki_nacl_corr, kaki_5HT_corr);
disp(['Kaki (NaCl; n=' num2str(size(kaki_nacl_corr, 1)) ') ------------------------'])
disp(stats_kaki_corr.pair(1).table)
disp(['Kaki (5HT; n=' num2str(size(kaki_5HT_corr, 1)) ') ------------------------'])
disp(stats_kaki_corr.pair(2).table)
disp(['Kaki (NaCl vs 5HT; n=' num2str(size(kaki_5HT_corr, 1)+size(kaki_nacl_corr, 1)) ') ---------'])
disp(stats_kaki_corr.table)

stats_mango_corr = pair_tests(mango_nacl_corr, mango_5HT_corr);
disp(['Mango (NaCl; n=' num2str(size(mango_nacl_corr, 1)) ') ------------------------'])
disp(stats_mango_corr.pair(1).table)
disp(['Mango (5HT; n=' num2str(size(mango_5HT_corr, 1)) ') ------------------------'])
disp(stats_mango_corr.pair(2).table)
disp(['Mango (NaCl vs 5HT; n=' num2str(size(mango_5HT_corr, 1)+size(mango_nacl_corr, 1)) ') ---------'])
disp(stats_mango_corr.table)

stats_all_corr = pair_tests([kaki_nacl_corr; mango_nacl_corr], [kaki_5HT_corr; mango_5HT_corr]);
disp(['Both (NaCl; n=' num2str(size(mango_nacl_corr, 1)+size(kaki_nacl_corr, 1)) ') ------------------------'])
disp(stats_all_corr.pair(1).table)
disp(['Both (5HT; n=' num2str(size(mango_5HT_corr, 1)+size(kaki_5HT_corr, 1)) ') ------------------------'])
disp(stats_all_corr.pair(2).table)
disp(['Both (NaCl vs 5HT; n=' num2str(size(mango_nacl_corr, 1)+size(mango_5HT_corr, 1)+size(kaki_nacl_corr, 1)+size(kaki_5HT_corr, 1)) ') ---------'])
disp(stats_all_corr.table)


% Perform statistical tests for Information data for each monkey (and
% together)
disp("Stats for Mutual Information Data")
stats_kaki_info = pair_tests(kaki_nacl_info, kaki_5HT_info);
disp(['Kaki (NaCl; n=' num2str(size(kaki_nacl_info, 1)) ') ------------------------'])
disp(stats_kaki_info.pair(1).table)
disp(['Kaki (5HT; n=' num2str(size(kaki_5HT_info, 1)) ') ------------------------'])
disp(stats_kaki_info.pair(2).table)
disp(['Kaki (NaCl vs 5HT; n=' num2str(size(kaki_5HT_info, 1)+size(kaki_nacl_info, 1)) ') ---------'])
disp(stats_kaki_info.table)

stats_mango_info = pair_tests(mango_nacl_info, mango_5HT_info);
disp(['Mango (NaCl; n=' num2str(size(mango_nacl_info, 1)) ') ------------------------'])
disp(stats_mango_info.pair(1).table)
disp(['Mango (5HT; n=' num2str(size(mango_5HT_info, 1)) ') ------------------------'])
disp(stats_mango_info.pair(2).table)
disp(['Mango (NaCl vs 5HT; n=' num2str(size(mango_5HT_info, 1)+size(mango_nacl_info, 1)) ') ---------'])
disp(stats_mango_info.table)

stats_all_info = pair_tests([kaki_nacl_info; mango_nacl_info], [kaki_5HT_info; mango_5HT_info]);
disp(['Both (NaCl; n=' num2str(size(mango_nacl_info, 1)+size(kaki_nacl_info, 1)) ') ------------------------'])
disp(stats_all_info.pair(1).table)
disp(['Both (5HT; n=' num2str(size(mango_5HT_info, 1)+size(kaki_5HT_info, 1)) ') ------------------------'])
disp(stats_all_info.pair(2).table)
disp(['Both (NaCl vs 5HT; n=' num2str(size(mango_nacl_info, 1)+size(mango_5HT_info, 1)+size(kaki_nacl_info, 1)+size(kaki_5HT_info, 1)) ') ---------'])
disp(stats_all_info.table)
    
    
% Collect the FR data for correlations and information for each monkey
kaki_corr = met(met(:,1)==0, 7:8);
kaki_info = met(met(:,1)==0, 9:10);

mango_corr = met(met(:,1)==1, 7:8);
mango_info = met(met(:,1)==1, 9:10);

% Perform statistical tests for Correlation data for each monkey (and
% together)
disp('Statistics for FR Control - Correlation data:')
stats_kaki_corr = pair_tests(kaki_corr);
disp(['Kaki (n=' num2str(size(kaki_corr, 1)) ') ------------------------'])
disp(stats_kaki_corr.pair(1).table)

stats_mango_corr = pair_tests(mango_corr);
disp(['Mango (n=' num2str(size(mango_corr, 1)) ') ------------------------'])
disp(stats_mango_corr.pair(1).table)

stats_all_corr = pair_tests([kaki_corr; mango_corr]);
disp(['Both (n=' num2str(size([kaki_corr;mango_corr], 1)) ') ------------------------'])
disp(stats_all_corr.pair(1).table)

% Perform statistical tests for Information data for each monkey (and
% together)
disp('Statistics for FR Control - Information data:')
stats_kaki_info = pair_tests(kaki_info);
disp(['Kaki (n=' num2str(size(kaki_info, 1)) ') ------------------------'])
disp(stats_kaki_info.pair(1).table)

stats_mango_info = pair_tests(mango_info);
disp(['Mango (n=' num2str(size(mango_info, 1)) ') ------------------------'])
disp(stats_mango_info.pair(1).table)

stats_all_info = pair_tests([kaki_info; mango_info]);
disp(['Both (n=' num2str(size([kaki_info;mango_info], 1)) ') ------------------------'])
disp(stats_all_info.pair(1).table)
end