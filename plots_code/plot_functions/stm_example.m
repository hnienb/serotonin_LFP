function stm_example
%%
% visualize example trace and prediction
%

% Use the cleaned version with 68 sessions
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-3), '/');
data_path = [dir_path, '/resources/Data/c2s/data/ka_0278_c1_sortLH_all.grating.ORxRC_5HT/'];
addpath(data_path)

prep = load([data_path '/preprocessed_base_cv10.mat']);

range = 1000:1200;

figure;
plot(range, prep.data{2}.calcium(range), '-k')
for i = 1:length(range)
    if prep.data{2}.spikes(range(i)) > 0
        hold on;
        plot(range(i)*[1 1], 0.1*[-1 1] -1.5, '-k')
    end
end
xlim([range(1) range(end)])
set(gca, 'box', 'off', 'tickdir', 'out')
axis off