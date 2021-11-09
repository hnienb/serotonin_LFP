function c2s_example
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
data_path = [dir_path, '/resources/Data/c2s/data/ka_0173_c1_sortLH_all.grating.ORxRC_NaCl/'];
addpath(data_path)

addpath(genpath([dir_path, '/external_libraries/']))

cb = cbrewer('div', 'PiYG', 10);
prep = load([data_path '/preprocessed_base_cv10.mat']);
prep = prep.data;
pred = load([data_path '/predictions_base_cv10.mat']);
pred = pred.data;

cols = {[cb(end, :); cb(1,:)], [cb(end-2, :); cb(3,:)]};
thre = 0;
movebin = 200;
x = 1:movebin;
for i = 1:2
    start = 1;
    figure;
    while start + movebin < length(prep{i}.calcium)  
        % trial
        range = start:start+movebin-1;

        % original data
        yorig = double(prep{i}.spikes(range));
        % model prediction
        ypred = pred{i}.predictions(range);
        % correlation
        r = corrcoef(yorig, ypred);
        if thre < r(1,2)
            % update threshold
            thre = r(1,2);

            % plot
            clf
            plot(x, prep{i}.calcium(range), '-', 'color', cols{i}(1,:))
            hold on;
            yorig = yorig/max(yorig);
            plot(x, yorig -2, '-', 'color', cols{i}(1,:))
   
            ypred = ypred/max(ypred);
            plot(x, ypred -3, '-', 'color', cols{i}(end,:))
            set(gca, 'box', 'off', 'tickdir', 'out')
            xlim([x(1) x(end)])
            axis off
        end
        start = start + movebin;
    end
    text(-40, 0, 'LFP', 'fontsize', 6, 'color', cols{i}(1,:))
    text(-40, -2, 'original spike rate', 'fontsize', 6, 'color', cols{i}(1,:))
    text(-40, -3, 'predicted spike rate', 'fontsize', 6, 'color', cols{i}(2,:))
    xlim([-40 x(end)])
end