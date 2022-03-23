function extract_c2s
%%
% extract info from c2s analysis
%

% path =======================================
path = mfilename( 'fullpath' );

if ispc    % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

rootpath = strjoin(parts(1:end-2), '/');
datapath = [rootpath, '/resources/Data/c2s/data/'];
savepath = [rootpath, '/resources/Data/c2s/'];
dirs = dir(datapath);
dirs(1:2) = [];
lend = length(dirs);

% animal, drug, corr (base), corr (drug), MI (base), MI (drug), corr (hFR), corr (lFR), MI (hFR), MI (lFR)
met = nan(lend, 10); 
mfnames = {'correlation', 'MI'};
fieldnames = {'correlations', 'info'};
lenm = length(mfnames);
for i = 1:length(dirs)
    disp(['ses ' num2str(i) ': ' dirs(i).name])
    if contains(dirs(i).name, 'ka_')
        met(i, 1) = 0;
    elseif contains(dirs(i).name, 'ma_')
        met(i, 1) = 1;
    end
    if contains(dirs(i).name, '5HT')
        met(i, 2) = 1;
    elseif contains(dirs(i).name, 'NaCl')
        met(i, 2) = 0;
    end
    
    % performance
    for m = 1:lenm
        % baseline
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_base_cv10.mat']);
        met(i, 3+2*(m-1)) = model_perf(data, 1, fieldnames{m});
        
        % drug
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_drug_cv10.mat']);
        met(i, 4+2*(m-1)) = model_perf(data, 1, fieldnames{m});
        
        % FR control
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_lowFR_cv10.mat']);
        met(i, 7+2*(m-1)) = model_perf(data, 1, fieldnames{m});
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_highFR_cv10.mat']);
        met(i, 8+2*(m-1)) = model_perf(data, 1, fieldnames{m});
    end
end

save([savepath, '/met_cv10.mat'], 'met', '-v7')

function mpf = model_perf(data, row, name)
if strcmp(name, 'info')
    mpf = nanmean(data.info(row, :)./data.entropy(row, :), 2);
else
    mpf = nanmean(data.(name)(row, :), 2);
end