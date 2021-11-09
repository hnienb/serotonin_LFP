 function anaT = analysis_table(lfps, splittype)
%%
% generate a table to see how variables are related to one another
% sessions x variables x stimulus type
%

if nargin < 2; splittype = 'drug'; end

% single pair or more? ===========================
switch splittype
    case 'drug'
        pairnames = {'base', 'NaCl'; 'base', '5HT'};
    case 'sc'
        pairnames = {'small sc', 'large sc'; 'small sc', 'large sc'};
    case 'pupil'
        pairnames = {'small ps', 'large ps'; 'small ps', 'large ps'};
end
if isfield(lfps, 'lfplist')
    datast = lfps.LFP_prepro;
    lists = [lfps.goodunit', lfps.is5ht', lfps.animal'];
else
    datast{1} = lfps;
    lists = [1, lfps.is5ht, lfps.ismango];
end

disp(['5HT: ' num2str(sum(lists(:,2)==1)) ' pairs'])
disp(['NaCl: ' num2str(sum(lists(:,2)==0)) ' pairs'])
lenses = size(lists, 1);

% exclude no data sessions
if strcmp(datast{1}.stm.param, 'rc')
    outs = zeros(lenses, 1);

    for i = 1:lenses
        if datast{i}.cond(1).sta.nspk==0 || datast{i}.cond(2).sta.nspk==0
            outs(i) = 1;
        end
    end
    datast(outs==1) = [];
    lists(outs==1, :) = [];
    lenses = size(lists, 1);
end


% File path information
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-3), '/');

addpath(genpath([dir_path '/external_libraries/CircStat2012a']))
addpath(genpath([dir_path '/external_libraries/chronux_2_12/chronux_2_12']))

% For saving data table
data_path = [dir_path, '/resources/Data/tables/'];
if ~exist(data_path, 'dir'); mkdir(data_path); end


% stimulus type ============================
stmdur = 0.45; % Stimulus duration in the experimental design

% Removed other cases since only working with reverse correlation data
switch datast{1}.stm.param
    case 'rc' % Reverse correlation
        stmidx = ones(lenses, 1);
        stmdur = 2;
end

% for oher structs
lenses = size(lists, 1);

% remove no-data sessions
datast(stmidx(:,1)==0) = [];
lists(stmidx(:,1)==0, :) = [];
stmidx(stmidx(:,1)==0, :) = [];

% preset params
ss = size(stmidx, 2);

% variable names
varnames = {'fr base', 'fr drug', 'r rate', ...
    'sta amp base', 'sta amp drug', 'd sta amp', 'sta Tmin base', 'sta Tmin drug', 'd sta Tmin', 'theta pow base', 'theta pow drug', ...
    'd theta pow', 'alpha pow base', 'alpha pow drug', 'd alpha pow', 'beta pow base', 'beta pow drug', ...
    'd beta pow', 'gamma pow base', 'gamma pow drug', 'd gamma pow', 'low-freq pow base', 'low-freq pow drug', 'd low-freq pow', ...
    'broadband pow base', 'broadband pow drug', 'd broadband pow', 'theta res base', 'theta res drug', 'd theta res', ...
    'alpha res base', 'alpha res drug', 'd alpha res', 'beta res base', 'beta res drug', ...
    'd beta res', 'gamma res base', 'gamma res drug', 'd gamma res', 'low-freq res base', 'low-freq res drug', 'd low-freq res', ...
    'broadband res base', 'broadband res drug', 'd broadband res', 'coh theta base', 'coh theta drug', 'd coh theta', ...
    'coh alpha base', 'coh alpha drug', 'd coh alpha', 'coh beta base', 'coh beta drug', 'd coh beta', ...
    'coh gamma base', 'coh gamma drug', 'd coh gamma','coh low-freq base', 'coh low-freq drug', 'd coh low-freq',...
    'coh broadband base', 'coh broadband drug', 'd coh broadband', 'phase theta base', 'phase theta drug', 'd phase theta', ...
    'phase alpha base', 'phase alpha drug', 'd phase alpha', 'phase beta base', 'phase beta drug', 'd phase beta', ...
    'phase gamma base', 'phase gamma drug', 'd phase gamma','phase low-freq base', 'phase low-freq drug', 'd phase low-freq',...
    'phase broadband base', 'phase broadband drug', 'd phase broadband', ...
    'lfp signal base', 'lfp signal drug', 'd lfp signal', 'lfppow ratio base', 'lfppow ratio drug', 'd lfppow ratio', ...
    'low-freq ps weight', 'low-freq dps weight', 'low-freq drug weight', 'low-freq drugxps weight', 'low-freq drugxdps weight',...
    'ntr base', 'ntr drug'};

lenv = length(varnames);

mdly = 8;
mdlout = {[6:11]};
lenm = length(mdly);
cv = 2;
lam = 0.084;

% TODO FIX THE BANDNAMES TO MATCH THE ACTUAL VALUES - CHECK MANUSCRIPT

bandnames = {'theta (3-7)', 'alpha (8-13)', 'beta (14-24)', 'gamma (36-48)', ...
    'low-freq (3-10)', 'broadband (3-48)'};
% TODO: Why are the actual band values different bandnames?
bandrange = {[3, 7], [8, 12], [15, 25], [40, 48], [3, 10], [3, 48]};
lenb = length(bandnames);
    

at = nan(lenses, lenv, ss);

for i = 1:lenses
    for s = 1:ss
        % firing rate  =====================
        % fr base
        at(i, 1, s) = datast{i}.cond(1).spk_tu{1}(1);
        % fr drug
        at(i, 2, s) = datast{i}.cond(2).spk_tu{1}(1);
        % r rate
        at(i, 3, s) = datast{i}.cond(2).spk_tu{1}(s, 1)/datast{i}.cond(1).spk_tu{1}(s, 1);
        
        % stLFP ================================
        % sta amp base
        sta_t = linspace(-datast{i}.wnd, datast{i}.wnd, ...
            size(datast{i}.cond(1).sta.mean, 2));

        sta_analysiswnd = sta_t >= -0.05 & sta_t <= 0.05;
        sta_analysis_time = sta_t(sta_analysiswnd);
        [maxv, maxt] = max(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_analysiswnd));
        [minv, zerotb] = min(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_analysiswnd));

        zerotb = sta_analysis_time(zerotb);

        at(i, 4, s) = maxv - minv;
        % sta amp drug
        [maxv, maxt] = max(datast{i}.cond(2).sta.mean(stmidx(i, s), sta_analysiswnd));
        [minv, zerotd] = min(datast{i}.cond(2).sta.mean(stmidx(i, s), sta_analysiswnd));

        zerotd = sta_analysis_time(zerotd);

        at(i, 5, s) = maxv - minv;
        % d sta amp
        at(i, 6, s) = at(i, 4, s) - at(i, 5, s);
        % sta Tmin base
        at(i, 7, s) = zerotb;
        % sta Tmin drug
        at(i, 8, s) = zerotd;
        % d sta Tmin
        at(i, 9, s) = at(i, 7, s) - at(i, 8, s);
        
        % spectral power & response & coherence & phase & lfp SD & LFP vs PS or FR ======
        freq = datast{i}.cond(1).spectrogram.f{1};
        spe_t = datast{i}.cond(1).spectrogram.t{1};
        lent = length(spe_t);
        spe_t = linspace(-0.1, stmdur, lent);
        S1 = datast{i}.cond(1).spectrogram.S{stmidx(i, s)};
        S2 = datast{i}.cond(2).spectrogram.S{stmidx(i, s)};
        sz = size(datast{i}.cond(1).spectrogram.S{stmidx(i, s)});
        if sz(2) < sz(1)
            S1 = datast{i}.cond(1).spectrogram.S{stmidx(i, s)}';
            S2 = datast{i}.cond(2).spectrogram.S{stmidx(i, s)}';
        end
        freqc = datast{i}.cond(1).coherence.f{stmidx(i, s)};
        coh1 = datast{i}.cond(1).coherence.C{stmidx(i, s)};
        coh2 = datast{i}.cond(2).coherence.C{stmidx(i, s)};
        phase1 = datast{i}.cond(1).coherence.phi{stmidx(i, s)};
        phase2 = datast{i}.cond(2).coherence.phi{stmidx(i, s)};
        for b = 1:lenb
            % frequency range
            frange = bandrange{b}(1) <= freq & bandrange{b}(2) >= freq;
            frangec = bandrange{b}(1) <= freqc & bandrange{b}(2) >= freqc;
            
            % overall power =============
            % base 
            at(i, 10 + 3*(b-1), s) = nanmean(nanmean(10*log10(S1(frange, spe_t > datast{i}.window{end}(1)))));
            % drug 
            at(i, 11 + 3*(b-1), s) = nanmean(nanmean(10*log10(S2(frange, spe_t > datast{i}.window{end}(1)))));
            % d
            at(i, 12 + 3*(b-1), s) = at(i, 10 + 3*(b-1), s) - at(i, 11 + 3*(b-1), s);
            
            % response power ===========
            % base 
            at(i, 28 + 3*(b-1), s) = at(i, 10 + 3*(b-1), s) - ...
                nanmean(nanmean(S1(frange, spe_t <= 0)));
            % drug 
            at(i, 29 + 3*(b-1), s) = at(i, 11 + 3*(b-1), s) - ...
                nanmean(nanmean(S2(frange, spe_t <= 0)));
            % d
            at(i, 30 + 3*(b-1), s) = at(i, 28 + 3*(b-1), s) - at(i, 29 + 3*(b-1), s);            
            
            % coherence ===============
            % base 
            at(i, 46 + 3*(b-1), s) = nanmean(coh1(frangec));
            % drug
            at(i, 47 + 3*(b-1), s) = nanmean(coh2(frangec));
            % d (difference)
            at(i, 48 + 3*(b-1), s) = at(i, 46 + 3*(b-1), s) - at(i, 47 + 3*(b-1), s);     
            
            % phase ==================
            % base 
            at(i, 64 + 3*(b-1), s) = maxcoh_rad(coh1(frangec), phase1(frangec));
            % drug
            at(i, 65 + 3*(b-1), s) = maxcoh_rad(coh2(frangec), phase2(frangec));
            % d
            at(i, 66 + 3*(b-1), s) = at(i, 64 + 3*(b-1), s) - at(i, 65 + 3*(b-1), s);     
        end  
        % LFP res base 
        at(i, 82, s) = nanstd(datast{i}.cond(1).mat{s}(:, 6));
        % LFP res drug
        at(i, 83, s) = nanstd(datast{i}.cond(2).mat{s}(:, 6));
        % d LFP res
        at(i, 84, s) = at(i, 82, s) - at(i, 83, s);  
        % LFP power ratio base 
        at(i, 85, s) = at(i, 10 + 3*(5-1), s)/at(i, 10 + 3*(4-1), s);
        % LFP power ratio drug
        at(i, 86, s) = at(i, 11 + 3*(5-1), s)/at(i, 11 + 3*(4-1), s);
        % d LFP power ratio
        at(i, 87, s) = at(i, 85, s) - at(i, 86, s);  
        
        % lasso GLM weights =====================
        % x and Y
        mat0 = datast{i}.cond(1).mat{stmidx(i, s)};
        mat2 = datast{i}.cond(2).mat{stmidx(i, s)};
        ntr0 = size(mat0, 1);
        ntr2 = size(mat2, 1);
        X = [[mat0(:, 3); mat2(:, 3)], ... % pupil size
            [mat0(:, 4); mat2(:, 4)], ... % pupil size derivative
            [zeros(ntr0, 1); ones(ntr2, 1)], ... % base or drug
            [zeros(ntr0, 1); ones(ntr2, 1)].*[mat0(: , 3); mat2(:, 3)], ... % interaction
            [zeros(ntr0, 1); ones(ntr2, 1)].*[mat0(: , 4); mat2(:, 4)], ... % interaction
            [mat0(:, 5); mat2(:, 5)], ... % LFP res
            [mat0(:, 6); mat2(:, 6)], ... % ddLFP
            [mat0(:, 7); mat2(:, 7)], ... % low-freq
            [mat0(:, 8); mat2(:, 8)], ... % gamma
            [mat0(:, 1); mat2(:,1)]/stmdur, ... % spike counts --> firing rate (su)
            [mat0(:, 2); mat2(:,2)]/stmdur, ... % spike counts --> firing rate (mu)
            ]; 
        
        % k-fold cross-validation
        cvidx1 = crossvalind('Kfold', ntr0, cv);
        cvidx2 = crossvalind('Kfold', ntr2, cv);
        cvidx = [cvidx1; cvidx2];
        
        % full model fitting
        for m = 1:lenm
            y = X(:, mdly(m));
            predictors = zscore(X(:, ~ismember(1:size(X, 2), mdlout{m})));
            
            % model prediction 
            beta = nan(cv, size(predictors, 2));
            for k = 1:cv
                beta(k, :) = lassoglm(predictors(cvidx==k, :), y(cvidx==k), 'normal', 'lambda', lam);
            end
            beta = mean(beta, 1);
            for k = 1:length(beta)
                at(i, 88 + (k-1), s) = beta(k); 
            end
        end
        
        % the number of trials
        at(i, 93, s) = ntr0;
        at(i, 94, s) = ntr2;
    end    
end
% fillnan
for i = 1:size(at, 2)
    nans = isnan(at(:, i));
    at(nans, i) = nanmedian(at(:, i));
end

% model performance (decoding)
try    
    load([dir_path '/resources/Data/c2s/met_cv10.mat'])
    if strcmp(splittype, 'drug')
        idx = [3, 4];
    else
        idx = [7, 8];
    end
    at = [at, met(:, idx), met(:, idx(1))-met(:,idx(2))];
    varnames{lenv+1} = 'mprf base';
    varnames{lenv+2} = 'mprf drug';
    varnames{lenv+3} = 'd mprf';
catch
    disp('No model performance added')
end

% autosave
anaT.table = at;
anaT.varnames = varnames;
anaT.lists = lists;
anaT.pairnames = pairnames;

% For saving data table
save([data_path '/anaT_' splittype '.mat'], 'anaT')
disp('table saved!')


function rad = maxcoh_rad(coh, rad)
lenv = length(coh);
if lenv < 4
    rad = circ_mean(rad);
else
    [~, maxi] = max(coh);
    if maxi==1
        idx = [1 2 3];
    elseif maxi==lenv
        idx = [lenv-2 lenv-1 lenv];
    else
        idx = [maxi-1 maxi maxi+1];
    end
    rad = circ_mean(rad(idx));
end