 function anaT = analysis_table(lfps, splittype)
%%
% generate a table to see how variables are related to one another
% output dims = sessions x variables
%

% default values
if nargin < 2; splittype = 'drug'; end

% single pair or more? ===========================
switch splittype
    case 'drug'
        pairnames = {'base', 'NaCl'; 'base', '5HT'};
    case 'sc'
        pairnames = {'small sc', 'large sc'; 'small sc', 'large sc'};
end

datast = lfps.LFP_prepro;
lists = [lfps.is5ht, lfps.animal];
lenses = size(lists, 1);

% exclude sessions with no spiking data
outs = zeros(lenses, 1);
for i = 1:lenses
    if datast{i}.cond(1).sta.nspk==0 || datast{i}.cond(2).sta.nspk==0
        outs(i) = 1;
    end
end
datast(outs==1) = [];
lists(outs==1, :) = [];
lenses = size(lists, 1);

% File path information
path = mfilename( 'fullpath' );
if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

% Path finding for external libraries
dir_path = strjoin(parts(1:end-2), '/');

% For saving data table
data_path = [dir_path, '/resources/Data/tables/'];
if ~exist(data_path, 'dir'); mkdir(data_path); end

% stimulus type ============================
stmidx = ones(lenses, 1);
stmdur = 2;                 % Stimulus duration in the experimental design

% preset params
ss = size(stmidx, 2);

% variable names
varnames = {'fr base', 'fr drug', 'r rate', 'sta amp base', 'sta amp drug', 'd sta amp', ...
    'gamma pow base', 'gamma pow drug', 'd gamma pow', 'low-freq pow base', 'low-freq pow drug', 'd low-freq pow', ...
    'coh gamma base', 'coh gamma drug', 'd coh gamma', 'coh low-freq base', 'coh low-freq drug', 'd coh low-freq',...
    'ntr base', 'ntr drug'};
lenv = length(varnames);

% LFP bands
bandnames = {'gamma (30-48)', 'low-freq (3-10)'};
bandrange = {[30, 48], [3, 10]}; 
lenb = length(bandnames);
    
% Initialization
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
        sta_t = linspace(-datast{i}.wnd, datast{i}.wnd, size(datast{i}.cond(1).sta.mean, 2));
        
        % Select data within window
        sta_analysiswnd = sta_t >= -0.05 & sta_t <= 0.05;
       
        % Determine the min and max of stLFP
        maxv = max(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_analysiswnd));
        minv = min(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_analysiswnd));
        at(i, 4, s) = maxv - minv;

        % sta amp drug
        maxv = max(datast{i}.cond(2).sta.mean(stmidx(i, s), sta_analysiswnd));
        minv = min(datast{i}.cond(2).sta.mean(stmidx(i, s), sta_analysiswnd));
        at(i, 5, s) = maxv - minv;

        % d sta amp
        at(i, 6, s) = at(i, 4, s) - at(i, 5, s);
        
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

        for b = 1:lenb
            % frequency range
            frange = bandrange{b}(1) <= freq & bandrange{b}(2) >= freq;
            frangec = bandrange{b}(1) <= freqc & bandrange{b}(2) >= freqc;
            
            % overall power =============
            % base 
            at(i, 7 + 3*(b-1), s) = nanmean(nanmean(10*log10(S1(frange, spe_t > datast{i}.window{end}(1)))));
            % drug 
            at(i, 8 + 3*(b-1), s) = nanmean(nanmean(10*log10(S2(frange, spe_t > datast{i}.window{end}(1)))));
            % d
            at(i, 9 + 3*(b-1), s) = at(i, 7 + 3*(b-1), s) - at(i, 8 + 3*(b-1), s);          
            
            % coherence ===============
            % base 
            at(i, 13 + 3*(b-1), s) = nanmean(coh1(frangec));
            % drug
            at(i, 14 + 3*(b-1), s) = nanmean(coh2(frangec));
            % d (difference)
            at(i, 15 + 3*(b-1), s) = at(i, 13 + 3*(b-1), s) - at(i, 14 + 3*(b-1), s);  
        end   
        
        % the number of trials
        at(i, 19, s) = datast{i}.cond(1).ntr;
        at(i, 20, s) = datast{i}.cond(2).ntr;
    end    
end

% fill nan values 
for i = 1:size(at, 2)
    nans = isnan(at(:, i));
    at(nans, i) = nanmedian(at(:, i));
end

% autosave
anaT.table = at;
anaT.varnames = varnames;
anaT.lists = lists;
anaT.pairnames = pairnames;

% For saving data table
save([data_path '/anaT_' splittype '.mat'], 'anaT')
disp('table saved!')