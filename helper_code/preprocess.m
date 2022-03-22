function ex = preprocess(ex, filt)
%%
% preprocess LFP and eye data
% INPUT: ex ... ex-file after 'loadCluster.m'
%        filt ... 0, no filter; 1, filtering
% OUTPUT: ex-file with preprocessed LFP and eye data
%

if nargin < 2; filt = 1; end

%% define variables
stimdur = getStimDur(ex); % stimulus presentation duration

% stimulus per trials
label_seq = label4StmPerTr(ex);

% sampling frequency
Fs = 1000;              

% time relative to stimulus onset to filter and interpolate LFP
t_off = -0.1;

% loop to get data ==================================
ex.Trials = ex.Trials(label_seq > 0);
label_seq = label_seq(label_seq > 0);
N = length(ex.Trials);

% remove trials with no lfp or eye data
oks = ones(1, N);
for n = 1:N
    if isempty(ex.Trials(n).LFP)
        oks(n) = 0;
    end
end
ex.Trials = ex.Trials(oks==1);
N = length(ex.Trials);

% initialization
lfps = [];
lfps_timing = ones(N, 2);

for n = 1:N    
    % label seq ==========
    ex.Trials(n).label_seq = label_seq(n);
    
    % LFP ================================
    % timing information for LFP
    t_frame = ex.Trials(n).Start - ex.Trials(n).TrialStart; % time of frame onsets
    t_lfp = ex.Trials(n).LFP_ts - t_frame(1) ; % time rel:stimulus onset
    time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
    ex.Trials(n).LFP_prepro_time = time - t_frame(1);
        
    % reduce the lfp signal to the period of stimulus presentation
    lfps_temp = interp1(t_lfp, ex.Trials(n).LFP, ex.Trials(n).LFP_prepro_time);
    lfps = [lfps, lfps_temp];
    
    % start and end of the trial
    if n > 1
        lfps_timing(n, 1) = lfps_timing(n-1, 2) + 1;
    end
    lfps_timing(n, 2) = lfps_timing(n, 1) + length(lfps_temp) - 1;
end

% remove nans ===============================
lfps = nan_interp(lfps);

% filtering ==================================
if filt
    % filter parameters
    notchf = [49 51];     % notch filter frequency
    notchord = 2;         % filter order

    bpf = [3 90];         % bandpass filter cutoff frequency
    bpord = [8 10];       % filter order
    
    % define notch filter
    d = designfilt('bandstopiir', 'FilterOrder', notchord, ...
                   'HalfPowerFrequency1', notchf(1), 'HalfPowerFrequency2', ...
                   notchf(2), 'DesignMethod','butter','SampleRate', Fs);
    lfps = filtfilt(d, lfps);

    % define bandpass filter for LFPs (Nauhaus et al., 2009)
    [b_high, a_high] = butter(bpord(1), bpf(1)/(Fs/2), 'high');
    [b_low, a_low] = butter(bpord(2), bpf(2)/(Fs/2), 'low');

    % LFP =======================
    lfps = filtfilt(b_high, a_high, lfps);
    lfps = filtfilt(b_low, a_low, lfps);
end

for n = 1:N
    ex.Trials(n).LFP_prepro = lfps( lfps_timing(n,1):lfps_timing(n,2) );
end

% remove unnecessary fields
rmf = 'LFP';
ex.Trials = rmfield(ex.Trials, rmf);

% Make all the LFP data the same length (sanity-check)
lens = length(ex.Trials(end).LFP_prepro_time);
for n = 1:N
    ex.Trials(n).LFP_prepro = ex.Trials(n).LFP_prepro(1:lens);
end


function avgstimdur = getStimDur(ex)
% returns the averaged and rounded stimulus presentation duration across
% trials
try
    t_frame = cellfun(@(x, y) x-y, {ex.Trials.Start}, {ex.Trials.TrialStart}, ...
        'UniformOutput', 0); % time of frame onsets
    stimdur = cellfun(@(x) x(end)+mean(diff(x)) - x(1), t_frame);

    % also round the average stimulus duration to 2 digits precision
    avgstimdur = round(mean(stimdur), 2);
catch
    if ex.exp.StimPerTrial == 4
        avgstimdur = 0.45;
    elseif ex.exp.StimPerTrial == 1
        avgstimdur = 2;
    end
end
