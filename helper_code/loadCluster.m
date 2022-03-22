function ex = loadCluster(fname, varargin)
% all loading and preprocessing goes here

%% define variables and parse input
rew_flag = 1;    % rewarded trials only
filtering = 1;
savepath = nan;
j = 1;

while j<=length(varargin)
    switch varargin{j}
        case 'reward'
            rew_flag = varargin{j+1};
        case 'filtering'
            filtering = varargin{j+1};
        case 'save'
            savepath = varargin{j+1};
    end
    j=j+2;
end

%% load ex file
load(fname, 'ex');

%% mark trials at the beginning of fixation
counter = 1;            % current trial - modified within the loop
for i = 1:length(ex.Trials)
    ex.Trials(i).FixCounter = counter; % number of stimulus during fixation
    if ex.Trials(i).RewardSize > 0 || counter == 4 || ex.Trials(i).st == 0
        counter = 0; % reset counter for new fixation period
    end     
    counter = counter+1;
end

%% add blank
addBlank;

%% adapt the stimulus information
if strcmp(ex.exp.e1.type, 'or')
    % collapse orientations to [0 180], ignore blanks for this
    trials = mod([ex.Trials.or], 180)'; trials = num2cell(trials);
    idx = [ex.Trials.or]<=360; 
    [ex.Trials(idx).or] = deal(trials{idx});
elseif strcmp(ex.exp.e1.type, 'co')
    % take contrast values and replace 0 contrast with blank value
    co = [ex.Trials.co]'; 
    co(co==0) = ex.exp.e1.blank;  
    co = num2cell(co);
    [ex.Trials.co] = deal(co{:});
end

% when analysing flashed grating experiments, the varying entries in
% Trials.or are misleading.
if contains(fname, 'RC')
    [ex.Trials(:).or] = deal(1);
end

%% exclude trials with fixation break
if rew_flag; ex.Trials = ex.Trials([ex.Trials.Reward] == 1); end

%% compute and add the spike count and resutling rate during stimulus presentation
if ~contains(fname, 'lfp'); addspkRate; end

%% if the combining of LFP and spike ex was successful, preprocess the LFP signal
if isfield(ex.Trials, 'LFP')
    ex = preprocess(ex, filtering);
end

%% autosave
if ~isnan(savepath)
    save(savepath, 'ex')
    disp([savepath ': preprocessed and saved!'])
end