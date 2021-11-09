function ex = loadCluster(fname, varargin)
% all loading and preprocessing goes here


%% define variables and parse input
rew_flag = 1; % rewarded trials only
ocul = nan; % occularity condition
lfp_flag = true; % combines the spike and LFP file and preprocesses the lfp
filtering = 1;
savepath = nan;
j = 1;

while j<=length(varargin)
    switch varargin{j}
        case 'reward'
            rew_flag = varargin{j+1};
        case 'ocul'
            ocul = varargin{j+1};
        case 'loadlfp'
            lfp_flag = varargin{j+1};
        case 'filtering'
            filtering = varargin{j+1};
        case 'save'
            savepath = varargin{j+1};
    end
    j=j+2;
end

%% load ex file
if mean(ismember('gpfs0', cd))==1
    fname = strrep(fname, '\', '/');
    fname = strrep(fname, 'Z:/', '/gpfs01/nienborg/group/');
end
load(fname, 'ex');
if contains(fname, 'c1')
    cn = 1;
elseif contains(fname, 'c2')
    cn = 2;
elseif contains(fname, 'c3')
    cn = 3;
end

%% try to load the corresponding lfp file and add the signal to each Trial
if lfp_flag
    try
        exwlfp = load(strrep(fname, ['c' num2str(cn)], 'lfp'), 'ex'); % ex with lfp information       
    catch
        sorter_i = strfind(fname, 'sort');
        sorter = fname(sorter_i+4:sorter_i+5);
        if strcmp(sorter, 'HN')
            fname(sorter_i+4:sorter_i+5) = 'LH';
        end
        exwlfp = load(strrep(fname, ['c' num2str(cn)], 'lfp'), 'ex'); % ex with lfp information      
    end
    exwlfp = exwlfp.ex;
        
        for i =1:length(ex.Trials)
            ex.Trials(i).LFP = exwlfp.Trials(i).LFP;
            ex.Trials(i).LFP_ts = exwlfp.Trials(i).LFP_ts;
        end
        ex.LPF = exwlfp.LFP;
%     catch
%         fprintf('no lfp file was found \n%s \n\n', getFname(ex));
%     end
end


%% mark trials at the beginning of fixation
counter = 1; % current trial - modified within the loop
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
    co(co==0)= ex.exp.e1.blank;  
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
if ~contains(fname, 'lfp'); addspkRate;  end

%% reduce to trials with specified ocular condition
if ~isnan(ocul); ex.Trials = ex.Trials([ex.Trials.me] == ocul); end

%% if the combining of LFP and spike ex was successfull, preprocess the LFP
% signal immediately
if isfield(ex.Trials, 'LFP')
%     ex = frequAnalysis(ex, varargin{:}); 
    ex = preprocess(ex, filtering);
end

%%
% autosave
if ~isnan(savepath)
    save(savepath, 'ex')
    disp([savepath ': preprocessed and saved!'])
end