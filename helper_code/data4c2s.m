function [para, data] = data4c2s(ex0, ex2, lfpfield, z)
%%
% format LFP data from the RC experiments such that they can be used in
% 'c2s' analysis
%
% source code for c2s : 'https://github.com/lucastheis/c2s'
%
% INPUT: ex0, ex2 ... ex file after 'loadCluster.m'
%        lfpfield ... 'LFP_prepro', 'iLFP_prepro'
%        z ... z-score (0 or 1)
%
% OUTPUT: para ... struct; session information
%         data ... cell array; something as follows to match the convention of 'c2s'
% data{1} =
%      calcium: [1x71985 double] .... in our case, preprocessed LFP
%       spikes: [1x71985 uint16] .... spikes
%          fps: 99.9998          .... sampling rate 
%     cell_num: 1                .... 1, baseline; 2, drug
%

% inputs =================
if nargin < 3; lfpfield = 'LFP_prepro'; end
if nargin < 4; z = 1; end

% stimulus type ==============
[stimparam, vals] = getStimParam(ex0);
lenv = length(vals);
if lenv==1 && strcmp(stimparam, 'or')
    para.stm.param = 'rc';
else
    para.stm.param = stimparam;
end
para.stm.vals = vals;

% drug & animal ==================
try
    para.ismango = contains(ex2.fileName, 'ma');
    para.is5ht = contains(ex2.fileName, '5HT');
catch
    para.ismango = nan;
    para.is5ht = nan;
end
    
% preset parameters ==================
fs = 1000;
para.ts = ex0.Trials(end).LFP_prepro_time; 
para.window = {[0 2]};        
stmdur = 2;
para.wnd = 0.07;


% LFP & spikes as a function of stimulus
exs = {ex0, ex2};
for d = 1:2
    % initialization
    para.cond(d).lfpfull = cell(1, lenv);
    para.cond(d).spk = cell(3, lenv);
    para.cond(d).spkc = cell(3, lenv);
    para.cond(d).ntr = zeros(1, lenv);
    para.cond(d).trials = cell(1, lenv);
    para.cond(d).spk_tu = cell(1, 3); 
    
    % stimulus types
    for i = 1:lenv    
        % trials with the unique stimulus
        para.cond(d).trials{i} = exs{d}.Trials([exs{d}.Trials.(stimparam)] == vals(i));
        para.cond(d).ntr(i) = length(para.cond(d).trials{i});        

        % firing rate per trial in analysis window
        [para.cond(d).spk{1, i}, para.cond(d).spkc{1, i}]  = ...
            getSpks(para.cond(d).trials{i}, [para.window{end}(1), ...
            stmdur - para.window{end}(2)]);    

        % firing rate per trial during stimulus presentation
        [para.cond(d).spk{2, i}, para.cond(d).spkc{2, i}]  = getSpks(para.cond(d).trials{i}, [0 0]);

        % elicited spikes (mean, SD, n)
        dur = para.window{end}(2) - para.window{end}(1);
        para.cond(d).spk_tu{1}(i, :) = [mean(para.cond(d).spkc{1, i})/dur, ...
            std(para.cond(d).spkc{1, i})/dur, sum(para.cond(d).spkc{1, i})];        
        para.cond(d).spk_tu{2}(i, :) = [mean(para.cond(d).spkc{2, i})/stmdur, ...
            std(para.cond(d).spkc{2, i})/stmdur, sum(para.cond(d).spkc{2, i})];    
        
        % preprocessed LFP
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).lfpfull{i}(n, 1:length(para.cond(d).trials{i}(n).(lfpfield))) = ...
                para.cond(d).trials{i}(n).(lfpfield);
        end
    end
end

% zscoring
if z==1
    lfpall = []; 
    for d = 1:2
        for i = 1:lenv    
            for n = 1:para.cond(d).ntr(i)
                lfpall = [lfpall, para.cond(d).trials{i}(n).(lfpfield)];
            end
        end
    end
    me = nanmean(lfpall);
    sd = nanstd(lfpall);
else
    me = zeros(1, 1); sd = ones(1, 1);
end
for d = 1:2
    lenlfp = length(para.cond(d).trials{end}(end).(lfpfield));
    for i = 1:lenv    
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).trials{i}(n).(lfpfield) = (...
                para.cond(d).trials{i}(n).(lfpfield) - me(1))/sd(1);
            para.cond(d).lfpfull{i}(n, 1:lenlfp) = para.cond(d).trials{i}(n).(lfpfield);
        end
    end
end
    
% generate 'data' for c2s analysis
data = cell(1, 2);
c = 1;
j = 1;      % Non-thinned data within analysis window
for d = 1:2 % baseline or drug

    for i = 1:lenv 
        for n = 1:para.cond(d).ntr(i)
            % sampling rate
            data{c}.fps = fs;

            % baseline or drug
            data{c}.cell_num = d;

            % spike times (ms)
            data{c}.spike_times = 1000*para.cond(d).spk{j, i}{n}';

            % LFP
            lfp = para.cond(d).lfpfull{i}(n, para.ts >= para.window{1}(1)...
                & para.ts <= para.window{1}(2));
            data{c}.calcium = lfp;
            c = c + 1;
        end
    end
end

% remove irrelevant fields
rmf = {'ts', 'window', 'cond'};
for r = 1:length(rmf)
    para = rmfield(para, rmf{r});
end

function spk = thinning(spk, fr_ratio)
%%
% perform thinning of spikes
% INPUT: spk ... spike times (1D)
%              fr_ratio ... ratio of mean firing rate across two conditions
%                              (u_small/u_large)
%              wnd ... window to apply thinning
%

if isempty(spk)
    return
end

r = rand(1, length(spk));
spk(r <= 1 - fr_ratio) = []; 