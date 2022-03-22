function [spk, spkc] = getSpks(trials, wnd, fieldname)
% spikes within the stimulus presentation time
% INPUT: trials ... ex.Trials
%        wnd ... [start, end] of analysis window
%        fieldname ... 'Spikes' or 'oSpikes'
%

% time offset for spikes (350ms after stimulus onset as default)
if nargin < 2; wnd = [0.35, 0]; end
if nargin < 3; fieldname = 'Spikes'; end
    
awnd_strt = wnd(1); %<- 350ms after stimulus onset (default)
awnd_end = wnd(2); %<- 000ms before stimulus end (default)
N = length(trials);
spk = cell(1, N);
spkc = zeros(1, N);
for i = 1:N
    t_strt = trials(i).Start - trials(i).TrialStart;
    t_end = t_strt(end) + mean(diff(t_strt));
    
    spk_tr = trials(i).(fieldname)(trials(i).(fieldname) >= t_strt(1)+awnd_strt & ...
        trials(i).(fieldname) <= t_end-awnd_end) - t_strt(1);
    spk{i} = round(spk_tr*1000)/1000;
    spkc(i) = length(spk_tr);
end
