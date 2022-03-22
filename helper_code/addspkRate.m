
%
% add...
% - ...number of spikes
% - ...number of presented frames
% - ...spike rate
%
% @author Corinna Lorenz
% created 22.09.2015
% changed 22.09.2015
% adapted 11.04.2017 - added the option to compute the spikes/cycle/seconds
%


for n = 1:length(ex.Trials)
    
    % only consider correctly executed trials
    if ex.Trials(n).Reward == 1
        t_spks = ex.Trials(n).Spikes;   % time of spike occurance
        f_strt = ex.Trials(n).Start;    % relative time of frame onset
        t_strt = f_strt - ex.Trials(n).TrialStart; % aligned time of stimulus frames onset
        
        frame_dur = mean(diff(t_strt)); % average frame duration
        t_end = t_strt(end)+frame_dur;  % time of stimulus ending
        
        ex.Trials(n).stimDuration   = t_end - t_strt(1); % stimulus duration
        ex.Trials(n).spkCount       = length(find( t_spks>=t_strt(1) & t_spks<=t_end)); % spikes occuring during stimulus presentation
        ex.Trials(n).spkRate        = ex.Trials(n).spkCount / (t_strt(end) - t_strt(1));  % resulting spike rate 
        
    else
        
        ex.Trials(n).stimDuration   = -1;
        ex.Trials(n).spkCount       = -1;
        ex.Trials(n).spkRate        = -1;
        
    end
 
end
% for debugging
% fprintf('tf %1.0f ==> number of full cycles = %1.0f \n', tf, n_cycles);

clearvars t_spks t_strt t_end n frame_dur