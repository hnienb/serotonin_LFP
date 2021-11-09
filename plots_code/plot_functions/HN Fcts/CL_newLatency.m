function res = CL_newLatency(res, p_flag, ml_lat_flag)
% res = CL_newLatency(res, p_flag, ml_lat_flag)
% 
% 
% This function computes two latency metrics 
% res.latFP - the parametric ML latency introduced by Friedmann and Priebe
%               (see friedmanpriebe.mat),
% res.lat   - and the non-parametric latency defined as the time at
%               half-max repsonse deviation relative to the baseline
%               variation after frame onset
% 
% In both cases it uses the deviation across SDFs as a signature of the
% response dynamic. Addtionally,
% res.dur   - is the duration of the response and computed as the time
%               between the response onset, i.e. the non-parametric latency
%               and the second time, the response deviation crosses the
%               half-max threshold.
% 
% 
% In case where there are two stimulus dimensions, the latency is estimated
% for the maxmimum value of the second stimulus dimension, assuming that
% this is contrast.
% 
% If ml_lat_flag is false, the  ML estimate is skipped
% If ml_lat_flag and p_flag is true, the ML estimate is plotted.
% 
% 
%@CL 22.01.2016
%@CL adapted 31.05.2017     added additional input arguments


% estimate the latency using the deviation across stimulus selective
% response components.% If there is a second stimulus dimension, assume it
% is contrast and restrict the analysis to the SDFs for orientation with
% maximum contrast.
[~,co_idx]= max(res.sdfs.y(1,:));

% add the blank response to the orientation selective components
if ~isempty(res.sdfs.extras)
    res.vars2(:) = var([horzcat(res.sdfs.s{:, co_idx}), res.sdfs.extras{1}.sdf], 0, 2);
else
    res.vars2(:) = var(horzcat(res.sdfs.s{:, co_idx}), 0, 2);
end

% call the latency computations
[res.sdfs.lat2hmax, res.sdfs.dur, res.sdfs.latFP] = ...
    CL_newLatency_helper(res.vars2(:), res.times, p_flag, ml_lat_flag);

res.lat = res.sdfs.lat2hmax;
res.dur = res.sdfs.dur;
res.latFP = res.sdfs.latFP;

end




function [lat2hmax, dur, latfp, pPoisson] = CL_newLatency_helper(vars, times, p_flag, ml_lat_flag)



% initialize variables
lat2hmax    = -1; % non-parametric responte latency estimate
dur         = -1; % response duration
latfp       = -1; % maximum likelihood estiamte of the response latency 
pPoisson    = 0;  % p-value for the ML estimate  

sd      = sqrt(vars);  % deviation acriss SDFs
noise   = mean(sd(200:400));  % baseline variability


% the response peak must be at least 4 times the baseline noise to qualify
% as response signature. this criteria is arbitrary.
if max(sd)>mean(noise)*4
    
    sd2 = sd-noise;     % normalization for baseline variability
    
    % first time the response deviation crosses the half-max threshold
    idx = find( sd2 >= (max(sd2)/2), 1, 'first');  
    lat2hmax = times(idx)/10;
    
    % last time the response deviation crosses the half-max threshold
    idx = find( sd2 >= (max(sd2)/2), 1, 'last');  
    dur = times(idx)/10 - lat2hmax;
    
    try
        if ml_lat_flag
            % some additional arguments improve the ML
            [latfp, ~, pPoisson] = friedmanpriebe(round(sd(200:end).*100), ...
                'minTheta', 250, 'responseSign', 0, 'graphics', p_flag);
        else
            latfp = -1;
        end
    catch
       disp('there were problems computing the latency using the ML algorithm') 
    end
    
    % convert to ms
    latfp = latfp/10;
    
end

end

