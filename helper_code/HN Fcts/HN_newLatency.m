function res = HN_newLatency(res)
%hn:  modified latency computation I used for SFN poster (Hruba et al. 2015)
% 11/11/15

vars = sqrt(res.vars);

if max((vars))>mean((vars(1:25)))*5
    vars = vars-mean(vars(1:25));
    idx = find((vars)>=max((vars))/2);
    lat = idx(1);
    res.latencyToHalfMax = res.delays(lat);
    res.dur = res.delays(idx(end)) - res.latencyToHalfMax ; % CL added 11.11.2015
end

% added by CL 19.02.2016

if ~isempty(res.sdfs.extras)
    res.vars2 = var([horzcat(res.sdfs.s{:}), res.sdfs.extras{1}.sdf], 0, 2);
else
    res.vars2 = var(horzcat(res.sdfs.s{:}), 0, 2); 
end