function stlfp = getSTA(lfptrace, lfptime, spktime, wnd, fs)
%%
% generic function to get spike-triggered averaging LFP
% note that the input 'lfptrace' must be a 0-mean vector
%

if isempty(spktime)
    stlfp = nan;
    return
end

% 0-mean vector
lfptrace = lfptrace - mean(lfptrace(lfptime >= spktime(1)-wnd & lfptime <= spktime(end)+wnd));

% initialization
nspk = length(spktime);
ncol = length(-wnd:(1/fs):wnd);
stlfp = nan(nspk, ncol);
lent = length(lfptrace);

% STA
for i = 1:nspk
    % time range
    tspk = find(lfptime <= spktime(i), 1, 'last');
    tstrt = tspk - wnd*fs;
    tend = tstrt + ncol - 1;
    
    % take all 
    if tstrt < 1
        continue
    elseif tend > lent
        continue
    else
        seg = lfptrace(tstrt:tend);
    end
    
    if ~isempty(seg)
        stlfp(i, :) = seg;
    else
        continue
    end
end

stlfp(any(isnan(stlfp), 2), :) = [];

