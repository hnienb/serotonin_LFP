function res = HN_computeLatencyAndNetSpk_PhaseSelectivity(res,ex, varargin)

% takes the output from PlotRevCorAny_tue, 'res' and ex file to compute
% net-spikes and latency of the orientation selective response.
% If res=[], this function calls PlotRevCorAny_tue.
%
% 10/12/15  hn:   wrote it


lat_flag = 0;
k = 1;
while k <length(varargin)
    switch varargin{k}
        case 'lat_flag'
          lat_flag = varargin{k+1};
    end
    k = k+2;
end


if isempty(res)
    res = HN_PlotRevCorAny_tue(ex,'times',[-200:1600],'sdfw',40, 'noplot', 'exp2', 'phase_seq');
end
   
if isfield(res, 'vars')
    if lat_flag
%     res =  HN_newLatency(res);  
     res =  CL_newLatency(res);  % changed @CL 22.01.2016 
    end
  
        
% set to new function to compute 2D stimuli
if size(res.sdfs.n, 1) > 1 && size(res.sdfs.n, 2) > 1
    res = getNetSpks2D(res, ex);
else
    res = getNetSpks1D(res, ex);
end
end


end


function res = getNetSpks2D(res, ex)

% for the 2D matrix, blanks are left added 'manually', i.e. not added to
% the vector

nspk_extra = 0;
n_extra = 0;
mn_resp_extra = 0;

idx = find([res.times]>=res.delays(1) & [res.times]<=res.delays(end));
for i = 1: length(res.sdfs.extras)
    nspk_extra(i) = res.sdfs.extras{i}.nspk;
    n_extra(i) = res.sdfs.extras{i}.n;
    mn_resp_extra(i) = mean(res.sdfs.extras{i}.sdf(idx)) * n_extra;
end


% compute net spikes-------------------------------------------------------
y = res.sdfs.nspk./res.sdfs.n;
N = res.sdfs.n;

% mean # of spikes per frame across all frames
% note: nspk gives us all spikes that are used in the analysis, i.e. spikes
% are counted several times because of overlapping analysis windows
mn = (sum(sum(res.sdfs.nspk)) + sum(nspk_extra))/...
    (sum(sum(N)) + sum(n_extra));
y1 = y-mn;  % this now gives me the deviation around the mean response
y1_extra = nspk_extra./n_extra -mn;

% mean number of spikes/frame via the SDF
for i = 1:size(res.y,1)
    for j = 1:size(res.y,2)
        mn_resp(i, j) = mean(squeeze(res.y(i,j,:)))*res.sdfs.n(i, j);
    end
end


% mean response in spikes/sec
mn_resp = (sum(sum(mn_resp)) + sum(mn_resp_extra) )/(sum(sum(N)) + sum(n_extra));
mn_resp = mn_resp/ex.setup.refreshRate;


% added by CL on Feb, 17th
if isfield(ex.stim.vals,'RCperiod') && ex.stim.vals.RCperiod>1
    mn_resp = mn_resp / ex.stim.vals.RCperiod;
end



% format the output-------------------------------------------------------

fieldNames = {'calctime','nbad','figa','labela','labelb','bestdelay','figb','timeoff','alltriggers'};
for n = 1:length(fieldNames)
    if isfield(res,fieldNames{n})
        res = rmfield(res,fieldNames{n});
    end
end

% net spikes per frame
res.mean_resp = mn_resp;
res.netSpikesPerFrame = y1+mn_resp;
res.netSpikesPerFrameBlank = y1_extra+mn_resp;
res.or = res.sdfs.x;

end





function res = getNetSpks1D(res, ex)
% compute net spikes-------------------------------------------------------
y = res.sdfs.nspk./res.sdfs.n;
N = res.sdfs.n;
for n = 1:length(res.sdfs.extras)
    y = [y; res.sdfs.extras{n}.nspk/res.sdfs.extras{n}.n];
    N = [N; repmat(res.sdfs.extras{n}.n, 1, size(y,2))];
end

% mean # of spikes per frame across all frames (I believe)
% note: nspk gives us all spikes that are used in the analysis, i.e. spikes
% are counted several times because of overlapping analysis windows
mn = sum(y.*N)/sum(N);
y1 = y-mn;  % this now gives me the deviation around the mean response


% I couldn't wrap my head around how to compute the mean number of spikes/frame
% without going via the SDF, so I do this via the sdf;
for nd = 1:size(res.y,1)
    %             for nd = 1:size(res.y,2) %hn version
    mn_resp(nd) = mean(squeeze(res.y(nd,1,:)))*res.sdfs.n(nd);
end
idx = find([res.times]>=res.delays(1) & [res.times]<=res.delays(end));

% add values for extra stimuli (blank)
extra_x=[];
for n = 1: length(res.sdfs.extras)
    mn_resp =[mn_resp, mean(res.sdfs.extras{n}.sdf(idx)) * res.sdfs.extras{n}.n];
    extra_x(n) = n*1000+min(res.sdfs.x(1));
end
% mean response in spikes/sec
mn_resp = sum(mn_resp)/sum(N);
mn_resp = mn_resp/ex.setup.refreshRate;


% added by CL on Feb, 17th
if isfield(ex.stim.vals,'RCperiod') && ex.stim.vals.RCperiod>1
    mn_resp = mn_resp / ex.stim.vals.RCperiod;
end

% net spikes per frame
nspk = y1+mn_resp;
x = [res.sdfs.x; extra_x];

% format the output-------------------------------------------------------
fieldNames = {'calctime','nbad','figa','labela','labelb','bestdelay','figb','timeoff','alltriggers'};
for n = 1:length(fieldNames)
    if isfield(res,fieldNames{n})
        res = rmfield(res,fieldNames{n});
    end
end

res.mean_resp = mn_resp;
res.netSpikesPerFrame = nspk(1:length(nspk)-length(res.sdfs.extras));

res.netSpikesPerFrameBlank = nspk(length(nspk)-length(res.sdfs.extras)+1:end);
%         res.netSpikesPerFrameBlank = nspk(length(nspk)-length(res.sdfs.extras):end);

res.or = x(1:length(nspk)-length(res.sdfs.extras));

end



