function [raster, sdf, nsac, nspikes, spikes, times, details] = ...
    CL_trigsdf(Trials, width, times, varargin)

%[raster, sdf, nsac, nspikes, spikes, times, details] = ...
%     CL_trigsdf(Trials, width, times, varargin)
% %suppresses command line output
% 
%former
% [sdf, nsac, nspikes, spikes, times] = trigsdf(Trials, width, times, ...)
%
% trigsdf takes a vector list of Trials, spikes
% and calculates a spike density function smoohted with 
% a Gaussian of S.D. 'width'.
% 'times' is a vector of timestamp values at which the density
% function is to be evaluated.
%
% sdf returns the smoothed function.
% nsac returns the number of sweeps making up this average.
% nspikes is the total number of spikes used.
% spikes is an unsmoothed function (a PSTH constructed with a bin
% width of 0.1ms) from which a tradional PSTH can be built.
%
% trigsdf(Trials, width, times, 'exp')
% Uses exponential of time constant width (if flag = 'exp') for smoothing.
%
% In each Trial structure, a vector of Trigger times may be
% specified. In this case, the spikes are aligned relative to the
% trigger time. If Trials.Trigger does not exist, all the trials
% are aligned at t=0. If more than one trigger time is specified,
% then than trial will contribute N times to the average (if there
% are N Trigger times.
%
% trigsdf(Trials, width, 'period', x)
% builds a psth of length period, to show average periodic
% modulation over period x. Here, 'times' specifies the time range from
% each trial that will be included. The sdf returned is evaluated at  0.1ms
% intervals
%
% trigsdf(Trials, width, 'kernel', x)
% uses the kernel x for smoothing


%
% trigsdf(Trials, width)

% build up a PSTH with 0.1ms bins. Making the bins shorter might
% make this ultrafast.

%history
% 10/05/15  hn: extended to adjust spike times for tue
% 26/01/16  cl: extended to derive single stimuli trigered spikes, variable
%           is called 'raster', it is basically:
%           sum(raster) = spikes



skipevent = 0;
nskip = 0;
period = 0;
flag = 'gauss';
clipping = 0; %% if 1, remove a full kernel width at each end.
nvar = nargin - 3;
j = 1;
lfpmin = [];


while j  <= nvar
  str = varargin{j};
  if strmatch('kernel',str)
    j = j+1;
    kl = varargin{j};
  elseif strncmpi('events',str,3)
      for k = 1:length(Trials)
          Trials(k).Trigger = Trials(k).Events{:,2};
      end
  elseif strncmpi('clip',str,3)
      clipping = 1;
  elseif strmatch('lfpmin',str)
    j = j+1;
    lfpmin = varargin{j};
    j = j+1;
    samplerate = varargin{j};
  elseif strmatch('period',str)
    j = j+1;
    period = round(varargin{j});
  elseif strncmpi('skipev',str,5)
      skipevent= 1;
  elseif strncmpi('trig',str,3)
    j = j+1;
    triggers = varargin{j};
    for k = 1:length(Trials)
        Trials(k).Trigger = round(triggers(k).trigger);
    end
  elseif strmatch(str,{'gauss','box','raw','halfg','exp'})
      flag = str;
  end
  j = j+1;
end


if length(Trials) == 0
    sdf = 0;
    nsac = 0;
    nspikes = 0;
    spikes = [];
    return;
end

nbins = 1+times(end) - times(1);
if period
    spikes = zeros(period*2,1);
    raster = zeros(period*2,1);
else
    spikes = zeros(nbins,1);
    raster = zeros(nbins,1);
end

if(~isfield(Trials,'Trigger'))
  [Trials.Trigger] = deal(0);
end

nsac = 0;
allspks = [];
spkcount = [];
raster = [];    rr=0; % index for raster
%bin is (spk - trigger) - times(1), = spk - (times(1)+trigger)

for j = 1:length([Trials.Trial])
    if ~isempty(Trials(j).Trigger) & ~isempty(Trials(j).Spikes)
% if skipevent is set, then spikes occuring with the first stimulus period
% after and event are NOT included. - allows the mean PSTH uncontaminated
% by skips to be calculated
      if isempty(lfpmin)
          sidx = 1:length(Trials(j).Spikes);
      else
          lfptimes = Trials(j).lfpo + round(Trials(j).Spikes .* samplerate);
          sidx = find(Trials(j).LFP(lfptimes) > lfpmin(lfptimes));
      end
        if skipevent & ~isempty(Trials(j).Events)
            skiper = Trials(j).Events{1,2};
            if skiper >= times(1) & skiper <= times(end)
                nskip = 1;
            else
                nskip = 0;
            end
            idx = find(Trials(j).Spikes < skiper | Trials(j).Spikes > skiper+period);
            txs = repmat(round(times(1) + Trials(j).Trigger),length(idx),1);
            sxs = repmat(Trials(j).Spikes(idx),1,length(Trials(j).Trigger));
        else
            txs = repmat(round(times(1) + Trials(j).Trigger),length(Trials(j).Spikes),1);
            sxs = repmat(Trials(j).Spikes(sidx),1,length(Trials(j).Trigger));
        end
%used to use fix(sxs-txs), but now 1bin = 1 clock tick, no need.
        spk = sxs - txs;
        idx = find(spk > 0 & spk <= nbins);
%        spkcount(j) = length(idx);
        if period
%      idx = find(spk >= 0 & spk <= (times(end) - times(1)));
%leave fix in here - period may not be an integer
            spk = fix(mod(spk(idx),period))+1;
            spk = [spk; spk + period];
            idx = 1:length(spk);
        end
% This would become nbins * spk(idx)/duration
%    bins = ceil(nbins .* spk(idx) ./ nbins);
% 

% The for loop is faster than spikes(spk(idx)) = spikes(spk(idx))+1
% And
        if(~isempty(idx))
            for k = length(idx):-1:1
                spikes(spk(idx(k))) = spikes(spk(idx(k))) +1;
            end
            % added @CL 26.01.2016
            for k = 1:size(spk, 2)
                rr = rr+1;
                idx2 = spk(:,k) > 0 & spk(:,k) <= nbins;
                raster(spk(idx2, k), rr) = 1;
                
            end
            
        end
%        allspks = [allspks; spk(idx)];
        if(period)
            nsac = nsac + length(Trials(j).Trigger) * (times(end)-times(1))/period;
            if nskip
                nsac = nsac - nskip;
            end
        else
        end
        
        
        
    end
%Do this AFTER test for Trigger, Spikes. If spikes is empty, still count
%trigger
        nsac = nsac + length(Trials(j).Trigger);
        
   
   
end

details.spkc = spkcount;
details.allspks = allspks;
if(period)
%  spikes = [spikes(1:period); spikes(1:period)+period];
end


if(strmatch('exp',flag))
  x = 0:width*3;
  kl = exp(-x/width);
  rates = conv(spikes,kl);
  fullsdf = [rates(1:end-width*3) .* 10000/(width * nsac)];
elseif(strmatch('gauss',flag))
  x = -width*3:width*3;
  kl = exp(-(x.^2)/(2 * width^2));
  rates = conv(spikes,kl);
  if clipping
      fullsdf = [rates(6*width:end-6*width) .* 10000/(sqrt(2 * pi) * ...
          width * nsac)];
      clipping = 6 * width;
      tclip(1) = 3 * width;
  else
      fullsdf = [rates(3*width:end-3*width) .* 10000/(sqrt(2 * pi) * ...
          width * nsac)];
  end
elseif(strmatch('halfg',flag))
  x = -0:width*3;
  kl = exp(-(x.^2)/(2 * width^2));
  rates = conv(spikes,kl);
  
  fullsdf = [rates(1:end-3*width) .* 2 * 10000/(sqrt(2 * pi) * ...
						  width * nsac)];
elseif(strmatch('box',flag))
  kl = ones(width,1);
  rates = conv(spikes,kl);
  if clipping
      tclip(1) = width/2;
      clipping = width;
  end
  if nsac == 0
      fullsdf = [rates(width:end-width+1)];
  elseif clipping
      fullsdf = [rates(width:end-width+1)] .* 10000/(width * nsac);
                      clipping = width;
  else
      fullsdf = [rates(1:end-width+1) .* 10000/(width * nsac)];
  end
elseif(strmatch('kernel',flag))
  kl = varargin{1};
  rates = conv(spikes,kl);
  w = floor(length(kl)/2);
  fullsdf = [rates(w:end-w) .* 10000/nsac];
elseif(strmatch('raw',flag))
  fullsdf = spikes;  
end
if clipping
    idx = (times - times(1));
    idx = idx(find(idx > clipping)) -clipping;
    times = times(1) + tclip(1) + idx;
else
    idx = 1 + (times - times(1));
end
if period
    w = round(period/2);
    idx = [period:period+w period+1-w:period-1];
end

idx = idx(find(idx <= length(fullsdf)));
sdf = fullsdf(idx);
nspikes = sum(spikes);
% wt = conv(ones(1,length(spikes)-length(kl)+1),kl);
% %wt = wt(1:end-length(kl)+1);
% 
% irates = rates(length(kl)+1:end-length(kl)+1)/length(kl);
% 
% disp('mean sdf')
% mean(irates)
% disp('mean spike count')
% mean(spikes'.*wt/sum(wt)*length(irates)*length(spikes)./length(irates))

