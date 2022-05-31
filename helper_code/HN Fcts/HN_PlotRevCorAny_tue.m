function [result, Expt] = HN_PlotRevCorAny_tue(Expt, varargin)

% result = PlotRevCorAny(Expt....)
% Plots up data for reverse correlation experiments. By default, a half Gaussian
% window with sd 2ms is used to smooth the sdfs, the variance betweeen sdfs
% as a function of time is estimated and time slices through the tuning
% curves are shown from the time when variance reaches 10% of peak to the time
% it falls back to 10% of peak
%
% PlotRevCorAny(Expt, 'sdfw', x) sets the Gaussian sd to x ticks (x/10 ms).
%
% PlotRevCorAny(Expt, 'slices', x) shows tuning curves for each of the
% delay times in vector x. 
%
% PlotRevCorAny(Expt, 'times', x) evaluates the sdfs for times x (relative to stimsulus transition)
% (default -100:10:1500)
%
% For all of these arguments, times are in clock ticks, not ms 
%
% PlotRevCorAny(Expt, 'yval', Y) uses only stimuli for which the value of
% exp2 is in Y
%
% Second order plots come in several types:
% PlotRevCorAny(Expt, 'secondorder', [4 1]);   
%            is a standard second order plot, showing data for all
%            combinations
% PlotRevCorAny(Expt, 'secondorder', [4 n]);   
%            shows a reduced data set, showing data for all stimuli at t = 0,
%            for each stimulus at t = -1 whose index is n or greater.
%            this can be combined with a third element in the specification
%            vector to limit the number of comparions:
% PlotRevCorAny(Expt, 'secondorder', [4 n m]);
%            this shows comparisons for only m+1 values at t = -1. So
%
% PlotRevCorAny(Expt, 'secondorder', [4 3 0]); Shows only the responses for
% where the third stimulus preceded the stimulus at t = 0
%
% Other uses of the 'secondorder' vector
%  [5 n m] is much the same as 4, but selects a limited number of stimuli
%  at t=0, showing all precedints. So [4 3 0] shows all responses to the
%  third stimulus at t=0, with a separate line for each preceding stimulus
%  value.
%
%  [1 1]  for AC expts. Plot cases where preceding = AC and preceding =
%  Corr separately
%  [2 1] Another special case, plotting preceding disp < 0 and >=0
%  separately.
%  [3 1] Plot separate RCs for preceding x < current x and preceding x >
%  current x   [3 1 0] shows only precedingx < current, [3 2 0] shows only
%  preceding x > current.
%
% ...'collapse' ignores one dimension of a 2-D expt. 
%     ....'collapse',1 collapses over values of expt1, 2 collapses expt2
% Add uncorr to plot
% make plots for AC -> and COR ->
% get default slices from variance wrt time.
%
% For DCPD expt (CP for time varying disparity). use
% PlotRevCorAny(Expt,'yval',0,'psych'...

% hn notes
% change 06/20/07 hn
% x values are binned with 0.01 resolution

% history
% 2/14/16   hn: changed code to correct for timing bug with RCperiod in
%           VisStim versions preceding 1.0.20

times = -200:10:1600;
np = 1;
sdfw = 20;
plotsum = 0;
plotbest = 0;
slices = [500 600 700 800 900 1000];
xrange = [];
labela = 'ACRCFiga';
labelb = 'ACRCFigb';
secondorder = 0;
nloops = 1;
legendpos = 0;
result = [];
showplot = 1;
needvar = 0;
smtype = 'halfg';
autoslice = 1;
legendstyle = 0;
label_this_yval = 1;
interp = 0;
summarize = 0;
nmin = 30;
getpsych = 0;
npsych = 1;
profile = 0;
iskip = 0; % bruce had this set to 1
aloop = 1;
collapse = [0 0];
setyvals = [];
type = [];
btype = [];
extraidx = {};
minplottime = 0;
autonmin = 0;
noextras = 0;
latency=[]; %0.1ms
duration =[]; %0.1ms
all_fr=[];
all_spk=[];
fl_fr = [];
% figh1 = figure;
% figh2 = figure;
type = 'or_seq';
if profile
    tic;
end

showpcolor = 0;
findarg = {};
if ischar(Expt)  % Named a file, not passed an Expt structure
    file = Expt;
    load(file);
end

j = 1; 
while j < nargin
    if strncmpi(varargin{j},'sdfw',4)
        j = j+1;
        sdfw = varargin{j};
    elseif strncmpi(varargin{j},'autoslices',4)
        autoslice = 1;
    elseif strncmpi(varargin{j},'slices',4)
        j = j+1;
        slices = varargin{j};
        autoslice = 0;
    elseif strncmpi(varargin{j},'add',3)
        plotsum = 1;
    elseif strncmpi(varargin{j},'box',3)
        smtype = 'box';
        minplottime = times(1) + sdfw;
    elseif strncmpi(varargin{j},'collapse',5)
        j = j+1;
        collapse(varargin{j}) = 1;
    elseif strncmpi(varargin{j},'range',3)
        j = j+1;
        xrange = varargin{j};
    elseif strncmpi(varargin{j},'mint',3)
        j = j+1;
        minplottime = varargin{j};
    elseif strncmpi(varargin{j},'exp2',4)
        j = j+1;
        btype = varargin{j};
    elseif strncmpi(varargin{j},'exp',3)
        j = j+1;
        type = varargin{j};
    elseif strncmpi(varargin{j},'best',4)
        plotbest = 1;
    elseif strncmpi(varargin{j},'figa',4)
        j = j+1;
        labela = varargin{j};
    elseif strncmpi(varargin{j},'figb',4)
        j = j+1;
        labelb = varargin{j}; 
    elseif strmatch(varargin{j},'legendpos')
      legendpos = varargin{j+1};
      j = j+1;
    elseif strncmpi(varargin{j},'nmin',4)
      nmin = varargin{j+1};
      if ischar(nmin) & strmatch(nmin,'auto')
          autonmin = 1;
      end
      j = j+1;
    elseif strcmpi(varargin{j},'latency') % in0.1ms
        j = j+1;
        latency = varargin{j};
    elseif strcmpi(varargin{j},'duration') % in 0.1ms
        j = j+1;
        duration = varargin{j};
    elseif strncmpi(varargin{j},'noplot',5)
        showplot = 0;
    elseif strncmpi(varargin{j},'noextra',5)
        noextras = 1;
    elseif strncmpi(varargin{j},'pcolor',5)
        showpcolor = 2;
        if length(varargin) > j & isnumeric(varargin{j+1})
            showpcolor = varargin{j+1};
            j = j+1;
        end
    elseif strncmpi(varargin{j},'profile',5)
        profile = 1;
    elseif strncmpi(varargin{j},'psych',5)
        getpsych = 1;
    elseif strncmpi(varargin{j},'upchoice',3)
        getpsych = 2;
    elseif strncmpi(varargin{j},'downchoice',4)
        getpsych = 3;
    elseif strncmpi(varargin{j},'interp',5)
        interp = 1;
        if length(varargin) > j & isnumeric(varargin{j+1})
            interp = varargin{j+1};
            j = j+1;
        end
    elseif strncmpi(varargin{j},'summary',6)
        summarize = 1;
    elseif strncmpi(varargin{j},'secondorder',6)
        j = j+1;
        secondorder = varargin{j}(1);
        aloop = varargin{j}(2);
        if length(varargin{j}) > 2
            nloops = aloop+varargin{j}(3);
        elseif ismember(secondorder,[3 4 5])
            nloops = 0;
        elseif ismember(aloop,1)
            nloops = 2;
        end
    elseif strncmpi(varargin{j},'select',6)
        j = j+1;
        idx = eval(['find(' varargin{j} ');']);
        Expt.Trials = Expt.Trials(idx);
        size(Expt.Trials)
    elseif strncmpi(varargin{j},'skip',4)
        j = j+1;
        iskip = varargin{j};
    elseif strncmpi(varargin{j},'splitce',6)
        secondorder = 1;
        aloop = 1;
        nloops = 2;
    elseif strncmpi(varargin{j},'splitdx',6)
        secondorder = 2;
        aloop = 1;
        nloops = 2;
    elseif strncmpi(varargin{j},'times',5)
       j = j+1;
       times = varargin{j};
    elseif strncmpi(varargin{j},'ytrial',4)
        j = j+1;
        yv = Expt.Stimvals.et;
        for k = 1:length(Expt.Trials)
            yvals(k) = Expt.Trials(k).(yv)(end);
        end
        idx = find(ismember(yvals,varargin{j}));
        Expt.Trials = Expt.Trials(idx);
        size(Expt.Trials)
    elseif strncmpi(varargin{j},'trials',6)
        j = j+1;
        Expt.Trials = Expt.Trials(varargin{j});
    elseif strncmpi(varargin{j},'yval',4)
        j = j+1;
        setyvals = varargin{j};
    end
    j = j+1;
end

% added CL 23.06.2016 to avoid figure popping
if showplot
    figh1 = figure;
    figh2 = figure;
end

% select only completed trials
itr = find(abs([Expt.Trials.Reward])>0);
Expt.Trials = Expt.Trials(itr);

Expt = AdjustSpikeTimes(Expt);

if isfield(Expt.stim.vals,'RCperiod') && Expt.stim.vals.RCperiod>1
    Expt = Adjust4RCperiod(Expt);
end

%btype = 'phase_seq';


result.calctime(1) = now;
% colors =  mycolors;
% excolors = mycolors;
% excolors{4} = [0.6 0.4 0]; %% dotted yellow invisible
% colors = {colors{:} colors{:}}; % just in case second order makes many lines
if ~isfield(Expt.Header,'RCparams')
    if isfield(Expt.Header,'StoreErr')
        for j = 1:length(Expt.Trials)
            lens(j) = length(Expt.Trials(j).Start);
        end
        nf = max(lens);
        for j = 1:length(Expt.Trials)
            if(length(Expt.Trials(j).Start) < nf)
                Expt.Trials(j).Start(end:nf) = NaN;
            end
        end
        fprintf('NaN Paddding for Store Bug\n');
    else
        lens =[];
        for j = 1:length(Expt.Trials)
            lens(j) = length(Expt.Trials(j).Start);
            id = find(Expt.Trials(j).Start < 0);
            if ~isempty(id)
                Expt.Trials(j).Start(id) = NaN;
            end
        end
        if isempty(lens)
            disp ('no trials')
            return
        else
        nf = min(lens);
        end
        nf = floor(prctile(lens,10));
        Expt.Trials = Expt.Trials(find(lens >= nf));
        result.excluded_trials = find(lens<nf);
        for j = 1:length(Expt.Trials)
            Expt.Trials(j).Start = Expt.Trials(j).Start(1:nf);
        end
    end

    starts = [Expt.Trials.Start];
    if isfield(Expt.Header,'Name')
            Expt.Header.Name = strrep(Expt.Header.Name,'\','/');
    elseif isfield(Expt.Header,'fileName')
            Expt.Header.Name = strrep(Expt.Header.fileName,'\','/');
    else Expt.Header.Name = '';
    end

if isempty(type)
    if ischar(Expt.Stimvals.et)
        type = Expt.Stimvals.et;
    else
        if Expt.Stimvals.et == 263
            type = 'dO';
        elseif Expt.Stimvals.et == 237
            type = 'dO';
        end
    end
end
if isempty(btype)
    if isfield(Expt,'Stimvals') 
        if ischar(Expt.Stimvals.e2)
            btype = Expt.Stimvals.e2;
        else
            if Expt.Stimvals.e2 == 112
                btype = 'ce';
            elseif Expt.Stimvals.e2 == 117
                btype = 'ce';
            end
        end
    end
end
    if strcmp(btype,'Pd') & strcmp(type,'Dc')
        type = 'Pd';
        btype = 'Dc';
    end
    if strcmp(btype,'or') & strcmp(type,'Dc')
        type = 'or';
        btype = 'Dc';
    end

    if strcmp(btype,'dx') & strcmp(type,'Dc')
        type = 'dx';
        btype = 'Dc';
    end

    if strcmp(btype,'e0') | collapse(2)
        if ~strcmp(Expt.Stimvals.e3,'e0')
            btype = Expt.Stimvals.e3;
        else
            btype = [];
        end
    end
if strcmp(btype,'ph')
    btype = [];
end

    if strcmp(btype,'ce')
        linestyles = {'--', '-', ':','-',':',};
        label_this_yval = 2;
    elseif getpsych == 1
        linestyles = {'--', '-', ':','-.','-'};
        label_this_yval = 2;
    else
        linestyles = {'-', '--', ':','-.','-'};
    end


    %
    % first make sure all trials have the same length vectors
    % describing the stimulus

    res.spksum = 0;
    for j = 1:length(Expt.Trials)
        if isempty(type)
            type = 'dx';
            disp('type set to dx')
        end
        Expt.Trials(j).eval = Expt.Trials(j).(type)(end);
        len(j) = eval(['length(Expt.Trials(j).' type ');']);
        res.spksum = res.spksum + length(Expt.Trials(j).Spikes);
    end

    if res.spksum == 0
        fprintf('No Spikes in %s\n',Expt.Header.Name);
        return;
    end

    if isfield(Expt.Header,'StoreErr')
        tlen = max(len);
        Nf = tlen;
    else
        tlen = floor(mean(len)+0.1);
        if nf < tlen & tlen -nf < tlen/10;
            tlen = nf;
        end
%what is stored in Trials.Nf is the number of frames completed _before_
%the last frame (becuase it is send then). So Add one.
        if isfield(Expt.Trials,'Nf')
            Nf = 1+round(mean([Expt.Trials.Nf]));
        elseif isfield(Expt,'Stimvals')
            Nf = Expt.Stimvals.Nf;
        else 
            for n = 1:length(Expt.Trials);
                nfs(n) = length([Expt.Trials(n).Start]);
            end
            Nf = round(mean(nfs));
        end
        if Nf < tlen & Nf > 0
            tlen = Nf;
        end
    end

    if profile
        fprintf('Filling trials at %.2f\n',toc);
    end
    if isfield(Expt,'Stimvals')
        if isfield(Expt.Stimvals,'fz') & ~isempty(Expt.Stimvals.fz)
            frameperiod = 10000/Expt.Stimvals.fz;
        else
            frameperiod = 104;
        end
    else frameperiod = 1/Expt.setup.refreshRate*10000;
    end
    if ~isfield(Expt,'Stimvals')
%         fprintf('2')
        for n = 1:length(Expt.Trials)  
            % why do we need this?--- old code
 %           if length(eval(['Expt.Trials(n).' type ]))>length(Expt.Trials(n).Start);
%                 eval(['Expt.Trials(n).' type ...
%                     '(1:length(Expt.Trials(n).Start),1)= [Expt.Trials(n).' ...
%                     type '(1:length(Expt.Trials(n).Start))];']);
%------------------------------------------------------------------
                % now new: make sure we have column vectors
                if eval(['size([Expt.Trials(n).' type '],2)'])>1
                eval(['Expt.Trials(n).' type ...
                    '= [Expt.Trials(n).' ...
                    type '(1:length(Expt.Trials(n).Start))];']);
                else
                eval(['Expt.Trials(n).' type ...
                    '= transpose([Expt.Trials(n).' ...
                    type '(1:length(Expt.Trials(n).Start))]);']);
                end
%            end
            Expt.Trials(n).Trial = n;
            %Expt.Trials(n).Start = (Expt.Trials(n).Start'-Expt.Trials(n).TrialStart)*10000;
        end
    end

    for j = 1:length(Expt.Trials)
        if ~isempty(btype) & length(Expt.Trials(j).(btype)) == 1
            Expt.Trials(j).(btype)(1:tlen,1) = Expt.Trials(j).(btype);
        end
        if len(j) < tlen
            Expt.Trials(j).(type)((len(j)+1):tlen) = NaN;
            if isfield(Expt.Trials,'ce')
                Expt.Trials(j).ce((len(j)+1):tlen) = NaN;
            end
            if ~isempty(btype)
                Expt.Trials(j).(btype)((len(j)+1):tlen) = NaN;
            end
            if isfield(Expt.Trials,'st_seq')
                Expt.Trials(j).st_seq((len(j)+1):tlen) = NaN;;
            end
            if isfield(Expt.Trials,'me')
                Expt.Trials(j).me((len(j)+1):tlen) = NaN;
            end
        elseif len(j) > tlen
            Expt.Trials(j).(type) = Expt.Trials(j).(type)(1:tlen);
            if ~isempty(btype)
                Expt.Trials(j).(btype) = Expt.Trials(j).(btype)(1:tlen);
            end
            if isfield(Expt.Trials,'ce') & length(Expt.Trials(j).ce) > 1
                Expt.Trials(j).ce  = Expt.Trials(j).ce(1:tlen);
            end
            if isfield(Expt.Trials,'st_seq') && length(Expt.Trials(1).st_seq)>1
                Expt.Trials(j).st_seq = Expt.Trials(j).st_seq(1:tlen);
            end
            if isfield(Expt.Trials,'me') & length(Expt.Trials(j).me)>=tlen
                Expt.Trials(j).me = Expt.Trials(j).me(1:tlen);
            elseif isfield(Expt.Trials,'me') & length(Expt.Trials(j).me) ==1
                Expt.Trials(j).me = ones(tlen,1)*Expt.Trials(j).me;
            end
        end
        id = find(isnan(Expt.Trials(j).Start));
        id = id(find(id <= tlen));
        for k = 1:length(id)
            Expt.Trials(j).(type)(id(k)) = Expt.Trials(j).(type)(id(k)-1);
        end
        id = find((Expt.Trials(j).(type)) > 1000);
        Expt.Trials(j).(type)(id) = -10000;
        tstart = Expt.Trials(j).Start(1) - Expt.Trials(j).TrialStart;
        tend = Expt.Trials(j).Start(tlen) - Expt.Trials(j).TrialStart;
        Expt.Trials(j).Count = sum(Expt.Trials(j).Spikes > tstart & Expt.Trials(j).Spikes < tend);
        Expt.Trials(j).TrialDur = tend-tstart;
    end
    rcparams.btype = btype;
    rcparams.type = type;
else
    btype = Expt.RCparams.btype;
    type = Expt.RCparams.type;
    starts = [Expt.Trials.Start];
end


if(autoslice)
    needvar = 1;
end


if profile
    fprintf('Done at %.2f\n',toc);
end

xvals = unique([Expt.Trials.(type)]);
    
xvals = xvals(find(~isnan(xvals)));
if strcmpi(type,'dx') | strcmpi(type,'Pd')|strcmpi(type,'dO') % hn change: bin data with 0.01 resolution
    xvals = round(xvals*100)/100;
    xvals = unique(xvals);
end
if ~isempty(btype)
    Expt.Trials(1).(btype);
    yid = find(~isnan([Expt.Trials.(btype)]));
    ymat = [Expt.Trials.(btype)];
    %    yvals = unique([Expt.Trials.(btype)](yid));
    yvals = unique(ymat(yid));
else
    yvals = [];
end

if ~isempty(setyvals)
    if sum(yvals == setyvals) == 0 %% no  matches
        [a,b] = min(abs(yvals - setyvals(1)));
        yvals = yvals(b);
    else
        yvals = setyvals;
    end
end

if ~isempty(xrange)
    idx = find(xvals >= xrange(1) & xvals <= xrange(2));
    xvals = xvals(idx);
end

if nloops == 0
    if ismember(secondorder, [4 5])
        nloops = length(xvals);
    elseif secondorder == 3
        nloops = aloop+1;
    else
        nloops = 0;
    end
end


endoffset = 5; %% 0.5ms
if isfield(Expt,'RCparams')
    tx = [Expt.Trials.(type)];
    ty = [Expt.Trials.(btype)];
    durs = [Expt.Trials.dur];
else
    for j = 1:length(Expt.Trials)
        a = Expt.Trials(j).(type);
        tx(:,j) = a(1,:)'; % total hack
        if ~isempty(yvals)
            ty(:,j) = Expt.Trials(j).(btype);
        end
        if isfield(Expt.Trials(j),'End')
        durs(:,j) = [diff(Expt.Trials(j).Start); endoffset + Expt.Trials(j).End(end) - Expt.Trials(j).Start(end)];
        Expt.Trials(j).dur = durs(1:tlen,j);

        else
            durs(:,j) = [diff(Expt.Trials(j).Start);diff([Expt.Trials(j).Start(end) Expt.Trials(j).Start(end-1)])];
        Expt.Trials(j).dur = durs(1:tlen,j);
        end
    end
end


% fprintf(' start: %1.2f \n',Expt.Trials(2).Start(1))
%
%
% check Stimulus durations and make sure that they are all the same.
% Exclude any outliers (more that 30% of frame period away from mean)
%

if isfield(Expt,'RCparams')
    stimdur = Expt.RCparams.stimdur;
else
    if profile
        fprintf('getting duriation at %.2f\n',toc);
    end
    stimdur = mode(durs(:));
    rcparams.stimdur = stimdur;
    if profile
        fprintf('Done at %.2f\n',toc);
    end
end

result.nbad = length(find(abs(durs(2:end,:) - stimdur) > frameperiod/3));
result.nframes = prod(size(tx));
result.frameperiod = frameperiod;
result.stimdur = stimdur;
result.Nf = Nf;
result.Stimvals.et = type;
result.Stimvals.e2 = btype;
result.sdfw = sdfw;
    
badtimes = find(abs(durs - stimdur) > frameperiod/3 | isnan(durs));
good(badtimes) = 0;
good(1,:) = 1;
good = ones(size([Expt.Trials.dur]));
if(iskip)
    good(1:iskip,:) = 0;
end

if autonmin
    txvals = unique(tx(:));
    tyvals = unique(ty(:));
    nstim = length(txvals) * length(tyvals);
    nmin = sum(good(:))./(nstim * 3);
    if getpsych
        nmin = nmin/3;
    end
end

result.nmin = nmin;

if profile
    fprintf('Start Respdir at %.2f\n',toc);
end

if getpsych & isfield(Expt.Trials,'RespDir')
    respdirs = [Expt.Trials.RespDir];
    id = find(respdirs == -1);
    %    good(:,id) = 1;
    for j = 1:length(respdirs)
        if respdirs(j) == 1
            good(find(good(:,j)),j) = 2;
        elseif respdirs(j) == -1
            good(find(good(:,j)),j) = 1;
        elseif respdirs(j) == 0
            good(:,j) = 0;
        end
    end
end

if summarize
    result.xv = unique(tx(find(~isnan(tx))));
    for j = 1:length(result.xv)
        result.nx(j) = length(find(tx == result.xv(j)));
    end
    if exist('ty')
        result.yv = unique(ty(find(~isnan(ty))));
    end
end


if showplot
    result.figa = figh2;
    result.labela = labela;
    result.labelb = labelb;
    hold off;
end
h = [];

if size(xvals,1) > size(xvals,2)
    xvals = xvals';
end
if isempty(yvals)
    yvals = 0;
    ty = zeros(size(tx));
end

sdfs.extras = {};
isextra = zeros(size(tx));
nextra = 0;
nstims = prod(size(tx));
if isfield(Expt.Trials,'ce') & size([Expt.Trials.ce],1) > 1
    idx = find([Expt.Trials.ce] == 0);
    if ~isempty(idx) & length(idx) < nstims/length(yvals)
        nextra = nextra+1;
        extraidx{nextra} = idx;
        extralabel{nextra} = 'Uncorr';
    else
        uidx = idx;
    end
end

if isfield(Expt.Trials,'st_seq')
    idx = find([Expt.Trials.st_seq] == 0);
    if ~isempty(idx)
        nextra = nextra+1;
        extraidx{nextra} = idx;
        extralabel{nextra} = 'Blank';
    end
end

% if isfield(Expt.Trials,'me')
%     idx = find([Expt.Trials.me] == -1);
%     if ~isempty(idx)
%         nextra = nextra+1;
%         extraidx{nextra} = idx;
%         extralabel{nextra} = 'Left';
%     end
%     idx = find([Expt.Trials.me] == 1);
%     if ~isempty(idx)
%         nextra = nextra+1;
%         extraidx{nextra} = idx;
%         extralabel{nextra} = 'Right';
%     end
% end
if noextras
    nextra = 0;
    extraidx = {};
end
allextras = [];
for j = 1:nextra
    allextras = [allextras; extraidx{j}];
    isextra(extraidx{j}) = j;
end

tmpy = ty;
tmpy(allextras) = NaN;
yvals = unique(tmpy(find(~isnan(tmpy))));
if isempty(yvals) %% all seem to be extras
    yvals = unique(ty);
end
if ~isempty(setyvals)
    if sum(yvals == setyvals) == 0 %% no  matches
        [a,b] = min(abs(yvals - setyvals(1)));
        yvals = yvals(b);
    else
        yvals = setyvals;
    end
end

if ~isfield(Expt,'RCparams')
    Expt.RCparams = rcparams;
end

loopctr = 1;

if profile
    fprintf('Structure ready at %.2f\n',toc);
end


%
% good((i,j) == 2 means RespDir == 1, means Downward eye movement.
% good((i,j) == 1 means RespDir == -1, means upward eye movement.
%

if size(good,1) > size(tx,1)
    good = good(1:size(tx,1),:);
end
ntg = 1;
if getpsych == 1
    respdirs = [-1 1];
    npsych = 2;
elseif getpsych == 2
    findarg = {'RespDir' 1};
    npsych = 1;   
elseif getpsych == 3
    findarg = {'RespDir' 2};
    respdirs = [1 -1];
    npsych = 1;  
end
pid = find(times > minplottime);

for j = 1:length(Expt.Trials)
    alltriggers{j} = [];
end
for loop = aloop: nloops;
    nx = 1;
    if secondorder
        findarg{1} = 'secondorder';
        findarg{2}(1) = secondorder;
        findarg{2}(2) = loop;
        if ismember(secondorder,[4 5])
            findarg{2}(3) = xvals(loop);
        end
    end
    for x = xvals;
        ny = 1;
        for iy = 1:length(yvals) * npsych %% not fussy about row/column vector

            y = yvals(mod(iy,length(yvals))+1);
            if npsych > 1
                choiceval = ceil(iy/length(yvals));
                [Expt, tidx, n,all_fr,all_spk,fl_fr] = FindTrials(Expt, x, y, tx, ty, good, ...
                    isextra, findarg{:},'RespDir',choiceval,'latency',latency,...
                    'duration',duration,'all_fr',all_fr,'all_spk',all_spk,'fl_fr',fl_fr);
            else
                [Expt, tidx, n,all_fr,all_spk,fl_fr] = FindTrials(Expt, x, y, tx, ty, good,...
                    isextra, findarg{:},'latency',latency,'duration',...
                    duration,'all_fr',all_fr,'all_spk',all_spk,'fl_fr',fl_fr);
                choiceval = 1;
            end
            %tidx is a list of Trials, containing at least one of the desired stim
            %types. n is the total number of stim presentations.
            for k =1:length(tidx)
                    alltriggers{tidx(k)} = [ alltriggers{tidx(k)}; Expt.Trials(tidx(k)).Trigger'];
            end
            if ~isempty(tidx) & (n > nmin);
                [sdf, n,nspk, psth] = HN_trigsdfa_hn(Expt.Trials(tidx),sdfw,times,smtype);
                sdfs.x(nx, ny, loopctr) = x;
                sdfs.y(nx, ny, loopctr) = y;
                sdfs.n(nx, ny, loopctr) = n;
                sdfs.nspk(nx,ny,loopctr) = nspk;
                if secondorder
                    sdfs.z(nx, ny, loopctr) = xvals(loop);
                end
                sdfs.s{nx, ny, loopctr} = sdf;
                sdfs.psth{nx,ny,loopctr} = psth;

                if ~showplot
                    
                elseif ny > 1
                    if( ~plotsum)
                        hn = plot(times(pid)/10,sdf(pid),'color',colors{nx},'linestyle',linestyles{mod(ny-1,4)+1});
                        hold on;
                    end
                elseif ny == 1
                    if ~plotsum
                        hn = plot(times(pid)/10,sdf(pid),'color',colors{nx},'linestyle',linestyles{loopctr});
                        hold on;
                    end
                elseif  ~isempty(tidx)
                    if ~plotsum
                        hn = plot(times(pid)/10,sdf(pid),'--','color','k');
                        hold on;
                    end
                end
                if (ny == label_this_yval | legendstyle ==1) & showplot
                    h(np) = hn;
                    if secondorder
                        labels{np} = sprintf('%.2f->%.2f n = %d',xvals(loop),x,n);
                    else
                        if choiceval > 1
                            labels{np} = sprintf('%.2f n = %d,%d',x,n,sdfs.n(np));
                        else
                            labels{np} = sprintf('%.2f n = %d',x,n);
                        end
                    end
                    np = np+1;
                end
                ny = ny+1;
            else
%                 if ny == 2
%                     disp('in ny ==2')
%                     sdfs.n(nx,ny) = n;
%                 end
            end
        end
        %ny only incremented if data found. ny == 1 means no data found.
        if ny > 1
            nx = nx+1;
        end
    end
    loopctr = loopctr+1;
end



if profile
    fprintf('Main loop done at %.2f\n',toc);
end
nex = 1;
goodextras = 0;
for loop = aloop: nloops;
for j = 1:nextra * npsych
    if secondorder
        findarg{1} = 'secondorder';
        findarg{2}(1) = secondorder;
        findarg{2}(2) = loop;
        if ismember(secondorder,[4 5])
            findarg{2}(3) = xvals(loop);
        end
    end

    if npsych > 1
        lid = 1+ mod(j-1,npsych);
        choiceval = ceil(j/nextra);
        [Expt, tidx, n,all_fr,all_spk,fl_fr] = FindTrials(Expt, x, y, tx, ty, good, isextra, ...
            'extra', lid, findarg{:},'RespDir',choiceval,...
            'latency',latency,'duration',duration,'all_fr',all_fr,...
            'all_spk',all_spk,'fl_fr',fl_fr);
    else
        [Expt, tidx, n,all_fr,all_spk,fl_fr] = FindTrials(Expt, x, y, tx, ty, good, isextra,...
            'extra', j, findarg{:},'latency',latency,'duration',duration,...
            'all_fr',all_fr,'all_spk',all_spk,'fl_fr',fl_fr);
        lid = j;
        choiceval = loop;
    end
    for k =1:length(tidx)
        alltriggers{tidx(k)} = [ alltriggers{tidx(k)}; Expt.Trials(tidx(k)).Trigger'];
    end
    if ~isempty(tidx) & (n > nmin);
        if choiceval == 1
            goodextras = goodextras+1;
        end
        [sdf, n,nspk ,psth] = HN_trigsdfa_hn(Expt.Trials(tidx),sdfw,times,smtype);
        if showplot
            h(np+nex-1) = plot(times(pid)/10,sdf(pid),':','color',excolors{nex - (choiceval-1)*goodextras},'linestyle',linestyles{3+choiceval-1},'linew',2);
            labels{np+nex-1} = [extralabel{lid} sprintf(' n = %d',n)];
        end
        sdfs.extras{lid,choiceval}.sdf = sdf;
        sdfs.extras{lid,choiceval}.n = n;
        sdfs.extras{lid,choiceval}.label = extralabel{lid};
        sdfs.extras{lid,choiceval}.nspk = nspk;
        sdfs.extras{lid,choiceval}.psth = psth;
        nex = nex+1;
    end
end
np = np+nex;
end

result.alltriggers = alltriggers;
if profile
    fprintf('Extras done at %.2f\n',toc);
end

np = 1;

bestvar = 0;



if (plotbest  | needvar) & isfield(sdfs,'s')
    nv = 1;
    id = find(times < 200);
    bestnv = NaN;
    for j = id(end):2:length(sdfs.s{1,1})
        x = [];
        y = [];
        z = [];
        for k = 1:size(sdfs.s,1)
            x(k) = sdfs.x(k,1);
            for co = 1:size(sdfs.s,2)
                y(k,co) = 0;
                nc = size(sdfs.s,3);
% for second order stimuli, just use first order kernels for variance
% estimate, so collapse across dimension 3 of sdfs
                for zi = 1:nc
                    if length(sdfs.s{k,co,zi}) >= j
                        y(k,co) = y(k,co) + sdfs.s{k,co,zi}(j)/nc;
                    end
                end
            end
        end
        dvar = sum(var(y));
        vars(nv) = dvar;
        delays(nv) = times(j);
        sampleid(nv) = j;
        result.x(:,nv) = x;
        for co = 1:size(sdfs.s,2)
            result.y(:,:,nv) = y;
        end
        if(dvar > bestvar)
            bestvar = dvar;
            bestj = j;
            bestx = x;
            besty = y;
            bestnv = nv;
        end
        nv = nv + 1;
    end
    result.delays = delays;
    result.bestdelay = bestnv;
    result.vars = vars;
    result.timeoff = times(1);
    if showplot & plotbest & bestnv > 0
        figure(figh2)
        hold on;
        for co = 1:size(sdfs.s,2)
            result.h(co) = plot(bestx,besty(:,co),'-','color',colors{np});
            hold on;
            np = np+1;
        end
        title(sprintf('%s at %d',Expt.Header.Name,bestj));
    end
    if showplot
        figure(figh2)
        plot(delays./10,sqrt(vars).*3,'k','linewidth',2);
    end
else
    bestnv = 1;
end

w = 100/(times(2)-times(1));

if exist('vars') & w > length(vars)/2;
    w = floor(length(vars)/2);
end
if bestnv <= w
    result.varratio(1) = 1;
    result.varratio(2) = 0;
    result.varratio(3) = 0;
elseif bestnv+w > length(vars)
    result.varratio(2) = mean(vars(end-2*w:end));
    result.varratio(3) = mean(vars(1:w*2));
    result.varratio(1) = result.varratio(2)/result.varratio(3);
elseif isnan(bestnv)
    result.varratio(1:3) = NaN;
else
    result.varratio(2) = mean(vars(bestnv-w:bestnv+w));
    result.varratio(3) = mean(vars(1:w*2));
    result.varratio(1) = result.varratio(2)/result.varratio(3);
end
if bestnv > w
    res.varratio(4) = prctile(vars,30);
end
if showpcolor
    hold off;
    startt = find(times > sdfw/2);
    PlotPcolor(sdfs,times/10,[],interp,startt(1));
    h = [];
end
if showplot
    figure(figh1)
    if(legendpos < 6 & ~isempty(h))
        hid = find(ishandle(h) & h > 0);
        legend(h(hid),labels{hid},legendpos);
    end
    title(sprintf('%s V%.1f',Expt.Header.Name,result.varratio(1)));
    ylabel('Rate (sp/s)');
    xlabel('Time (ms)');

    figure(figh1)
    hold on;
end
h = [];
labels = {};
if ~plotbest & isfield(sdfs,'s') & bestnv > 0
    if showplot
        result.figb = figh1;
    end
    if autoslice
        [mv, mj] = max(result.vars);
        minv = min(result.vars);
        th = minv + (mv-minv)/10;
        id = find(result.vars > th);
        step = range(result.delays(id))/6;
        slices = round(result.delays(id(1)):step:result.delays(id(end)));
        result.slices = slices;
        step = range(id)/6;
        slices = round(id(1):step:id(end));
        %cslices are coarse time samples, at the peak and at half height.
        id = find(result.vars > minv +(mv-minv)/2);
        cslices(2) = result.bestdelay;
        cslices(1) = id(1);
        cslices(3) = id(end);
    else
        %
        %
        result.slices = slices;
        for j = 1:length(slices)
            [diffs, id] = min(abs(slices(j) - times));
            sliceid(j) = id;
            sampleid(id) = id;
            delays(id) = times(id);
        end
        slices = sliceid;
    end

    if profile
        fprintf('Slices ready at %.2f\n',toc);
    end

    %
    % id() keeps track of the data sample points where var was evaluated.
    % so sampleid(id(k)) is the data sample in the sdf associated with vars(k),
    % at time delays(id(k))
    %
    %
    if showplot
        subplot(1,1,1);
        hold off;
        for j = slices;
            x = [];
            y = [];
            z = [];
            sn = sampleid(j);
            for k = 1:size(sdfs.s,1)
                x(k) = sdfs.x(k,1);
                for co = 1:size(sdfs.s,2)
                    if isempty(sdfs.s{k,co})
                        y(k,co) = NaN;
                    else
                        y(k,co) = sdfs.s{k,co}(sn);
                    end
                end
            end
            for co = 1:size(y,2)% size(sdfs.s,2)
                if size(linestyles,2)>=co
                h(np) = plot(x,y(:,co),'color',colors{np},'linestyle',linestyles{co});
                else
                h(np) = plot(x,y(:,co),'color',colors{np});
                end
                hold on;
            end
            for co = 1:size(sdfs.extras)
                if ~isempty(sdfs.extras{co})
                    plot([min(x) max(x)],[sdfs.extras{co}.sdf(sn) sdfs.extras{co}.sdf(sn)],'--','color',colors{np});
                end
            end
            labels{np} = sprintf('dT = %.0fms',delays(j)/10);
            np = np+1;
        end
        if legendpos < 6
            legend(h,labels,0);
        end
        if getpsych
            title(sprintf('%s Solid = Up',Expt.Header.Name));
        else
            title(sprintf('%s',Expt.Header.Name));
        end
        if secondorder ==4
            zmax = 0;
            if showpcolor == 2
                subplot(2,2,1);
                hold off;
                zmax(1) = PlotPcolor(sdfs, times, result.delays(cslices(2))-300,interp);
                hold on;
                subplot(2,2,2);
                zmax(2) = PlotPcolor(sdfs, times, result.delays(cslices(1)),interp);
                subplot(2,2,3);
                zmax(3) = PlotPcolor(sdfs, times, result.delays(cslices(2)),interp);
                subplot(2,2,4);
                hold off;
                zmax(4) = PlotPcolor(sdfs, times, result.delays(cslices(3)),interp);
                for j = 1:4;
                    subplot(2,2,j);
                    caxis([0 max(zmax)]);
                    colorbar;
                end
            elseif showpcolor == 1
                hold off;
                PlotPcolor(sdfs, times, result.delays(result.bestdelay),interp);
            else
                subplot(1,1,1);
                h = [];
                labels = {};
                hold off;
                if(autoslice)
                    sn = sampleid(cslices(2));
                else
                    sn = sampleid(slices);
                end
                for gr = 1:length(sn)
                    if length(sn) > 2
                        subplot(2,2,gr);
                    elseif length(sn) > 1
                        subplot(2,1,gr);
                    else
                        subplot(1,1,1);
                    end
                    for pre = 1:size(sdfs.s,3)
                        for k = 1:size(sdfs.s,1)
                            x(k) = sdfs.x(k,1,pre);
                            for co = 1:size(sdfs.s,2)
                                if isempty(sdfs.s{k,co,pre})
                                    y(k,co) = NaN;
                                else
                                    y(k,co) = sdfs.s{k,co,pre}(sn(gr));
                                end
                            end
                        end
                        for co = 1:size(sdfs.s,2)
                            h(pre) = plot(x,y(:,co),'color',colors{pre},'linestyle',linestyles{co});
                            hold on;
                        end
                        labels{pre} = sprintf('pre = %.3f',sdfs.z(1,1,pre));
                    end
                    if(legendpos < 6 & ~isempty(h) && gr == 1)
                        legend(h,labels,legendpos);
                    end
                    title(sprintf('At %d',times(sn(gr))));
                end
            end
        end
        

        if profile
            fprintf('Done at %.2f VR %.1f,%.1f,%.1f\n',toc,result.varratio(1:3));
        end
    end
end
if isfield(sdfs,'s')
    result.sdfs = sdfs;
    result.times = times;
    result.delaysamples = sampleid;
end
result.all_fr = all_fr;
result.all_spk = all_spk;
result.fl_fr = fl_fr;
for j = 1:length(extraidx)
    tx(extraidx{j}) = -(1000+j);
end
for j = 1:length(Expt.Trials)
    Expt.Trials(j).Pd = tx(:,j);
end
result.calctime(2) = now;
% Real end of PlotRevCorAny function
function zmax = PlotPcolor(sdfs, times, t,interp,startt)

if nargin < 5
    startt = 10;
end

if isempty(t)   %% Plot first order resp w.r.t t
for j = 1:size(sdfs.x,1)
    for k = startt:length(times)
        ik = k - startt+1;
        Y(j,ik) = sdfs.x(j,1,1);
        X(j,ik) = times(k);
        Z(j,ik) = sdfs.s{j,1,1}(k);
    end
end
else
id = find(times == t);
for j = 1:size(sdfs.x,1)
    for k = 1:size(sdfs.x,3)
        X(j,k) = sdfs.x(j,1,k);
        Y(j,k) = sdfs.z(j,1,k);
        if id > length(sdfs.s{j,1,k})
            Z(j,k) = NaN;
        else
            Z(j,k) = sdfs.s{j,1,k}(id);
        end
    end
end
end
[X,Y,Z] = fillpmesh(X,Y,Z);
pcolor(X,Y,Z);
colormap('hot');
if interp
    shading('interp');
else
    shading('flat');
end

title(sprintf('At %d ms',t/10));
zmax = max(caxis);

function [Expt, tidx,ns,all_fr,all_spk,fl_fr] = FindTrials(Expt, x, y, tx, ty, good, isextra, varargin)

secondorder =  0;
loop = 0;
extra = 0;
gval = 1;
j =1;
latency =[];
duration =[];
while j < nargin -6
    if strncmpi(varargin{j},'secondorder',4)
        j = j+1;
        loop = varargin{j}(2);
        secondorder = varargin{j}(1);
        if ismember(secondorder, [4 5])
            xval = varargin{j}(3);
        end
    elseif strncmpi(varargin{j},'xvals',4)
        j = j+1;
        xvals = varargin{j};
    elseif strncmpi(varargin{j},'RespDir',5)
        j = j+1;
        gval = varargin{j};
    elseif strcmpi(varargin{j},'latency')
        j = j+1;
        latency = varargin{j};
    elseif strcmpi(varargin{j},'duration')
        j = j+1;
        duration = varargin{j};
    elseif strcmpi(varargin{j},'all_fr')
        j = j+1;
        all_fr = varargin{j};
    elseif strcmpi(varargin{j},'all_spk')
        j = j+1;
        all_spk = varargin{j};
    elseif strcmpi(varargin{j},'fl_fr')
        j = j+1;
        fl_fr = varargin{j};
    elseif strncmpi(varargin{j},'extra',4)
        j = j+1;
        extra = varargin{j};
    end
    j = j+1;
end



tidx = [];
ns = 0;
% change by hn 06/20/07: bin x-values with  0.001 resolution
tx = round(tx*100)/100;
for j = 1:length(Expt.Trials)
    if extra
        idx = find(isextra(:,j) == extra & good(:,j) == gval);
        if secondorder == 1
            if loop == 1
                fidx = find(ty(:,j) == 1);
            else
                fidx = find(ty(:,j) == -1);
            end
            idx = intersect(idx,fidx+1);
        end
    else
% if secondorder ==1, build RC using only frames preceded by
% negative (fy == -1) one where the preceding frame was positive (fy==1)
% i.e. for AC expts, look at resp to Corr when preceded by Corr or AntiCorr
        if secondorder == 1
            if loop == 1
                fidx = find(ty(:,j) == 1);
            else
                fidx = find(ty(:,j) == -1);
            end
            idx = find(abs(tx(:,j)- x)<0.01 & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
% if secondorder == 2, build 2 RCs, one where the preceding frame had
% x < 0 and y ==1 ie crossed correlated disparity.
        elseif secondorder == 2
            if loop == 1
                fidx = find(tx(:,j)  < 0 & ty(:,j) ==1);
            else
                fidx = find(tx(:,j) > 0 & ty(:,j) ==1);
            end
            idx = find(abs(tx(:,j)-x)<0.01 & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
% if secondorder == 3, build RC only where the preceding frame had
% x < current (loop = 1) or x > current (loop = 2)
        elseif secondorder == 3
            if loop == 1
                fidx = find(tx(:,j)  < x);
            elseif loop == 2
                fidx = find(tx(:,j) > x);
            elseif loop == 3
                fidx = find(tx(:,j) == 0);
            end
            idx = find(abs(tx(:,j)- x)<0.01 & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
% if secondorder == 4, build RC for all preceding x vals;
% if aloop is set to n, build only for precding values xvals[aloop:nloop]
        elseif secondorder == 4
            fidx = find(abs(tx(:,j)-xval)<0.01);
            idx = find(abs(tx(:,j)-x)<0.01 & ty(:,j) == y & good(:,j) == gval);
            if ~isempty(idx) & ~isempty(fidx)
                idx = intersect(idx,fidx+1);
            end

% if secondorder == 5, build RC for all preceding x vals, similar
% to secondorder == 4. But now sort by  second value. So that
% if aloop is set to n, and nloop is only 1, show all responses to 
% stimulus n, with a different line for each preceding stimulus.
        elseif secondorder == 5
            fidx = find(abs(tx(:,j)-x)<0.01);
            idx = find(abs(tx(:,j)-xval)<0.01 & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
            idx = find(abs(tx(:,j)-x)<0.01 & ty(:,j) == y & good(:,j) == gval);
        else            
            idx = find(abs(tx(:,j)-x)<0.01 & ty(:,j) == y & good(:,j) == gval);
        end
        
        idx = idx(find(isextra(idx,j) == 0));
        
    end
    if ~isempty(idx)
        if ~isfield(Expt.Trials,'TrialStart')
            fprintf('%s not build with -trials, Cannot do RC\n');
            return;
        end
        if ~isempty(latency)
            trst = Expt.Trials(j).TrialStart;
            i_idx = find(([Expt.Trials(j).Start(idx)]-trst)>latency &...
                ([Expt.Trials(j).Start(idx)]-trst)<=(latency+duration));
            Expt.Trials(j).Trigger = Expt.Trials(j).Start(idx(i_idx))' -trst;
             ns = ns+length(i_idx);
             allfr = find([Expt.Trials(j).Start]-trst>latency &...
                ([Expt.Trials(j).Start]-trst)<=(latency+duration));
             all_fr(j) = length(allfr);
             fl_fr(j,1) = allfr(1);
             fl_fr(j,2) = allfr(end);
        else  % this is the code run for StimInducedCP(...'whole_trial')
            Expt.Trials(j).Trigger = Expt.Trials(j).Start(idx)' - ...
                Expt.Trials(j).TrialStart;
            ns = ns+length(idx);
            all_fr(j) = length(Expt.Trials(j).TrialStart);
            fl_fr(j,1) = 1;
            fl_fr(j,2) = length(Expt.Trials(j).Start);
        end
        tidx = [tidx j];       
    end
end
function [Expt] = AdjustSpikeTimes(Expt)
if ~isfield(Expt.Trials,'Spikes')
    for n = 1:length(Expt.Trials)
        Expt.Trials(n).Spikes = Expt.Trials(n).oSpikes;
    end
end
for n = 1:length(Expt.Trials)
    Expt.Trials(n).or_seq= mod(Expt.Trials(n).or_seq,180);
    % a total hack 
    if size(Expt.Trials(n).or_seq,1)>1
        disp('fixed or_seq field')
        Expt.Trials(n).or_seq = Expt.Trials(n).or_seq(1,:);
    end
    if length(Expt.Trials(n).Start>=1)
    Expt.Trials(n).Spikes = Expt.Trials(n).Spikes+Expt.Trials(n).TrialStart-Expt.Trials(n).Start(1);
    Expt.Trials(n).Start = (Expt.Trials(n).Start'-Expt.Trials(n).TrialStart)*10000;
    Expt.Trials(n).TrialStart = Expt.Trials(n).Start(1);
    Expt.Trials(n).Spikes = round(Expt.Trials(n).Spikes*10000);
    end
end
function [Expt, tidx,ns] = NewFindTrials(Expt, x, y, tx, ty, good, isextra, varargin)
%Bizarely, this is slower, even though the main step is done as on matrix;
secondorder =  0;
loop = 0;
extra = 0;
gval = 1;
j =1;
while j < nargin -6
    if strncmpi(varargin{j},'secondorder',4)
        j = j+1;
        loop = varargin{j}(2);
        secondorder = varargin{j}(1);
        if ismember(secondorder, [4 5])
            xval = varargin{j}(3);
        end
    elseif strncmpi(varargin{j},'xvals',4)
        j = j+1;
        xvals = varargin{j};
    elseif strncmpi(varargin{j},'RespDir',5)
        j = j+1;
        gval = varargin{j};
    elseif strncmpi(varargin{j},'extra',4)
        j = j+1;
        extra = varargin{j};
    end
    j = j+1;
end
tidx = [];
ns = 0;
[fidx, tidx] = find(tx == x & ty == y & good == gval & isextra == 0);
ns = length(tidx);
for j = 1:length(Expt.Trials)
    
    id = find(tidx == j);
    Expt.Trials(j).Trigger = Expt.Trials(j).Start(fidx(id))' - Expt.Trials(j).TrialStart;
end

%------------------
function Expt = Adjust4RCperiod(Expt);
% check version of VisStim (pre 1.0.20 we had a bug that shifted all
% stimuli by RCperiod-1 frames
newversion = 0;
if isfield(Expt,'Header') && isfield(Expt.Header,'VisStimVersion')
    VSversion = Expt.Header.VisStimVersion;
    % parse VS version
    dots = findstr(VSversion,'.');
    if length(dots) ==3
        if str2num(VSversion(1:dots(1)-1)) >1
            newversion = 1;
        elseif str2num(VSversion(dots(1)+1:dots(2)-1)) >0
            newversion = 1;
        elseif str2num(VSversion(dots(2)+1:dots(3)-1)) > 19
            newversion = 1;
        end
    end
end
if newversion
    for n=1:length(Expt.Trials)
        idx = [1:Expt.stim.vals.RCperiod:length(Expt.Trials(n).Start)];
        Expt.Trials(n).Start = Expt.Trials(n).Start(idx);
        Fields = {'phase_seq','st_seq','me_seq','sf_seq','or_seq','co_seq'};
        for n2 = 1:length(Fields)
            if isfield(Expt.Trials(n),Fields{n2}) && length(eval(['Expt.Trials(n).' Fields{n2}]))>1
                eval(['Expt.Trials(n).' Fields{n2} '= Expt.Trials(n).' Fields{n2} '(1,1:size(idx,2));'])
            end
        end
    end
else 
    % to remove the bug in the old versions we need to remove the first
    % video frame and stimulus in the sequence. (This first stimulus was only on for
    % one video frame, not for n = RCperiod number of frames.) Hence, if we
    % do not correct, all stimuli are shifted by RCperiod-1 frames (e.g.
    % latency for RCperiod = 3 is 20ms too short when we have a 100Hz
    % framerate).
    for n=1:length(Expt.Trials)
        % removing the time-stamp of the first frame
        Expt.Trials(n).Start = Expt.Trials(n).Start(2:end);
        idx = [1:Expt.stim.vals.RCperiod:length(Expt.Trials(n).Start)];
        Expt.Trials(n).Start = Expt.Trials(n).Start(idx);
        Fields = {'phase_seq','st_seq','me_seq','sf_seq','or_seq','co_seq'};
        for n2 = 1:length(Fields)
            if isfield(Expt.Trials(n),Fields{n2}) && length(eval(['Expt.Trials(n).' Fields{n2}]))>1
                % removing the first stimulus in the sequence
                eval(['Expt.Trials(n).' Fields{n2} '= Expt.Trials(n).' Fields{n2} '(1,2:size(idx,2)+1);'])
            end
        end
    end
end
        
