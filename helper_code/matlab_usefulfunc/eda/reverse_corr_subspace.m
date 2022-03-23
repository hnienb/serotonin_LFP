function rcsub = reverse_corr_subspace(stmMat, actMat, wnd, stm_samprate, plot_flag)
% generic function to perform a reverse correlation subspace analysis
% INPUT: 
% stmMat ... row: trials, column: stimulus values (e.g. hdx, or...)
% actMat ... row: trials, column: numeric vector (e.g. spikes, LFP...)
% wnd ... analysis window (ms)
% stm_samprate ... sampling rate of stimulus presentation  
% plot_flag ... 0, no plot; 1, plot
%
% OUTPUT:
% rcsub ... activity elicited by stimulus
%
% written by Katsuhisa (06.04.18)
% ++++++++++++++++++++++++++++++++++++++++++++++++++
 
if nargin<2; error('Provide stmMat & actMat!'); end
if nargin<3; wnd = 150; end
if nargin<4; stm_samprate = 100; end
if nargin<5; plot_flag = 0; end

% experiment setting -------------------------------
% the number of trials
ntr = size(stmMat, 1);

% stimulus duration
nframes = size(stmMat, 2);
stmdur = 1000*nframes/stm_samprate;
ndata = size(actMat ,2);

% equalize sampling rate of stmMat & actMat
% (note that I assume the sampling rate for actMat >= stmMat)
% actMatBin = nan(ntr, nframes);
% nperbin = floor(size(actMat, 2)/nframes);
% begin = 1;
% for n = 1:nframes
%     actMatBin(:,n) = nanmean(actMat(:, begin:begin+nperbin-1), 2);
%     begin = begin + nperbin;
% end
upstmMat = zeros(ntr, ndata);
for n = 1:ntr
    upstmMat(n,:) = interp1(1:nframes, stmMat(n,:), ...
                linspace(1, nframes, ndata), 'nearest');
end

% reverse correlation --------------------
nframeperwnd = floor(ndata*wnd/stmdur);

% initialize output structure 
unique_stm = unique(round(1000*stmMat(:))/1000);
nstm = length(unique_stm);
rcsub = struct('stm', [], 'res', [], 'res_sd', [], 'latency', []);
for n = 1:nstm
    rcsub.stm(n).val = unique_stm(n);
    rcsub.stm(n).idx = [];
end
% disp('accumulating responses...')
c = 1;
for i = 1:ntr
    stpos = 1;
    while stpos + nframeperwnd - 1 <= ndata
        stmidx = unique_stm == upstmMat(i, stpos);
        rcsub.stm(stmidx).idx = [rcsub.stm(stmidx).idx; c];
        c = c + 1;
        rcsub.res = [rcsub.res; actMat(i, stpos:stpos+nframeperwnd-1)];
        stpos = stpos + 1;
    end
end

% grand mean of responses
grandmean = 1000*mean(rcsub.res(:));

% response elicited by stimulus
m = nan(nstm, wnd);
for n = 1:nstm
    rcsub.stm(n).mean = 1000*mean(rcsub.res(rcsub.stm(n).idx, :), 1);
    rcsub.stm(n).sem = 1000*std(rcsub.res(rcsub.stm(n).idx, :), [], 1)...
        /sqrt(length(rcsub.stm(n).idx));    
    % upsample so that 1 frame becomes 1 ms 
     rcsub.stm(n).mean = interp1(1:nframeperwnd, rcsub.stm(n).mean, ...
                linspace(1, nframeperwnd, wnd), 'linear');
     rcsub.stm(n).sem = interp1(1:nframeperwnd, rcsub.stm(n).sem, ...
                linspace(1, nframeperwnd, wnd), 'linear');    
    % smoothing (4ms boxcar convolution)
    rcsub.stm(n).mean = boxcar_smooth(rcsub.stm(n).mean, 4);
    rcsub.stm(n).sem = boxcar_smooth(rcsub.stm(n).sem, 4);
    m(n, :) = rcsub.stm(n).mean;
    % metric
    [~, rcsub.stm(n).peak_t] = max(abs(rcsub.stm(n).mean));  % peak time (ms)
    rcsub.stm(n).peak = rcsub.stm(n).mean(rcsub.stm(n).peak_t); % peak value
    rcsub.stm(n).totalcounts = sum(rcsub.stm(n).mean)/wnd; % counts per frame
end

% latency estimate
sdvec  = std(m, [], 1);
rcsub.res_sd = sdvec - mean(sdvec(1:20)); % baseline correction
rcsub.latency = find(rcsub.res_sd > 0.5*max(rcsub.res_sd), 1, 'first'); % latency

% visualization -----------------------------
if plot_flag==1
    close all;
    col = jet(nstm);
    h = figure;
    % mean in each stimulus
    subplot(1,2,1)
    plot([0 wnd+1], grandmean.*[1 1], ':k')
    hold on;
    xlim([0 wnd+1])
    ylabel({'response', '(mean)'})
    xlabel('time (ms)')
    set(gca, 'XTick', [0 wnd+1], 'XTickLabel', [0 wnd])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')    
    pl = nan(1, nstm);
    leg = cell(1, nstm);
    for n = 1:nstm
        subplot(1,2,1)
        hold on;
        pl(n) = plot(1:wnd, rcsub.stm(n).mean, '-', 'color', col(n,:));
        leg{n} = ['stm:' num2str(unique_stm(n))];
    end
    legend(pl, leg, 'location', 'northeast'); legend('boxoff')
    % sd
    subplot(1,2,2)
    plot(1:wnd, rcsub.res_sd, '-r') 
    hold on;
    yy = get(gca, 'YLim');
    plot(rcsub.latency*[1 1], yy, ':k')
    xlim([0 wnd+1])
    ylabel('(SD)')
    xlabel('time (ms)')
    set(gca, 'XTick', [0 wnd+1], 'XTickLabel', [0 wnd])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    set(h, 'Name', ['ntr = ' num2str(ntr)], 'NumberTitle', 'off')
end