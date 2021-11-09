function lfps_pair_batch_visualize(lfps, analysis, splittype)

if nargin < 2; analysis = {'all'}; end
if nargin < 3; splittype = 'drug'; end

% single pair or more? ===========================
ndrug = 1;
switch splittype
    case 'drug'
        ndrug = 2;
        pairnames = {'base', 'NaCl'; 'base', '5HT'};
    case 'sc'
        pairnames = {'small sc', 'large sc'; 'small sc', 'large sc'};
    case 'pupil'
        pairnames = {'small ps', 'large ps'; 'small ps', 'large ps'};
end
if isfield(lfps, 'lfplist')
    fnames = lfps.lfplist;
    datast = lfps.LFP_prepro;
    lists = [lfps.goodunit', lfps.is5ht', lfps.animal'];
else
    fnames = '';
    datast{1} = lfps;
    ndrug = 1;
    lists = [1, lfps.is5ht, lfps.ismango];
end

close all;
disp(['5HT: ' num2str(sum(lists(:,2)==1)) ' pairs'])
disp(['NaCl: ' num2str(sum(lists(:,2)==0)) ' pairs'])

animals = {'kaki', 'mango', 'both'};
if ndrug==1 
    lists(:,2) = 0;
end
lena = length(animals);
lenses = size(lists, 1);

% exclude no data sessions
if strcmp(datast{1}.stm.param, 'rc')
    outs = zeros(lenses, 1);

    for i = 1:lenses
        if datast{i}.cond(1).sta.nspk==0 || datast{i}.cond(2).sta.nspk==0
            outs(i) = 1;
        end
    end
    fnames(outs==1) = [];
    datast(outs==1) = [];
    lists(outs==1, :) = [];
    lenses = size(lists, 1);
end

% path =======================================
% TODO: MIGHT NEED TO REMOVE THIS
if ismember(1, contains(cd, 'gpfs0'))
    mypath = '/gpfs01/nienborg/group';
elseif ismember(1, contains(cd, '\\172.25.250.112'))
    mypath = '//172.25.250.112/nienborg_group';
else
    mypath = 'E:/';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

% stimulus type ============================
stmdur = 0.45;
switch datast{1}.stm.param
    case 'rc' 
        stmidx = ones(lenses, 1);
        stmlab = {'ORxRC'};
        stmdur = 2;
    case 'or' % TODO: CHECK IF SHOULD BE REMOVED AND OTHERS BENEATH
        stmidx = zeros(lenses, 2);
        stmlab = {'anti-pref', 'pref'};
        for i = 1:lenses
            [~, sortidx] = sort(datast{i}.cond(1).spk_tu{2}(:,1));
            mini = 2;
            maxi = length(sortidx);
            try
                while isnan(datast{i}.cond(1).sta.p{sortidx(mini)})...
                        | isnan(datast{i}.cond(2).sta.p{sortidx(mini)})
                    mini = mini + 1;
                end
                while isnan(datast{i}.cond(1).sta.p{sortidx(maxi)})...
                        | isnan(datast{i}.cond(2).sta.p{sortidx(maxi)})
                    maxi = maxi - 1;
                end
                stmidx(i, :) = [sortidx(mini), sortidx(maxi)];
            catch
                disp('')
            end
        end
    case 'co'
        stmidx = zeros(lenses, 3);
        stmlab = {'0.25', '0.5', '1'};
        cos = [0.25, 0.5, 1];
        for i = 1:lenses
            for c = 1:3
                [~, idx] = min(abs(datast{i}.stm.vals - cos(c)));
                stmidx(i, c) = idx;
            end
        end
    case {'sz', 'sf'} 
        stmidx = zeros(lenses, 2);
        stmlab = {'min res stm', 'max res stm'};
        for i = 1:lenses
            [~, sortidx] = sort(datast{i}.cond(1).spk_tu{2}(:,1));
            mini = sortidx(1);
            if mini == 1
                mini = sortidx(2);
            end
            stmidx(i, :) = [mini, sortidx(end)];
        end
end

% remove no-data sessions
fnames(stmidx(:,1)==0) = [];
datast(stmidx(:,1)==0) = [];
lists(stmidx(:,1)==0, :) = [];
stmidx(stmidx(:,1)==0, :) = [];
lenses = size(lists, 1);

nrow = floor(sqrt(lenses));
ncol = ceil(lenses/nrow);
ss = size(stmidx, 2);
cols = [0 0 0; 1 0 0];
j = 0;

%% 
% firing rate as a function of trial ======================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'fr'))==1
    figure(j +1);
    for s = 1:ss
        for i = 1:lenses
            subplot(nrow, ncol, i)
            for d = 1:2
                vec = datast{i}.cond(d).mat{end}(:,1);
                trt = length(vec);
                me = mean(vec);
                vec = locdetrend(vec, 1, [20, 1]);
%                 vec = vec - mean(vec);
                vec = vec - mean(vec) + me;
%                 vec(vec < 0) = 0;
%                 vec = filter(b1, a1, vec);
                plot(1:trt, vec, '.', 'color', cols(d, :))
                hold on;
                beta = glmfit([1:trt]', vec);
                plot(1:trt, glmval(beta, 1:trt, 'identity'), '-', 'color', cols(d, :))
                hold on;
            end
            
            % format
%             set(gca, 'XTick', [1 trt])
            set(gca, 'box', 'off', 'tickdir', 'out')
            if s==1
%                 title([animals{lists(i,3)+1} ', ' pairnames{lists(i,2)+1, 2}])
                    title(fname2title(fnames{i}{2}), 'fontsize', 7)
            end
        end
    end
    j = j + 1;   
end

%% 
% LFP time-series =============================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'lfp'))==1
    disp('LFP time-series analysis ---------------------------------')
    figure(j + 1);              
    for s = 1:ss
        for i = 1:lenses
            subplot(nrow, ncol, i)
            % plot
            for d = 1:2                
                me = nanmean(datast{i}.cond(d).lfpfull{stmidx(i, s)}, 1);
                sem = nanstd(datast{i}.cond(d).lfpfull{stmidx(i, s)}, [], 1)/...
                    sqrt(size(datast{i}.cond(d).lfpfull{stmidx(i, s)}, 1));
                fill_between(datast{i}.ts, me-sem, me+sem, cols(d, :), 0.5)
                hold on;
                plot(datast{i}.ts, me, '-', 'color', cols(d, :))
                hold on;
            end
            
            % format
            xlim([datast{i}.ts(1) datast{i}.ts(end)])
            yl = get(gca, 'YLim');
            set(gca, 'XTick', [0 datast{i}.ts(end)])
            plot([0 0], yl, '--', 'color', 0.5*[1 1 1])
            set(gca, 'box', 'off', 'tickdir', 'out')
            if s==1
%                 title([animals{lists(i,3)+1} ', ' pairnames{lists(i,2)+1, 2}])
                    title(fname2title(fnames{i}{2}), 'fontsize', 7)
            end
        end
    end
    set(gcf, 'Name', 'LFP time-series', 'NumberTitle', 'off')
    j = j + 1;   
end

%%
% spike-triggered average LFP ===================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'sta'))==1
    disp('STA analysis ----------------------------------------')
    figure(j+1);
    sta_t = linspace(-datast{end}.wnd, datast{end}.wnd, ...
        size(datast{end}.cond(1).sta.mean, 2));
    xrange = [-0.1 0.1];
    for s = 1:ss
        for i = 1:lenses
            subplot(nrow, ncol, i)
            % plot
            nspk = zeros(1, 2);
            for d = 1:2
                sta_temp = datast{i}.cond(d).sta.mean(stmidx(i, s), :);
                nspk(d) = datast{i}.cond(d).sta.nspk(stmidx(i, s)); 
                plot(sta_t, sta_temp, '-', 'color', cols(d, :))
                hold on; 
            end
            
            % format
            xlim(xrange)
            yl = get(gca, 'YLim');
            set(gca, 'XTick', [xrange(1) 0 xrange(end)])
            plot([0 0], yl, '--', 'color', 0.5*[1 1 1])
            set(gca, 'box', 'off', 'tickdir', 'out')
            xl = get(gca, 'XLim');
            for d = 1:2
                text(xl(1)+0.6*(xl(2)-xl(1)), yl(1)+(0.05 + 0.1*(d-1))*(yl(2)-yl(1)), ...
                    num2str(nspk(d)), 'color', cols(d,:), 'fontsize', 6)
            end
            if s==1
%                 title([animals{lists(i,3)+1} ', ' pairnames{lists(i,2)+1, 2}])
                title(fname2title(fnames{i}{2}), 'fontsize', 7)
            end
        end
    end
    set(gcf, 'Name', 'stLFP', 'NumberTitle', 'off')
    j = j + 1;
end


%%
% Spectrogram =============================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'spectrogram'))==1
    disp('Spectrogram analysis ----------------------------------------')
    for s = 1:ss
        for i = 1:lenses
            for d = 1:2
                sz = size(datast{i}.cond(d).spectrogram.S{stmidx(i, s)});
                if sz(2) > sz(1)
                    S = 10*log10(datast{i}.cond(d).spectrogram.S{stmidx(i, s)});
                else
                    S = 10*log10(datast{i}.cond(d).spectrogram.S{stmidx(i, s)})';
                end 
                spe_t = datast{end}.cond(1).spectrogram.t{1};
                lent = length(spe_t);
                spe_t = linspace(-0.1, stmdur, lent);
                Sabs = nanmean(S(:, spe_t > datast{end}.window{end}(1)), 2);
                subplot(nrow, ncol, i)
                plot(datast{i}.cond(d).spectrogram.f{stmidx(i, s)}, Sabs, '-', 'color', cols(d, :))
                hold on;
            end
        
            % format
            xlim([datast{i}.cond(d).spectrogram.f{stmidx(i, s)}(1), ...
                datast{i}.cond(d).spectrogram.f{stmidx(i, s)}(end)])
            set(gca, 'box', 'off', 'tickdir', 'out')
            set(gca, 'XScale', 'log')
            if s==1
    %                 title([animals{lists(i,3)+1} ', ' pairnames{lists(i,2)+1, 2}])
                    title(fname2title(fnames{i}{2}), 'fontsize', 7)
            end
        end
    end
    set(gcf, 'Name', 'PSD', 'NumberTitle', 'off')
    j = j + 1; 
end

%%
% coherence ==============================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'coherence'))==1
    disp('coherence analysis --------------------------------------')
    for s = 1:ss
        for i = 1:lenses
            for d = 1:2
                subplot(nrow, ncol, i)
                plot(datast{i}.cond(d).coherence.f{stmidx(i, s)}, ...
                    datast{i}.cond(d).coherence.C{stmidx(i, s)}, '-', 'color', cols(d, :))
                hold on;
            end
        
            % format
%             xlim([datast{i}.cond(d).coherence.f{stmidx(i, s)}(1), ...
%                 datast{i}.cond(d).coherence.f{stmidx(i, s)}(end)])
            xlim([3 48])
            set(gca, 'box', 'off', 'tickdir', 'out')
            if s==1
    %                 title([animals{lists(i,3)+1} ', ' pairnames{lists(i,2)+1, 2}])
                    title(fname2title(fnames{i}{2}), 'fontsize', 7)
            end
        end
    end
    set(gcf, 'Name', 'coherence', 'NumberTitle', 'off')
    j = j + 1;  
end

%%
% Tuning =================================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'tuning'))==1
    disp('tuning analysis ---------------------------')
    
%     addpath(genpath([mypath '/Katsuhisa/code/integrated/cbrewer']))
%     CT = cbrewer('div', 'RdBu', 13);
%     CT = flipud(CT);

    % data extraction
    if strcmp(datast{1}.stm.param, 'rc')        
        % load Corinna's analysis =================
        load([mypath '/serotonin_project/LFP_project/Data/others/smlinfo.mat'])
        
        indnames = {'mean', 'selectivity', 'snr2'};
        tunames = {'spikes'};
        
        % compute tuning parameters
        datast = encoding_tuning4rc(datast, smlinfo);
    else
        indnames = {'mean', 'reliability', 'selectivity', 'snr2', 'discriminability', 'entropy', 'c-entropy', 'MI'};
        tunames = {'spikes', 'theta', 'alpha', 'beta', 'slow_gamma', 'fast_gamma'};
    end
    lent = length(tunames);
    leni = length(indnames);
    tudata = cell(2, lent);
    for i = 1:lenses % sessions
        for d = 1:2 % base or drug
            for t = 1:lent % activity type
                for l = 1:leni % index
                    try
                        if strcmp(indnames{l}, 'selectivity') && ismember(1, contains({'or', 'rc'}, datast{1}.stm.param))
                            tudata{d, t}(i, l) = datast{i}.cond(d).tuning.(tunames{t}).unique.circularvariance;
                        elseif strcmp(indnames{l}, 'snr2') || strcmp(indnames{l}, 'mean')
                            tudata{d, t}(i, l) = mean(datast{i}.cond(d).tuning.(tunames{t}).(indnames{l})); 
                        elseif strcmp(indnames{l}, 'entropy')           
                            tudata{d, t}(i, l) = datast{i}.cond(d).tuning.(tunames{t}).metabcost(1);
                        elseif strcmp(indnames{l}, 'c-entropy')           
                            tudata{d, t}(i, l) = datast{i}.cond(d).tuning.(tunames{t}).metabcost(2);
                        elseif strcmp(indnames{l}, 'MI')           
                            tudata{d, t}(i, l) = datast{i}.cond(d).tuning.(tunames{t}).metabcost(3);
                        else
                            tudata{d, t}(i, l) = datast{i}.cond(d).tuning.(tunames{t}).(indnames{l});
                        end
                    catch
                        tudata{d, t}(i, l) = nan;
                    end
                end
            end
        end
    end       

    % tuning index ==================================
    for l = 1:leni % index
        figure(j + l + 1);
        for k = 1:ndrug % is5ht
            for t = 1:lent % activity type
                subplot(ndrug, lent, t + lent*(k-1))
                try
                    anim = 1*lists(lists(:,2)==k-1, 3);
                    X = tudata{1, t}(lists(:,2)==k-1, l);
                    Y = tudata{2, t}(lists(:,2)==k-1, l);
                    nans = isnan(X) | isnan(Y);
                    X(nans) = []; Y(nans) = []; anim(nans) = [];
                    unity_scatter(X, Y, anim);
                catch
                    continue
                end

                % format
                set(gca, 'box', 'off', 'tickdir', 'out')
                if k==1
                    title(tunames{t})
                end
                if t==1
                    if k==1
                        ylabel({'NaCl', '# pairs'})
                    else
                        xlabel(indnames{l})
                        ylabel({'5HT', '# pairs'})
                    end
                end
            end
        end
        set(gcf, 'Name', indnames{l}, 'NumberTitle', 'off')
    end
    j = j + l + 1;

%     % confusion matrix (spikes vs LFP) ===============================
%     if ~strcmp(datast{1}.stm.param, 'rc')
%         cmr = cell(2,2);
%         cmp = cell(2,2);
%         for k = 1:ndrug % is5ht
%             for d = 1:2 % base or drug
%                 cmr{k, d} = zeros(leni, lent);
%                 cmp{k, d} = zeros(leni, lent);
%                 for l = 1:leni % index
%                     for t = 1:lent
%                         [rr, pp] = corrcoef(tudata{d, 1}(lists(:,2)==k-1, l),...
%                             tudata{d, t}(lists(:,2)==k-1, l));
%                         cmr{k, d}(l, t) = rr(1, 2);
%                         cmp{k, d}(l, t) = pp(1, 2);
%                     end
%                 end            
%             end
%         end

%         % visualize
%         figure(j + l + 2);
%         clim = [0 0];
%         for k = 1:ndrug
%             % base & drug
%             for d = 1:2
%                 subplot(ndrug, 3, d + 3*(k-1))
%                 imagesc(1:lent, 1:leni, cmr{k, d})
%                 caxis([0 1])
%                 colormap(CT)
%                 colorbar('northoutside')
%                 set(gca, 'YTick', 1:leni, 'YTickLabel', indnames)
%                 set(gca, 'XTick', 1:lent, 'XTickLabel', tunames)
%                 xtickangle(45)
%                 ytickangle(45)
%                 set(gca, 'box', 'off', 'tickdir', 'out')
%                 if k==1 && d==1
%                     title('baseline')
%                 elseif k==1 && d==2
%                     title('drug')
%                 end
%                 if d==1
%                     if collapse==1
%                         ylabel({drugnames, 'spikes'})
%                     else
%                         ylabel({drugnames{k}, 'spikes'})
%                     end
%                 end
%             end
% 
%             % difference
%             subplot(ndrug, 3, 3 + 3*(k-1))
%             imagesc(1:lent, 1:leni, cmr{k, 1} - cmr{k, 2})
%             ctemp = caxis;
%             if ctemp(1) < clim(1)
%                 clim(1) = ctemp(1);
%             end
%             if ctemp(2) > clim(2)
%                 clim(2) = ctemp(2);
%             end
%             set(gca, 'box', 'off', 'tickdir', 'out')
%         end
%         clim = max(abs(clim))*[-1 1];
%         for k = 1:ndrug
%             subplot(ndrug, 3, 3 + 3*(k-1))
%             caxis(clim)
%             colormap(CT)
%             set(gca, 'YTick', 1:leni, 'YTickLabel', indnames)
%             set(gca, 'XTick', 1:lent, 'XTickLabel', tunames)
%             xtickangle(45)
%             ytickangle(45)
%             colorbar('northoutside')
%             if k==1
%                 title('\Delta')
%             end
%         end
%         set(gcf, 'Name', 'spike vs LFP in tuning', 'NumberTitle', 'off')
%         j = j + l + 3;
%     end
        
    % coherence vs tuning parameter    
    if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'coherence'))==1
        for l = 1:leni % index
            figure(j + l);
            for a = 1:lena
                for k = 1:ndrug % is5ht
                    if a < 3
                        cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                    else
                        cond = lists(:,2)==k-1;
                    end
                    freq =52; % Hard code
                    coh_base = zeros(sum(cond), length(freq));
                    coh_drug = coh_base;
                    for s = 1:ss
                        %temp = squeeze(coh_base{s}(cond, :, 1));
                        nanrows = any(isnan(temp), 2);
                        med = nanmedian(temp, 1);
                        temp(nanrows, :) = repmat(med, sum(nanrows), 1);
                        coh_base = coh_base + temp;
                        temp = squeeze(coh{s}(cond, :, 2));
                        nanrows = any(isnan(temp), 2);
                        med = nanmedian(temp, 1);
                        temp(nanrows, :) = repmat(med, sum(nanrows), 1);
                        coh_drug = coh_drug + temp;
                    end
                    coh_base = coh_base/ss;
                    coh_drug = coh_drug/ss;
                    for b = 1:lenb % activity type
                        % scatter & regression line
                        if b==1
                            frange = bandrange{b}(1) <= freq & ...
                                bandrange{b+1}(2) >= freq;
                        elseif b > 1 && b < lenb
                            frange = bandrange{b+1}(1) <= freq & ...
                                bandrange{b+1}(2) >= freq;
                        else
                            frange = bandrange{b+1}(1) <= freq & ...
                                80 >= freq;
                        end
                        subplot(lena*ndrug, lenb, b + lenb*(k-1) + ndrug*lenb*(a-1))   
                        X = [mean(coh_base(:, frange), 2); mean(coh_drug(:, frange), 2)];
                        Y = [tudata{1, 1}(cond, l); tudata{2, 1}(cond, l)];% spikes
                        ls_scatter(X, Y, [zeros(sum(cond), 1); ones(sum(cond), 1)], 1);                             

                        % format
                        if k==1 && a==1
                            title(bandnames{b})
                        end
                        if b==1
                            if k==1
                                ylabel({animals{a}, pairnames{k, 2}, indnames{l}})
                            else
                                if a==lena
                                    xlabel('spike-LFP coherence')
                                end
                                ylabel({animals{a}, pairnames{k, 2}, indnames{l}})
                            end
                        end
                    end
                end
            end
            set(gcf, 'Name', ['spike-LFP coherence vs ' indnames{l} ' of spikes'], 'NumberTitle', 'off')
        end

    end
    j = j + l;
end



%%
% subfunction
function datast = encoding_tuning4rc(datast, smlinfo)
% compute encoding indices from tuning curve
if length(datast) > 50
    sesidx = find(smlinfo.paramat(:,4)==0);
else
    sesidx = find(smlinfo.paramat(:,3)==1 & smlinfo.paramat(:,4)==0);
end
lenses = length(datast);
dlabs = {'', '_drug'};
for i = 1:lenses
    for d = 1:2
        % tuning values
        or = smlinfo.(['fitparam' dlabs{d}]){sesidx(i)}.val.or;
        mn = smlinfo.(['fitparam' dlabs{d}]){sesidx(i)}.val.mn;
        sem = smlinfo.(['fitparam' dlabs{d}]){sesidx(i)}.val.sem;
        
        % mean
        datast{i}.cond(d).tuning.spikes.mean = mn;
        
        % SNR2
        snr2 = (mn./sem).^2;
        snr2(isnan(snr2)) = 0;
        snr2(isinf(snr2)) = 100;
        datast{i}.cond(d).tuning.spikes.snr2 = snr2;
        
        % circular variance --- Ringach et al. (2002)
        or = or * pi /180;
        or = mod(or, 2*pi);
        if sum(mn < 0) > 0
            mn = mn + abs(min(mn));
        end
        % compute weighted sum of cos and sin of angles
        r = sum(mn.*exp(1i*or));
        % obtain length 
        r = abs(r)./sum(mn);
        datast{i}.cond(d).tuning.spikes.unique.circularvariance = 1 - r;
    end
end

function rad = maxcoh_rad(coh, rad)
sz = size(rad);
rad_in = rad;
rad = zeros(sz(1), 1);
[~, maxi] = max(coh);
if maxi==1
    idx = [1 2 3];
elseif maxi==sz(2)
    idx = [sz(2)-2 sz(2)-1 sz(2)];
else
    idx = [maxi-1 maxi maxi+1];
end
for i = 1:sz(1)
    rad(i) = circ_mean(rad_in(i, idx)');
end

function tlab = fname2title(fname)
% dotpos = strfind(fname, '.');
% tlab = [fname(1:3) fname(dotpos(end-1)+1:dotpos(end)-1)];
sti = strfind(fname, '_');
dots = strfind(fname, '.');
tlab = [fname(1:sti(2)-1) ':' fname(sti(end)+1:dots(end)-1)];