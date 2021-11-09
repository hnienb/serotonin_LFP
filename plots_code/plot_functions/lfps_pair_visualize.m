function lfps_pair_visualize(lfps, analysis, splittype)

if nargin < 2; analysis = {'all'}; end
if nargin < 3; splittype = 'drug'; end

% single pair or more? ===========================
ndrug = 1;
switch splittype
    case 'drug'
        ndrug = 2;
        pairnames = {'base', 'NaCl'; 'base', '5HT'};
    case 'sc'
        pairnames = {'large sc', 'small sc'; 'large sc', 'small sc'};
    case 'pupil'
        pairnames = {'small ps', 'large ps'; 'small ps', 'large ps'};
end
if isfield(lfps, 'lfplist')
    datast = lfps.LFP_prepro;
    lists = [lfps.goodunit', lfps.is5ht', lfps.animal'];
else
    datast{1} = lfps;
    ndrug = 1;
    lists = [1, lfps.is5ht, lfps.ismango];
end

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
    datast(outs==1) = [];
    lists(outs==1, :) = [];
    lenses = size(lists, 1);
end
disp(['analyzed ses:' num2str(lenses)])
% path =======================================

% TODO: FIX THIS OR REMOVE IT
if ismember(1, contains(cd, 'gpfs0'))
    mypath = '/gpfs01/nienborg/group';
elseif ismember(1, contains(cd, '\\172.25.250.112'))
    mypath = '//172.25.250.112/nienborg_group';
else
    mypath = 'Z:';
end


addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/cbrewer']))

% stimulus type ============================
stmdur = 0.45;
switch datast{1}.stm.param
    case 'rc'     
        stmidx = ones(lenses, 1);
        stmlab = {'ORxRC'};
        stmdur = 2;
    case 'or' % TODO: Check to see if I need to keep this (and the rest underneath)
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
        if length(datast{end}.stm.vals) > 1
            stmidx = zeros(lenses, 3);
            stmlab = {'0.25', '0.5', '1'};
            cos = [0.25, 0.5, 1];
            for i = 1:lenses
                for c = 1:3
                    [~, idx] = min(abs(datast{i}.stm.vals - cos(c)));
                    stmidx(i, c) = idx;
                end
            end
        else
            stmidx = ones(lenses, 1);
            stmlab = {'ORxRC (CO)'};
            stmdur = 2;
        end
    case {'sz', 'sf'} % sf = spatial frequency, size stim parameters
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
datast(stmidx(:,1)==0) = [];
lists(stmidx(:,1)==0, :) = [];
stmidx(stmidx(:,1)==0, :) = [];
lenses = size(lists, 1);

ss = size(stmidx, 2);
cols = [0 0 0; 1 0 0];
j = 0;

%%
% check thinning ===========================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'thinning'))==1
    figure(j+1);
    if isempty(datast{1}.cond(1).spk_tu{3})
        sidx = [2, 1];
        titles = {'stimulus duration', 'analysis window'};
    else
        sidx = [2, 1, 3];
        titles = {'stimulus duration', 'analysis window', 'thinned'};
    end
    nsidx = length(sidx);
    for s = 1:nsidx
        me = zeros(lenses, 2);
        % data extraction
        for i = 1:lenses
            me(i, 1) = mean(datast{i}.cond(1).spk_tu{sidx(s)}(:,1));
            me(i, 2) = mean(datast{i}.cond(2).spk_tu{sidx(s)}(:,1));
        end
        for d = 1:ndrug
            % plot
            subplot(ndrug, nsidx, s + nsidx*(d - 1))
            if lenses==1
                mex = me(1, 1);
                mey = me(1, 2);
                anim = 1;
            else
                mex = me(lists(:,2)==d-1, 1);
                mey = me(lists(:,2)==d-1, 2);
                anim = 1*lists(:,3);
                anim = anim(lists(:,2)==d-1);
            end
            unity_scatter(mex, mey, anim)

            % labels
            if d==1
                title(titles{sidx(s)})
            end
            if s==1
                if d==ndrug
                    xlabel({pairnames{d, 1}, 'firing rate'})
                end
                ylabel({pairnames{d, 2}, 'firing rate'})
            end
        end    
    end
    set(gcf, 'Name', 'spiking activity', 'NumberTitle', 'off')
    j = j + 1;
end

%%
% pupil  ===========================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'pupil'))==1
    figure(j+1);
    nf = 200;
    for s = 1:ss
        ps = zeros(lenses, nf, 2);
        dp = zeros(lenses, nf, 2);
        % data extraction
        for i = 1:lenses
            for d = 1:2
                % both eyes
                ps_temp = (datast{i}.cond(d).eye{stmidx(i, s), 5} + ...
                    datast{i}.cond(d).eye{stmidx(i, s), 6})/2;
                dp_temp = (datast{i}.cond(d).eye{stmidx(i, s), 7} + ...
                    datast{i}.cond(d).eye{stmidx(i, s), 8})/2;
                % zscore
                ps_temp = (ps_temp - mean(ps_temp(:)))/std(ps_temp(:));
                dp_temp = (dp_temp - mean(dp_temp(:)))/std(dp_temp(:));
                % mean
                ps(i, :, d) =  adjust_veclen(mean(ps_temp, 1), nf);
                dp(i, :, d) = adjust_veclen(mean(dp_temp, 1), nf);
            end
        end
        % visualize ===========================
        for a = 1:lena
            for k = 1:ndrug
                % time-series ==========================
                if a < 3
                    cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                else
                    cond = lists(:,2)==k-1;
                end
                
                % plot
                for d = 1:2
                    subplot(lena*ndrug, 2*ss, s + 2*ss*(k-1) + ndrug*2*ss*(a-1))
                    % ps
                    me = squeeze(nanmean(ps(cond, :, d), 1));
                    sem = squeeze(nanstd(ps(cond, : ,d), [], 1))/sqrt(sum(cond));
                    fill_between(1:nf, me-sem, me+sem, cols(d, :), 0.4)
                    hold on;
                    plot(1:nf, me, '-', 'color', cols(d, :))
                    hold on;
                    
                    subplot(lena*ndrug, 2*ss, s+1 + 2*ss*(k-1) + ndrug*2*ss*(a-1))
                    % dp
                    me = squeeze(nanmean(dp(cond, :, d), 1));
                    sem = squeeze(nanstd(dp(cond, : ,d), [], 1))/sqrt(sum(cond));
                    fill_between(1:nf, me-sem, me+sem, cols(d, :), 0.4)
                    hold on;
                    plot(1:nf, me, '-', 'color', cols(d, :))
                    hold on;
                end        
                
                % format
                xlim([0 nf])
                yl = get(gca, 'YLim');
                set(gca, 'XTick', [0 nf])
                plot([0 0], yl, '--', 'color', 0.5*[1 1 1])
                if k==1 && a==1
                    title(stmlab{s})            
                end
                if s==1
                    ylabel({animals{a}, pairnames{k, 2}, 'pupil'})
                    if k == ndrug && a==lena
                        xlabel('time after stimulus onset (sec)')  
                    end
                end
            end
        end
    end
    set(gcf, 'Name', 'pupil time-series', 'NumberTitle', 'off')
    j = j + 1;   
end

%% 
% LFP time-series =============================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'lfp'))==1
    disp('LFP time-series analysis ---------------------------------')
    ts = linspace(-0.1, stmdur, 1000*(stmdur+0.1));
    lent = length(ts);
    figure(j + 1);              
    for s = 1:ss
        % data extraction
        lfp = zeros(lenses, lent, 2);
        for i = 1:lenses
            for d = 1:2
                lfp_me = nanmean(datast{i}.cond(d).lfpfull{stmidx(i, s)}, 1);
                lfp(i, :, d) = interp1(datast{i}.ts, lfp_me, ts);
            end
        end

        % visualize ===========================
        for a = 1:lena
            for k = 1:ndrug
                % time-series ==========================
                if a < 3
                    cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                else
                    cond = lists(:,2)==k-1;
                end
                
                % plot
                subplot(lena*ndrug, ss, s + ss*(k-1) + ndrug*ss*(a-1))
                for d = 1:2
                    me = squeeze(nanmean(lfp(cond, :, d), 1));
                    sem = squeeze(nanstd(lfp(cond, : ,d), [], 1))/sqrt(sum(cond));
                    fill_between(ts, me-sem, me+sem, cols(d, :), 0.4)
                    hold on;
                    plot(ts, me, '-', 'color', cols(d, :))
                    hold on;
                end        
                
                % format
                xlim([ts(1) ts(end)])
                yl = get(gca, 'YLim');
                set(gca, 'XTick', [0 stmdur])
                plot([0 0], yl, '--', 'color', 0.5*[1 1 1])
                if k==1 && a==1
                    title(stmlab{s})            
                end
                if s==1
                    ylabel({animals{a}, pairnames{k, 2}, 'LFP (uV)'})
                    if k == ndrug && a==lena
                        xlabel('time after stimulus onset (sec)')  
                    end
                end
            end
        end
    end
    set(gcf, 'Name', 'LFP time-series', 'NumberTitle', 'off')
    j = j + 1;   
end


%%
% GLM analysis ================================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'glm'))==1
    if strcmp(splittype, 'drug') && strcmp(datast{1}.stm.param, 'rc')
        disp('GLM analysis -------------------------')
        paircol = cbrewer('qual', 'Paired', 12);
        varnames = {'ps', 'dps', 'drug', 'drug x ps', 'drug x dps', ...
            'lfp res', 'csd', 'low-freq', 'gamma', 'su', 'mu'};
%         mdlnames = {'su', 'mu', 'low-freq'};
%         mdly = [10, 11, 8];
        mdlnames = {'low-freq'};
        mdly = 8;
        mdlout = {[6:11]};
        lenm = length(mdly);
        lenp = zeros(1, lenm);
        cv = 3;
        cc = cell(1, lenm);
        w = cell(1, lenm);
%         lam = 0.084;
        for i = 1:lenses
            % x and Y
            mat0 = datast{i}.cond(1).mat{stmidx(i, end)};
            mat2 = datast{i}.cond(2).mat{stmidx(i, end)};
            ntr0 = size(mat0, 1);
            ntr2 = size(mat2, 1);
            X = [...
                [mat0(:, 3); mat2(:, 3)], ... % pupil size
                [mat0(:, 4); mat2(:, 4)], ... % pupil size derivaticc
                [zeros(ntr0, 1); ones(ntr2, 1)], ... % base or drug
                [zeros(ntr0, 1); ones(ntr2, 1)].*[mat0(: , 3); mat2(:, 3)], ... % interaction
                [zeros(ntr0, 1); ones(ntr2, 1)].*[mat0(: , 4); mat2(:, 4)], ... % interaction
                [mat0(:, 5); mat2(:, 5)], ... % LFP res
                [mat0(:, 6); mat2(:, 6)], ... % ddLFP
                [mat0(:, 7); mat2(:, 7)], ... % low-freq
                [mat0(:, 8); mat2(:, 8)], ... % gamma
                [mat0(:, 1); mat2(:,1)]/stmdur, ... % spike counts --> firing rate (su)
                [mat0(:, 2); mat2(:,2)]/stmdur, ... % spike counts --> firing rate (mu)
                ]; 

            % k-fold cross-validation
            rng(19891220);
            cvidx1 = crossvalind('Kfold', ntr0, cv);
            cvidx2 = crossvalind('Kfold', ntr2, cv);
            cvidx = [cvidx1; cvidx2];
            
            % fit GLM stepwise
%             L = [];
            for m = 1:lenm
                y = 10*log10(X(:, mdly(m)));
%                 if m < 3
%                     % firing rate model
%                     y = log(1+y);
%                 end
                predictors = zscore(X(:, ~ismember(1:size(X, 2), mdlout{m})));
                lenp(m) = size(predictors, 2);
                r_temp = nan(cv, lenp(m));
                w_temp = nan(cv, lenp(m));
                for v = 1:cv
                    for k = 1:lenp(m)
                        % model prediction (stepwise)
%                         [B, FitInfo] = lassoglm(predictors(cvidx==v, 1:k), y(cvidx==v), 'normal', 'lambda', lam);
%                         beta = [FitInfo.Intercept; B];                        
                        beta = glmfit(predictors(cvidx~=v, 1:k), y(cvidx~=v), 'normal');
                        ypred = glmval(beta, predictors(cvidx==v, 1:k), 'identity');

                        % correlation coefficient
                        rr = corrcoef(y(cvidx==v), ypred);
                        r_temp(v, k) = rr(1,2);
%                         r_temp(v, k) = varexp(y(cvidx==vec(~ismember(vec, v))), ypred);                        
                    end
%                     L = [L, FitInfo.LambdaMinDeviance];
                    % weight
                    w_temp(v, :) = beta(2:end);              
                end 
                [cvscore, cvi] = max(r_temp(:, end));
                cc{m}(i, :) = r_temp(cvi, :);      
                w{m}(i, :) = w_temp(cvi, :);
            end    
        end
%         median(L)
        % visualize
        c = 1;
        vars = cell(1, m);
        for m = 1:lenm
            for a = 1:lena
                for k = 1:ndrug
                    % bar plots ================================================
                    figure(j+c);
                    if a < 3
                        cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                    else
                        cond = lists(:,2)==k-1;
                    end
                   % individual correlation coefficients
                   subplot(lena*ndrug, 2*ss, 1 + 2*ss*(k-1) + ndrug*2*ss*(a-1))
%                    me = nanmean(squeeze(w(cond, :, m)), 1);
%                    sem = nanstd(squeeze(w(cond, :, m)), [], 1)/sqrt(sum(cond));

                   plot([0 lenp(m)], [0 0], ':k')
                   for p = 1:lenp(m)
                       pval = signrank(w{m}(cond, p));
                       if pval < 0.05/lenp(m)
                            pcol = [1 0 0];
                       else
                           pcol = [0 0 0];
                       end
%                        bar(p, me(p), 'FaceColor', barcol, 'EdgeColor', barcol)
%                        hold on;
%                        errorbar(p, me(p), sem(p), 'color', barcol, 'capsize', 0)
%                        hold on;
                        y = w{m}(cond, p);
                        x = p*ones(size(y));
                        hold on;
                        scatter(x, y, 20, 'o', 'markerfacecolor', pcol, 'markeredgecolor', pcol, ...
                            'markerfacealpha', 0.4, 'markeredgealpha', 0.1)
                        hold on;
                        plot([p-0.25 p+0.25], mean(y)*[1 1], '-', 'color', [0 0.5 0], 'linewidth', 2)
                   end
%                    hold on;
                   vars{m} = varnames(~ismember(1:lenp(m), mdlout{m}));
%                    violin(squeeze(w(cond, :, m)),'facecolor', pcol, 'edgecolor', [], 'mc', [], 'medc', 'g-')
%                    legend('off')
                   if a == 1 && k == 1
                       title(mdlnames{m})
                   end
                    ylabel({animals{a}, pairnames{k, 2}})
                    set(gca, 'XTick', 1:lenp(m), 'XTickLabel', cell(1, lenp(m)))
                    set(gca, 'box', 'off', 'tickdir', 'out')
                   if k == ndrug
                       if a == lena
                           set(gca, 'XTick', 1:lenp(m), 'XTickLabel', vars{m})
                       end
                       xtickangle(45)
                   end

                   % variance explained
                   subplot(lena*ndrug, 2*ss, 2 + 2*ss*(k-1) + ndrug*2*ss*(a-1))
                   waterfallchart4mat(cc{m}(cond, :))
                   me = nanmean(cc{m}(cond, :), 1);
%                    sem = nanstd(cc{m}(cond, :), [], 1)/sqrt(sum(cond));
                   for p = 1:lenp(m)
%                        barcol = paircol(3, :);
                       if p > 1
                           pval = signrank(cc{m}(cond, p-1), cc{m}(cond, p));
                           if pval < 0.05/lenp(m)
                               text(p, me(p)*1.1, '*')
%                                 barcol = paircol(4, :);
                           end
%                            if me(p-1) > me(p) 
%                                bar(p, me(p-1), 'FaceColor', barcol, 'EdgeColor', 'k')
%                                hold on;
%                                bar(p, me(p), 'FaceColor', 'w', 'EdgeColor', 'w')
%                            else
%                                bar(p, me(p), 'FaceColor', barcol, 'EdgeColor', 'k')
%                                hold on;
%                                bar(p, me(p-1), 'FaceColor', 'w', 'EdgeColor', 'w')
%                            end
%                            hold on;
%                            plot([p-1.4 p+0.4], me(p-1)*[1 1], '-k', 'linewidth', 0.25)
                       else
                           pval = signrank(cc{m}(cond, p));
                           if pval < 0.05/lenp(m)
                               text(p, me(p)*1.1, '*')
%                                 barcol = paircol(4, :);
                           end
%                            bar(p, me(p), 'FaceColor', barcol, 'EdgeColor', 'k')
                       end
                       disp([animals{a} ': ' vars{m}{p} ': p = ' num2str(pval*lenp(m))]) 
%                        hold on;
%                        errorbar(p, me(p), sem(p), 'color', barcol, 'capsize', 0)
%                        hold on;
                   end
                   if a == 1 && k == 1
                       title('correlation coefficient')
                   end
                   set(gca, 'XTick', 1:lenp(m), 'XTickLabel', [])
                   yy = get(gca, 'YLim');
                  yy(2) = 0.4; 
                  ylim([floor(0.7*me(1)) yy(2)])
                   set(gca, 'box', 'off', 'tickdir', 'out')
                   if k == ndrug
                       if a == lena
                           set(gca, 'XTick', 1:lenp(m), 'XTickLabel', vars{m})
                       end
                       xtickangle(45)
                   end                
                end
                set(gcf, 'Name', ['GLM: ' mdlnames{m}], 'NumberTitle', 'off') 

%                 % histogram of weights =============================================
%                figure(j+c+1);
%                if a < 3
%                   cond0 = lists(:,2)==0 & lists(:,3)==a-1;
%                   cond2 = lists(:,2)==1 & lists(:,3)==a-1;
%                else
%                   cond0 = lists(:,2)==0;
%                   cond2 = lists(:,2)==1;
%                end
%                for p = 1:lenp-1
%                    subplot(lena, lenp-1, p + (lenp-1)*(a-1))
%                    [~, pval] = ttest2(squeeze(w(cond0, p, m)), squeeze(w(cond2, p, m)));
%                    histogram(squeeze(w(cond0, p, m)), 'FaceColor', 'k')
%                    hold on;
%                    histogram(squeeze(w(cond2, p, m)), 'FaceColor', 'r')
%                    hold on;
%                    if a == 1
%                        title(varnames{2+p})
%                    end        
%                    if p==1
%                        ylabel(animals{a})
%                    end
%                    set(gca, 'box', 'off', 'tickdir', 'out')
%                    xx = get(gca, 'XLim');
%                    yy = get(gca, 'YLim');
%                    text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+0.9*(yy(2)-yy(1)), ...
%                        ['p=' num2str(pval)], 'fontsize', 6)
%                    axis([xx yy])
%                end 
            end
%             set(gcf, 'Name', 'GLM weight distribution', 'NumberTitle', 'off') 
            j = j + c + 1;   
        end
    end
end

%%
% spike-triggered accrage LFP ===================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'sta'))==1
    disp('STA analysis ----------------------------------------')
    figure(j+1);
    sta_t = linspace(-datast{end}.wnd, datast{end}.wnd, ...
        size(datast{end}.cond(1).sta.mean, 2));
    sta_analysiswnd = sta_t >= -0.05 & sta_t <= 0.05;
    sta_analysis_time = sta_t(sta_analysiswnd);
    pow_t = linspace(-datast{end}.wnd, datast{end}.wnd, ...
        length(datast{end}.cond(1).sta.t{end}));
    pow_analysiswnd = pow_t > -0.05 & pow_t < 0.05;
    freq = datast{end}.cond(1).sta.f{end};
    lenf = length(freq);
    yy = [0 0];
    xrange = [-0.07 0.07];
    cc = zeros(ss*2, 2);
    c = 0;
    for s = 1:ss
        % data extraction
        stas = zeros(lenses, size(datast{end}.cond(1).sta.mean, 2), 2);
        amp = zeros(lenses, 2);
        zerot = zeros(lenses, 2);
        pow = zeros(lenses, lenf, 2); 
        sz = size(datast{end}.cond(1).sta.p{end});
        staim = nan(lenses, sz(1), sz(2));
        nspk = zeros(lenses, 2);
        for i = 1:lenses
            for d = 1:2
                sta_temp = datast{i}.cond(d).sta.mean(stmidx(i, s), :);
                stas(i, :, d) = sta_temp;

                [amp(i, d), zerot(i, d)] = min(sta_temp(sta_analysiswnd));

                amp(i, d) = max(sta_temp(-0.1 < sta_t & sta_t < 0.1)) - ...
                    min(sta_temp(-0.1 < sta_t & sta_t < 0.1));
                zerot(i, d) = sta_analysis_time(zerot(i, d));
                pow(i, :, d) = nanmean(datast{i}.cond(d).sta.p{stmidx(i, s)}...
                    (:, pow_analysiswnd), 2);
                nspk(i, d) = datast{i}.cond(d).sta.nspk(stmidx(i, s)); 
            end
            staim(i, :, :) = (datast{i}.cond(1).sta.p{stmidx(i, s)} - datast{i}.cond(2).sta.p{stmidx(i, s)}).*sum(nspk(i,:));
        end

        % visualize
        figure
        for a = 1:lena
            % STA mean
            
            for k = 1:ndrug
                subplot(lena*ndrug, 4*ss, 1 + 4*(s - 1) + 4*ss*(k-1) + ndrug*4*ss*(a-1))
                if a < 3
                    cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                else
                    cond = lists(:,2)==k-1;
                end
                for d = 1:2
                    mat = squeeze(stas(cond, :, d));
                    nans = any(isnan(mat), 2);          
                    mat = mat(nans==0, :);
                    [me, sem] = weighted(mat, nspk(cond, d));
                    fill_between(sta_t, me - sem, me + sem, cols(d, :), 0.1)
                    hold on;
                    plot(sta_t, me, '-', 'color', cols(d, :))
                    hold on;            
                end
                yy_temp = get(gca, 'YLim');
                if yy_temp(1) < yy(1)
                    yy(1) = yy_temp(1);
                end
                if yy_temp(2) > yy(2)
                    yy(2) = yy_temp(2);
                end
                set(gca, 'box', 'off', 'tickdir', 'out')
                if k==1 && a==1
                    title(['STA in ' stmlab{s}])            
                end
                if s==1
                    ylabel({animals{a}, pairnames{k, 2}, 'LFP (uV)'})
                    if k == ndrug && a==lena
                        xlabel('time from spike (sec)')     
                    end
                end

                % zero-amplitude
                subplot(lena*ndrug, 4*ss, 2 + 4*(s - 1) + 4*ss*(k-1) + ndrug*4*ss*(a-1))
                try
                    unity_scatter(amp(cond, 1), amp(cond, 2))
                catch
                end
                if k==1 && a==1
                    title('STA amplitude')
                end

                % min-time
                subplot(lena*ndrug, 4*ss, 3 + 4*(s - 1) + 4*ss*(k-1) + ndrug*4*ss*(a-1))
                try
                    unity_scatter(zerot(cond, 1), zerot(cond, 2))
                catch
                end
                if k==1 && a==1
                    title('time of minima')
                end

                % delta PSD
                subplot(lena*ndrug, 4*ss, 4 + 4*(s - 1) + 4*ss*(k-1)+ ndrug*4*ss*(a-1))
                m = squeeze(nansum(staim(cond, :, :), 1))./sum(sum(nspk(cond, :)));
               
                frange = freq >= 3 & freq <= 48;
                imgt = linspace(-0.5, 0.5, size(m, 2));
                trange = imgt >= xrange(1) & imgt <= xrange(2);
                ff = freq(frange);
                imagesc(imgt(trange), ff, m(frange, trange));
                colormap(jet)
                hold on;
                cc(c+1, :) = caxis;
                c = c + 1;
                yl = get(gca, 'YLim');
                plot([0 0], yl, '--w')
                if k==ndrug && a==lena
                    xlabel('time from spike (sec)')
                    ylabel('frequency (Hz)')
                end
                set(gca, 'box', 'off', 'tickdir', 'out')

                % stats for sta power
                subplot(lena*ndrug, 4*ss, 3 + 4*(s - 1) + 4*ss*(k-1)+ ndrug*4*ss*(a-1))
          
            end
               
        end
    end
    % align y-axis & same color range
    cc = [min(cc(:,1)), max(cc(:,2))];
    for k = 1:lena*ndrug*ss*4
        if mod(k, 4)==1
            subplot(lena*ndrug, 4*ss, k)
            plot([0 0], yy, '--k')
            ylim(yy)
            xlim(xrange)
        end
        if mod(k, 4)==0
            subplot(lena*ndrug, 4*ss, k)
            caxis(cc)
            c = colorbar('northoutside');
            if ndrug > 1
                c.Label.String = '\Delta PSD (drug - base)';
            else
                c.Label.String = ['\Delta PSD (' pairnames{1, 1} ' - ' pairnames{1, 2} ')'];
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
    bandnames = {'low-freq (3-10)', 'alpha (8-12)', 'beta (15-25)', 'gamma (40-48)'};
    bandrange = {[3, 10], [8, 12], [15, 25], [40 48]};
    lenb = length(bandnames);
    freq = datast{end}.cond(1).spectrogram.f{1};
    if size(freq, 1) > size(freq, 2)
        freq = freq';
    end
    lenf = length(freq);
    spe_t = datast{end}.cond(1).spectrogram.t{1};
    lent = length(spe_t);
    spe_t = linspace(-0.1, stmdur, lent);
    yy = {[0 0], [0 0]};
%     step = 3:1:48; % Hz
    for s = 1:ss
        % data extraction
        ntr = zeros(lenses, 2);
        power = zeros(lenses, length(freq), 2);
        delta = zeros(lenses, length(freq), 2);
        Sall = zeros(lenses, lenf, lent, 2);
        for i = 1:lenses
            for d = 1:2
                ntr(i, d) = size(datast{i}.cond(d).mat{1}, 1);
                sz = size(datast{i}.cond(d).spectrogram.S{stmidx(i, s)});
                if sz(2) > sz(1)
%                     S = 10*log10(datast{i}.cond(d).spectrogram.S{stmidx(i, s)});
                    S = datast{i}.cond(d).spectrogram.S{stmidx(i, s)};
                else
%                     S = 10*log10(datast{i}.cond(d).spectrogram.S{stmidx(i, s)})';
                    S = datast{i}.cond(d).spectrogram.S{stmidx(i, s)}';
                end 
                if sum(isnan(S(:)))==0
                    S = imresize(S, [lenf, lent]);
%                     Sres = S - repmat(nanmean(S(:, spe_t < 0), 2), 1, lent);
%                     Sall(i, :, :, d) = Sres(:, end - length(spe_t)+1:end);
                    Sall(i, :, :, d) = S(:, end - length(spe_t)+1:end).*ntr(i, d);                    
                end
                Sabs = nanmean(S(:, spe_t > datast{end}.window{end}(1)), 2);
                power(i, :, d) = Sabs;
                delta(i, :, d) = Sabs - nanmean(S(:, spe_t < 0), 2);
            end
        end

        % visualize ===========================
        for k = 1:ndrug
            statistics.drug(k).occrall_pow = cell(lenb, 2);
            statistics.drug(k).response_pow = cell(lenb, 2);
        end
        for a = 1:lena
            for k = 1:ndrug
                % occrall power ====================
                figure(j + 1);
                subplot(lena*ndrug, 2*ss, 1 + 2*(s-1) + 2*ss*(k-1) + ndrug*2*ss*(a-1))
                if a < 3
                    cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                else
                    cond = lists(:,2)==k-1;
                end
                nans = [];
                mat = cell(1, 2);
                for d = 1:2
                    % power
                    mat{d} = squeeze(power(cond, :, d));
                    nans = [nans; find(any(isnan(mat{d}), 2))];
                end
                nans = unique(nans);
                for d = 1:2        
                    mat{d}(nans, :) = [];
                    [me, sem] = weighted(mat{d}, ntr(cond, d));
%                     me = mean(mat{d}, 1);
%                     sem = std(mat{d}, [], 1)/sqrt(lenses - length(nans));
                    fill_between(freq, me - sem, me + sem, cols(d, :), 0.4)
                    hold on;
                    plot(freq, me, '-', 'color', cols(d, :))
                    hold on;                                
                    
                    % for stats
%                     for b = 1:length(freq)
%                         if b == length(step)-1
%                             sz = step(2)-step(1)+1;
%                         else
%                             sz = step(2)-step(1);
%                         end
%                         statistics.drug(k).occrall_pow{b, d} = mat{d}(:, step(b) <= freq & ...
%                             step(b)+sz > freq);
%                         statistics.drug(k).occrall_pow{b, d} = mat{d}(:, b);
%                     end
                    for b = 1:lenb
                        statistics.drug(k).occrall_pow{b, d} = mat{d}(:, bandrange{b}(1) <= freq & ...
                            bandrange{b}(2) >= freq);
                    end
                end

%                 % stats
%                 sfreq = freq(freq >= 3 & freq <= 48);
                yy_temp = get(gca, 'YLim');
                yy_temp(2) = yy_temp(2)*1.1;
                hold on;
                plot([10 10], yy_temp, '-k')
%                 start = 3;
%                 for b = 1:length(sfreq)
%                     % ttest
%                     pval = signrank(nanmean(statistics.drug(k).occrall_pow{b, 1}, 2), ...
%                         nanmean(statistics.drug(k).occrall_pow{b, 2}, 2));
%                     if pval < 0.05/length(sfreq)      
%                           pb = plot([start (freq(b)+freq(b+1))/2], [1 1]*yy_temp(2), ...
%                               '-', 'color', 0.5*[0 1 0], 'linewidth', 2);
%                           pb.Color(4) = 0.4;
%                     end
%                     start = (freq(b)+freq(b+1))/2;
%                 end
                % for stats
                for b = 1:lenb
                    pval = signrank(nanmean(statistics.drug(k).occrall_pow{b, 1}, 2), ...
                        nanmean(statistics.drug(k).occrall_pow{b, 2}, 2));
                    disp(['cond ' num2str(k) ':' bandnames{b} ', ' num2str(a) ', p=' num2str(pval*(2*lenb))])
                    if pval < 0.05/(2*lenb)    
                        disp([animals{a} ', ' bandnames{b} ', p = ' num2str(pval*lenb)])
                          pb = plot(bandrange{b}, [1 1]*yy_temp(2), ...
                              '-', 'color', [0.8000    0.9216    0.7725], 'linewidth', 2);
                          pb.Color(4) = 0.4;
                    end
                end

                % format
                if yy_temp(1) < yy{1}(1)
                    yy{1}(1) = yy_temp(1);
                end
                if yy_temp(2) > yy{1}(2)
                    yy{1}(2) = yy_temp(2);
                end
%                 set(gca, 'YScale', 'log')
                set(gca, 'XScale', 'log')
                xlim([3 48])
                set(gca, 'box', 'off', 'tickdir', 'out', 'XScale', 'log')
                if k==1
                    title(stmlab{s})            
                end
                if s==1
                    ylabel({pairnames{k, 2}, 'power'})
                    if k == ndrug
                        xlabel('frequency (Hz)')     
                    end
                end

                % baseline corrected =========================
                subplot(lena*ndrug, 2*ss, 2 + 2*(s-1) + 2*ss*(k-1) + ndrug*2*ss*(a-1))
                nans = [];
                mat = cell(1, 2);
                for d = 1:2
                    % power
                    mat{d} = squeeze(delta(cond, :, d));
                    nans = [nans; find(any(isnan(mat{d}), 2))];
                end
                nans = unique(nans);
                for d = 1:2
                    mat{d}(nans, :) = [];
                    [me, sem] = weighted(mat{d}, ntr(cond, d));
%                     me = mean(mat{d}, 1);
%                     sem = std(mat{d}, [], 1)/sqrt(lenses - length(nans));
                    fill_between(freq, me - sem, me + sem, cols(d, :), 0.1)
                    hold on;
                    plot(freq, me, '-', 'color', cols(d, :))
                    hold on;            

                    % for stats
                    for b = 1:lenb
                        statistics.drug(k).response_pow{b, d} = mat{d}(:, bandrange{b}(1) <= freq & ...
                            bandrange{b}(2) >= freq);
                    end
                end

                % stats
                pbar = nan(1, length(freq));
                for b = 1:lenb
                    % ttest
                    [~, pval] = ttest(mean(statistics.drug(k).response_pow{b, 1}, 2), ...
                        mean(statistics.drug(k).response_pow{b, 2}, 2));
                    if pval < 0.05/lenb
                        pbar(bandrange{b}(1) <= freq & bandrange{b}(2) >= freq) = 1;
                    end
                end
                yy_temp = get(gca, 'YLim');
                plot(freq, pbar*yy_temp(2), '-', 'color', [0.8000    0.9216    0.7725], 'linewidth', 3)

                % format
                if yy_temp(1) < yy{2}(1)
                    yy{2}(1) = yy_temp(1);
                end
                if yy_temp(2) > yy{2}(2)
                    yy{2}(2) = yy_temp(2);
                end
                set(gca, 'XScale', 'log')        
                xlim([3 48])
                set(gca, 'box', 'off', 'tickdir', 'out')
                if k==1
                    title(stmlab{s})            
                end
                if s==1
                    ylabel({animals{a}, pairnames{k, 1}, '\Delta power'})
                    if k == ndrug
                        xlabel('frequency (Hz)')     
                    end
                end

                % average PSDs ==========================
                figure(j + 1 + s);
                clim = {zeros(1, 2), zeros(1, 2)};
                for d = 1:2
                    subplot(lena*ndrug, 3, d + 3*(k-1) + ndrug*3*(a-1))
%                     imagesc(spe_t, freq, squeeze(nanmean(Sall(cond, :, :, d), 1)));
                    h = pcolor(spe_t, freq, squeeze(nansum(Sall(cond, :, :, d), 1))./sum(ntr(cond, d)));
                    h.EdgeColor = 'none';
%                     set(gca, 'YScale', 'log')
                    colormap(jet);
                    cl = caxis;
                    if cl(1) < clim{1}(1)
                        clim{1}(1) = cl(1);
                    end
                    if cl(2) > clim{1}(2)
                        clim{1}(2) = cl(2);
                    end
                    hold on;
                    yl = get(gca, 'YLim');
                    plot([0 0], yl, '--w')
                end        

                subplot(lena*ndrug, 3, 3 + 3*(k-1) + ndrug*3*(a-1))
                S1 = squeeze(nansum(Sall(cond, :, :, 1), 1)./sum(ntr(cond, 1)));
%                 S1 = S1 - repmat(nanmean(S1(:, spe_t <= 0), 2), 1, length(spe_t));
                S2 = squeeze(nansum(Sall(cond, :, :, 2), 1)./sum(ntr(cond, 2)));
%                 S2 = S2 - repmat(nanmean(S2(:, spe_t <= 0), 2), 1, length(spe_t));
%                 imagesc(spe_t, freq, S1 - S2);
                h = pcolor(spe_t, freq, S1 - S2);
                h.EdgeColor = 'none';
%                 set(gca, 'YScale', 'log')
                colormap(jet)
                cl = caxis;
                if cl(1) < clim{2}(1)
                    clim{2}(1) = cl(1);
                end
                if cl(2) > clim{2}(2)
                    clim{2}(2) = cl(2);
                end
                ylim([3 48])
                hold on;
                yl = get(gca, 'YLim');
                plot([0 0], yl, '--w')
            end
        end

        % stats across conditions (ttest2)
        if ndrug>1
            for b = 1:lenb
                disp(['--- ' stmlab{s} '(' bandnames{b} ') -----------------------'])

                % occrall power
                Xcond = mean(statistics.drug(1).occrall_pow{b, 1}, 2) - mean(statistics.drug(1).occrall_pow{b, 2}, 2);
                Ycond = mean(statistics.drug(2).occrall_pow{b, 1}, 2) - mean(statistics.drug(2).occrall_pow{b, 2}, 2);
                [~, pval] = ttest2(Xcond, Ycond);
                disp(['occrall power (cond1 vs cond2): p(ttest2) = ' num2str(pval)])

                % response power
                Xcond = mean(statistics.drug(1).response_pow{b, 1}, 2) - mean(statistics.drug(1).response_pow{b, 2}, 2);
                Ycond = mean(statistics.drug(2).response_pow{b, 1}, 2) - mean(statistics.drug(2).response_pow{b, 2}, 2);
                [~, pval] = ttest2(Xcond, Ycond);
                disp(['response power (cond1 vs cond2): p(ttest2) = ' num2str(pval)])
            end
        end
        
        % format 
        disp(['clim ' num2str(clim{2})])
        for k = 1:lena*ndrug*3
            figure(j + 1 + s);
            subplot(lena*ndrug, 3, k)
%             c = colorbar('northoutside');
            if mod(k, 3) == 1
                caxis(clim{1})
%                 c.Label.String = 'PSD (base)';
            elseif mod(k, 3) == 2
                caxis(clim{1})
%                 c.Label.String = 'PSD (drug)';
            elseif mod(k, 3) == 0
                caxis(clim{2})
%                 c.Label.String = 'PSD: base - drug';
            end
            if k==1
                ylabel({animals{a}, pairnames{k, 1}, 'frequency (Hz)'})
                xlabel('time after stimulus onset (sec)')
            end
            if k==4
                ylabel({animals{a}, pairnames{2, 1}, 'frequency (Hz)'})
                xlabel('time after stimulus onset (sec)')
            end
            xl = get(gca, 'XLim');
            xlim([-0.1 xl(2)])
            set(gca, 'box', 'off', 'tickdir', 'out')
        end    
    end
    % align y-axis
    for k = 1:lena*ndrug*ss*2
        figure(j+1);
        subplot(lena*ndrug, 2*ss, k)    
        if mod(k, 2)==0
            ylim(yy{2})
        else
            ylim(yy{1})
        end
    end
    figure(j+1);
    set(gcf, 'Name', 'PSD', 'NumberTitle', 'off')
    for s = 1:ss
        figure(j+1+s);
        set(gcf, 'Name',  ['accrage spectrogram: ' ...
            datast{1}.stm.param ' = ' stmlab{s}], 'NumberTitle', 'off')
    end
    j = j + 1 + ss;
end

%%
% coherence ==============================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'coherence'))==1
    disp('coherence analysis --------------------------------------')
    % for circular statistics
    addpath(genpath([mypath '/Katsuhisa/code/integrated/CircStat2012a']))
    switch datast{1}.stm.param
        case 'rc'
            bandnames = {'theta (3-7)', 'alpha (8-12)', 'beta (15-25)', 'gamma (40-48)'};
            bandrange = {[3, 7], [8, 12], [15, 25], [40, 48]};
        otherwise
            bandnames = {'< theta (0-7)', 'alpha (8-13)', 'beta (14-24)', 'gamma (36-48)'};
            bandrange = {[0, 7], [8, 13], [14, 24], [36, 40], [40, 44], [44, 48]};
    end
    lenb = length(bandnames);
    xx = [3 10];
    yy = [0.1 0.4];
    freq = datast{end}.cond(1).coherence.f{end}';
    lenf = length(freq);
    analysisfrange = freq(freq >= xx(1) & freq <= xx(2));
    coh = cell(1, ss);
    for s = 1:ss
        % data extraction
        coh{s} = zeros(lenses, lenf, 2);
        phi = zeros(lenses, lenf, 2);
        for i = 1:lenses
            for d = 1:2
                coh{s}(i, :, d) = datast{i}.cond(d).coherence.C{stmidx(i, s)};
%                 phases = datast{i}.cond(d).coherence.phi{stmidx(i, s)};
%                 for j = 1:size(phases, 1)
%                     nans = isnan(phases(j, :));
%                     nonan = phases(j, nans==0);
%                     phases(j, nans) = circ_mean(nonan, [], 2);
%                 end                
%                 phi(i, :, d) = circ_mean(phases, [], 2);
                phi(i, :, d) = datast{i}.cond(d).coherence.phi{stmidx(i, s)};
            end
        end

        % visualize
        for k = 1:ndrug
            statistics.drug(k).coh = cell(lenf, 2);
            statistics.drug(k).phi = cell(lenf, 2);
        end
        for a = 1:lena
            for k = 1:ndrug
                if a < 3
                    cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                else
                    cond = lists(:,2)==k-1;
                end
                figure(j + 1);
                subplot(lena*ndrug, ss, s + ss*(k-1) + ndrug*ss*(a-1))
                matc = cell(1, 2);
                for d = 1:2
                    % coherence
                    matc{d} = squeeze(coh{s}(cond, :, d));
                end
                cohd = nan(d, lenf);
                for d = 1:2
                    cohd(d, :) = nanmean(matc{d}, 1);
                    sem = nanstd(matc{d}, [], 1)/sqrt(lenses);
                    fill_between(freq, cohd(d, :) - sem, cohd(d, :) + sem, cols(d, :), 0.1)
                    hold on;
                    plot(freq, cohd(d, :), '-', 'color', cols(d, :))
                    hold on; 
%                     singcol = jet(size(matc{d}, 1));
%                     for l = 1:size(matc{d}, 1)
%                         plot(freq, matc{d}(l, :), '-', 'color', singcol(l, :))
%                         hold on;
%                     end

                    % for stats
                    for b = 1:lenb
                        if b < lenb
                            frange = bandrange{b}(1) <= freq & bandrange{b}(2) >= freq;
                        else % gamma
                            frange = bandrange{b}(1) <= freq & 80 >= freq;
                        end
                        statistics.drug(k).coh{b, d} = matc{d}(:, frange);
                    end
                end

                % format        
                ylim(yy)
%                 yy_temp = get(gca, 'YLim');
%                 if yy_temp(1) < yy(1)
%                     yy(1) = yy_temp(1);
%                 end
%                 if yy_temp(2) > yy(2)
%                     yy(2) = yy_temp(2);
%                 end
                set(gca, 'box', 'off', 'tickdir', 'out')
%                 xlim([3 48])
                xlim([3 13])
%                 set(gca, 'XScale', 'log')
%                 set(gca, 'YScale', 'log')
                if k==1 && a==1
                    title(stmlab{s})            
                end
                if s==1
                    ylabel({animals{a}, pairnames{k, 2}, 'coherence'})
                    if k == ndrug && a==lena
                        xlabel('frequency (Hz)')     
                    end
                end

                figure(j + 1 + s);
                rad = cell(lenb, 2);
                for d = 1:2
                    % phase
                    matp = squeeze(phi(cond, :, d));
                    for b = 1:lenb
                        % plot
                        subplot(lena*ndrug, lenb, b + lenb*(k-1) + ndrug*lenb*(a-1))
%                         if b==lenb
%                            frange = bandrange{b}(1) <= freq & ...
%                                 bandrange{b}(2) >= freq;
%                            rad1 = maxcoh_rad(cohd(frange), matp(:, frange));
%                            frange = bandrange{b+1}(1) <= freq & ...
%                                 bandrange{b+1}(2) >= freq;                           
%                            rad2 = maxcoh_rad(cohd(frange), matp(:, frange));
%                            frange = bandrange{b+2}(1) <= freq & ...
%                                 bandrange{b+2}(2) >= freq;
%                            rad3 = maxcoh_rad(cohd(frange), matp(:, frange));
%                            rad{b, d} = circ_axial([rad1; rad2; rad3]);
%                         else
                            frange = bandrange{b}(1) <= freq & ...
                                bandrange{b}(2) >= freq;
                            rad1 = maxcoh_rad(cohd(frange), matp(:, frange));
                            rad{b, d} = circ_axial(rad1);
%                         end
%                         rad_nonan = rad{b, d}(~isnan(rad{b, d}));
        %                 rad{b, d}(isnan(rad{b, d})) = circ_median(rad_nonan);
%                         rad{b, d}(isnan(rad{b, d})) = circ_mean(rad_nonan);
%                         pol = polarhistogram(rad{b, d});
                        pol = circ_plot(rad{b, d},'hist',[],20,false,true,'linewidth',3,'color',cols(d,:));
%                         pol.BinWidth = 0.2;
%                         pol.DisplayStyle = 'stairs';
                        for l = 1:2
                            pol.Children(l).Color = cols(d, :);
                        end
                        hold on;                
                    end
                end
                % circular statistics
                statistics.drug(k).phi = rad;
                for b = 1:lenb
                    subplot(lena*ndrug, lenb, b + lenb*(k-1) + ndrug*lenb*(a-1))
        %             disp(['circular ww test: ' 'drug; ' drugnames{k} ' ' bandnames{b}])
                    [rx, ry] = nan_remove_pair(rad{b, 1}, rad{b, 2}, 'mean');
                    [ppval, s_table] = circ_wwtest(rx, ry);
%                     pdelta = circ_rad2ang(circ_mean(rad{b,1}) - circ_mean(rad{b,2}));
                    npval = circ_cmtest(rx, ry);
%                     ndelta = circ_rad2ang(circ_median(rad{b,1}) - circ_median(rad{b,2}));
                    if a==1 && k==1
                        title({[datast{end}.stm.param ' = ' stmlab{s} ' (' bandnames{b} ')'], ...
                            ['p(ww)=' num2str(ppval) ', p(cm)=' num2str(npval)]}, 'FontSize', 6)  
                    else
                        title(['p(ww)=' num2str(ppval) ', p(cm)=' num2str(npval)], 'FontSize', 6)  
                    end
        %             pval = circ_cmtest(rad{b, 1}, rad{b, 2});
        %             delta = circ_rad2ang(circ_median(rad{b,1}) - circ_median(rad{b,2}));               
        %             title({[datast{end}.stm.param ' = ' stmlab{s} ' (' bandnames{b} ')'], ...
        %                 ['p(cm)=' num2str(pval) ', \Delta=' num2str(delta) '^o']})  
                end
            end
            
            % stats
            figure(j + 1);
%             pbar = nan(1, length(analysisfrange));
%             for b = 1:length(analysisfrange)
%                 pval = signrank(squeeze(coh{s}(:,analysisfrange==analysisfrange(b),1)), ...
%                     squeeze(coh{s}(:,analysisfrange==analysisfrange(b),2)));
%                 if pval < 0.05/length(analysisfrange)
%                     hold on;
%                     plot
%                 end
%             end
%             for b = 1:length(freq)
                
%             for b = 1:lenb
%                 % ttest
%                 if ndrug == 1
%                     stats = pair_tests([mean(statistics.drug(1).coh{b, 1}, 2), ...
%                         mean(statistics.drug(1).coh{b, 2}, 2)]);
%                 else
%                     stats = pair_tests([mean(statistics.drug(1).coh{b, 1}, 2), ...
%                         mean(statistics.drug(1).coh{b, 2}, 2)], ...
%                         [mean(statistics.drug(2).coh{b, 1}, 2), ...
%                         mean(statistics.drug(2).coh{b, 2}, 2)]);
%                 end
% 
%                 % within (cond 1)
%                 subplot(lena*ndrug, ss, s+ ndrug*ss*(a-1))
%                 if stats.pair(1).signrank.p < 0.05/lenb
%                     hold on;
%                     if b < lenb
%                         plot(bandrange{b}, yy(2)*[1 1],  '-', 'color', 0.5*[1 1 1], 'linewidth', 2)
%                     else
%                          plot([bandrange{b}(1) 80], yy(2)*[1 1],  '-', 'color', 0.5*[1 1 1], 'linewidth', 2)
%                     end
%                 end
%                 
%                 if ndrug > 1
%                     % within (cond 2)
%                     subplot(lena*ndrug, ss, s + ss*(ndrug-1) + ndrug*ss*(a-1))
%                     disp([])
%                     if stats.pair(2).signrank.p < 0.05/lenb
%                         hold on;
%                         if b < lenb
%                             plot(bandrange{b}, yy(2)*[1 1],  '-', 'color', 0.5*[1 1 1], 'linewidth', 2)
%                         else
%                              plot([bandrange{b}(1) 80], yy(2)*[1 1],  '-', 'color', 0.5*[1 1 1], 'linewidth', 2)
%                         end
%                     end  
%                     
%                     % across (cond 1 vs 2)
%                     disp([bandnames{b} ': NaCl vs 5HT; p=' num2str(stats.ranksum.p)])
%                     if stats.ranksum.p < 0.05/lenb
%                         hold on;
%                         if b < lenb
%                             plot(bandrange{b}, (0.9*(yy(2)-yy(1)) + yy(2))*[1 1],  '-', 'color', 0.5*[1 0 0], 'linewidth', 2)
%                         else
%                              plot([bandrange{b}(1) 80], (0.9*(yy(2)-yy(1)) + yy(2))*[1 1],  '-', 'color', [1 0 0], 'linewidth', 2)
%                         end
%                     end 
%                 end
%             end
%             yy_temp = get(gca, 'YLim');
%             plot(freq, pbar*yy_temp(2), '-', 'color', 0.5*[1 1 1], 'linewidth', 3)
        end
    end
    
    % align y-axis
    for k = 1:lena*ndrug*ss
        figure(j+1);
        subplot(lena*ndrug, ss, k)    
        ylim(yy)
    end
    set(gcf, 'Name', 'spike-LFP coherence', 'NumberTitle', 'off')
    for s = 1:ss
        figure(j+1+s);
        set(gcf, 'Name', ['spike-LFP phase: ' datast{1}.stm.param ' = ' stmlab{s}], 'NumberTitle', 'off')
    end
    j = j + 1 + ss;
end


%%
% subfunction
function datast = encoding_tuning4rc(datast, smlinfo)
% compute encoding indices from tuning curcc
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
        
        % noise correlation
        datast{i}.cond(d).nc = smlinfo.paramat(sesidx(i), 9+d);
    end
end

function y = adjust_veclen(x, l)
lenx = length(x);
if lenx > l
    y = x(1:l);
else
    try
        y = [x, x(end)*ones(1, l-lenx)];
    catch
        y = [x; x(end)*ones(1, l-lenx)];
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

function mout = mymodel(x, base, predictors)
lenx = length(x);
s = 0;
for i = 1:lenx
    s = s + x(i)*predictors(:, i);
end
mout = base*exp(s);

function f = cost(x, y, base, predictors)
mout = mymodel(x, base, predictors);
f = 0;
for i = 1:size(predictors, 1)
    if mout(i) > 0
        f = f - (y(i)*log(mout(i)) - mout(i));
    end
%     f = f + abs(mout(i) - y(i));
end