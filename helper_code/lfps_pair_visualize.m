function lfps_pair_visualize(lfps, analysis, splittype)

% default values
if nargin < 2; analysis = {'all'}; end
if nargin < 3; splittype = 'drug'; end

% drug or FR control condition ===========================
ndrug = 1;
switch splittype
    case 'drug'
        ndrug = 2;
        pairnames = {'base', 'NaCl'; 'base', '5HT'};
    case 'sc'
        pairnames = {'large sc', 'small sc'; 'large sc', 'small sc'};
end

fnames = lfps.lfplist;
datast = lfps.LFP_prepro;
lists = [lfps.is5ht, lfps.animal];

animals = {'kaki', 'mango', 'both'};
if ndrug==1 
    lists(:,1) = 0;
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


% stimulus type ============================
stmidx = ones(lenses, 1);
stmlab = {'ORxRC'};
stmdur = 2;

ss = size(stmidx, 2);
cols = [0 0 0; 1 0 0];
j = 0;


%%
% spike-triggered average LFP ===================================
if sum(contains(analysis, 'all'))==1 || (sum(contains(analysis, 'sta') && ~contains(analysis, 'batch')))==1
    figure(j+1);

    sta_t = linspace(-datast{end}.wnd, datast{end}.wnd, size(datast{end}.cond(1).sta.mean, 2));
    sta_analysiswnd = sta_t >= -0.05 & sta_t <= 0.05;
    sta_analysis_time = sta_t(sta_analysiswnd);

    pow_t = linspace(-datast{end}.wnd, datast{end}.wnd, length(datast{end}.cond(1).sta.t{end}));
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
                amp(i, d) = max(sta_temp(-0.1 < sta_t & sta_t < 0.1)) - min(sta_temp(-0.1 < sta_t & sta_t < 0.1));
                
                zerot(i, d) = sta_analysis_time(zerot(i, d));
                pow(i, :, d) = nanmean(datast{i}.cond(d).sta.p{stmidx(i, s)}(:, pow_analysiswnd), 2);
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
                    cond = lists(:,1)==k-1 & lists(:,2)==a-1;
                else
                    cond = lists(:,1)==k-1;
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
% spike-triggered average LFP ===================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'sta_batch'))==1
    
    figure(j+1);
    
    sta_t = linspace(-datast{end}.wnd, datast{end}.wnd, size(datast{end}.cond(1).sta.mean, 2));
    xrange = [-0.1 0.1];
    
    nrow = floor(sqrt(lenses));
    ncol = ceil(lenses/nrow);
    
    for s = 1:ss
        for i = 1:lenses

            % plot
            subplot(nrow, ncol, i)
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

    freq = datast{end}.cond(1).spectrogram.f{1};
    if size(freq, 1) > size(freq, 2)
        freq = freq';
    end
    % Need to add a small epsilon value to 0 since log scale removes the plot
    freq_fill  = freq;
    if freq(1) == 0
        freq_fill(1) = freq_fill(1) + 10^-10;
    end
    lenf = length(freq);

    spe_t = datast{end}.cond(1).spectrogram.t{1};
    lent = length(spe_t);
    spe_t = linspace(-0.1, stmdur, lent);
    yy = {[0 0], [0 0]};

    for s = 1:ss
        % data extraction
        ntr = zeros(lenses, 2);
        power = zeros(lenses, length(freq), 2);
        delta = zeros(lenses, length(freq), 2);
        Sall = zeros(lenses, lenf, lent, 2);
        for i = 1:lenses
            for d = 1:2
                ntr(i, d) = datast{i}.cond(d).ntr;
                sz = size(datast{i}.cond(d).spectrogram.S{stmidx(i, s)});
                if sz(2) > sz(1)
                    S = datast{i}.cond(d).spectrogram.S{stmidx(i, s)};
                else
                    S = datast{i}.cond(d).spectrogram.S{stmidx(i, s)}';
                end 
                if sum(isnan(S(:)))==0
                    S = imresize(S, [lenf, lent]);
                    Sall(i, :, :, d) = S(:, end - length(spe_t)+1:end).*ntr(i, d);                    
                end
                Sabs = nanmean(S(:, spe_t > datast{end}.window{end}(1)), 2);
                power(i, :, d) = Sabs;
                delta(i, :, d) = Sabs - nanmean(S(:, spe_t < 0), 2);
            end
        end

        % visualize ===========================
        for a = 1:lena
            for k = 1:ndrug
                % occrall power ====================
                figure(j + 1);
                subplot(lena*ndrug, 2*ss, 1 + 2*(s-1) + 2*ss*(k-1) + ndrug*2*ss*(a-1))
                if a < 3
                    cond = lists(:,1)==k-1 & lists(:,2)==a-1;
                else
                    cond = lists(:,1)==k-1;
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
                                   
                    fill_between(freq_fill, me - sem, me + sem, cols(d, :), 0.4)
                    hold on;
                    plot(freq, me, '-', 'color', cols(d, :))
                    hold on;                 
                end
                
                % format
                yy_temp = get(gca, 'YLim');
                yy_temp(2) = yy_temp(2)*1.1;
                if yy_temp(1) < yy{1}(1)
                    yy{1}(1) = yy_temp(1);
                end
                if yy_temp(2) > yy{1}(2)
                    yy{1}(2) = yy_temp(2);
                end
                
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

                    fill_between(freq_fill, me - sem, me + sem, cols(d, :), 0.1)
                    hold on;
                    plot(freq, me, '-', 'color', cols(d, :))
                    hold on;            
                end

                % format
                yy_temp = get(gca, 'YLim');
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
                    
                    h = pcolor(spe_t, freq, squeeze(nansum(Sall(cond, :, :, d), 1))./sum(ntr(cond, d)));
                    h.EdgeColor = 'none';
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
                S2 = squeeze(nansum(Sall(cond, :, :, 2), 1)./sum(ntr(cond, 2)));
                
                h = pcolor(spe_t, freq, S1 - S2);
                h.EdgeColor = 'none';
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
        
        % format 
        for k = 1:lena*ndrug*3
            figure(j + 1 + s);
            subplot(lena*ndrug, 3, k)

            if mod(k, 3) == 1
                caxis(clim{1})
            elseif mod(k, 3) == 2
                caxis(clim{1})
            elseif mod(k, 3) == 0
                caxis(clim{2})
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
        set(gcf, 'Name',  ['average spectrogram: ' datast{1}.stm.param ' = ' stmlab{s}], ...
            'NumberTitle', 'off')
    end
    j = j + 1 + ss;
end

%% subfunction
function tlab = fname2title(fname)
    sti = strfind(fname, '_');
    dots = strfind(fname, '.');
    tlab = [fname(1:sti(2)-1) ':' fname(sti(end)+1:dots(end)-1)];