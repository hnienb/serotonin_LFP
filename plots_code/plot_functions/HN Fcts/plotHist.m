function p_tt = plotHist(val, i_sub, cax, ax_spec)
% plots a histogram


if nargin < 3
    cax = gca;
    ax_spec = '';
elseif nargin < 4
    delete(allchild(cax));
    ax_spec = '';
elseif nargin == 4
    delete(allchild(cax));
    
end

if strcmp(ax_spec, 'log')
    nbin = 10;   N = inf;
    while max(N) > floor(length(val)/4) && nbin<20
        nbin = nbin+1;
        [N, binrng] = histcounts(log(val), nbin);
    end
    binrng = exp(binrng);
else
    
    nbin = 10;   N = inf;
    while max(N) > floor(length(val)/4) && nbin<20
        nbin = nbin+1;
        [N, binrng] = histcounts(val, nbin);
    end
    
end


h1 = histogram(cax, val(i_sub), binrng, 'EdgeColor', 'w', 'FaceColor', 'r');
hold on;
h2 = histogram(cax, val(~i_sub) ,binrng, 'EdgeColor', 'w', 'FaceColor', 'k');



[~, p_tt] = ttest2(val(i_sub), val(~i_sub));
if any(~i_sub) && any(i_sub)
    p_wil = ranksum(val(i_sub), val(~i_sub));
    [~, p20b] = ttest(val(~i_sub), 0);  psignb = signrank(val(~i_sub));
else
    p_wil = nan;
    p20b = nan;  psignb = nan;
end
[~, p20r] = ttest(val(i_sub), 0);   psignr = signrank(val(i_sub));


if any(val <= 0)
    geomn_r = nan;
    geomn_b = nan;
else
    geomn_r = geomean(val(i_sub));
    geomn_b = geomean(val(~i_sub));
end

med_r = nanmedian(val(i_sub));
med_b = nanmedian(val(~i_sub));
qua_r =  quantile(val(i_sub), [.25 .75]);
qua_b =  quantile(val(~i_sub), [.25 .75]);

mn_r = nanmean(val(i_sub));
mn_b = nanmean(val(~i_sub));
std_r = nanstd(val(i_sub));
std_b = nanstd(val(~i_sub));



if strcmp(ax_spec, 'log')
    [~, p20r] = ttest(log(val(i_sub)), 0);   psignr = signrank(log(val(i_sub)));
    
    if any(~i_sub)     
        [~, p20b] = ttest(log(val(~i_sub)), 0);  psignb = signrank(log(val(~i_sub)));
    else
        p20b = nan;  psignb = nan;
    end
    
    title(cax, sprintf(['p_{tt}=%1.3f, p_{wilc}=%1.3f \n n_{red}=%1.0f ,' ...
        'plog_{ttvs0}=%1.3f,  plog_{wilvs0}=%1.3f,  '...
        '\\mu_{geom}=%1.2f ' ...
        '\n n_{black}=%1.0f , plog_{ttvs0}=%1.3f,  plog_{wilvs0}=%1.3f, \\mu_{geom}=%1.2f '], ...
        p_tt, p_wil, length(val(i_sub)), p20r, psignr, geomn_r, ...
        length(val(~i_sub)), p20b, psignb, geomn_b), 'FontSize', 8, 'FontWeight', 'bold');
    
else
    
    title(cax, sprintf(['p_{tt}=%1.3f, p_{wilc}=%1.3f \n n_{red}=%1.0f , p_{ttvs0}=%1.3f,  p_{wilvs0}=%1.3f,  '...
        '\\mu=%1.2f +/- %1.2f, med=%1.2f qu = %1.2f, %1.2f' ...
        '\n n_{black}=%1.0f , p_{ttvs0}=%1.3f,  p_{wilvs0}=%1.3f, \\mu=%1.2f +/- %1.2f, med=%1.2f qu = %1.2f, %1.2f'], ...
        p_tt, p_wil, length(val(i_sub)), p20r, psignr, mn_r, std_r,  med_r, qua_r(1), qua_r(2),...
        length(val(~i_sub)), p20b, psignb,  mn_b, std_b,  med_b, qua_b(1), qua_b(2)), 'FontSize', 8, 'FontWeight', 'bold');
    
    
end


max_y = max(get(cax, 'Ylim'));
plot(cax, med_r, max_y, 'vk', 'MarkerFaceColor', 'k');
text(med_r, max_y-max_y/10, sprintf('%1.2f', med_r), 'Parent', cax,'FontWeight', 'bold', 'FontSize', 8);

plot(cax ,med_b, max_y, 'vk', 'MarkerFaceColor', 'k');
text(med_b, max_y-max_y/5, sprintf('%1.2f', med_b), 'Parent', cax, 'FontWeight', 'bold', 'FontSize', 8 );


box off


end