function [data, ceyes] = eyepos_vs_ps(x, y, p)
% check whether eye positions are affected by pupil size
% after Choe et al., 2016
% INPUT:
% x ... x position of the eye
% y ... y position of the eye
% p ... pupil size
%
% OUTPUT:
% data ... output structure
% ceyes ... corrected eye traces
%

%%
% % trial average
% [xme, yme, pme] = trial_average(x, y, p, tr);

%%
% correlation
data.original = compute_correlation(x, y, p);

%% 
% linear regression
data.reg(1).xparas = glmfit(p', x', 'normal', 'link', 'identity', 'constant', 'on');
data.reg(1).yparas = glmfit(p', y', 'normal', 'link', 'identity', 'constant', 'on');

x_new = x' - glmval(data.reg(1).xparas, ...
    p', 'identity', 'constant', 'on');
y_new = y' - glmval(data.reg(1).yparas, ...
    p', 'identity', 'constant', 'on');
ceyes = [x_new, y_new];
data.reg(1).corrected = compute_correlation(x_new', y_new', p);

%%
% 2nd polynomial regression
data.reg(2).xparas = polyfit(p', x', 2);
data.reg(2).yparas = polyfit(p', y', 2);

x_new = x' - polyval(data.reg(2).xparas, p');
y_new = y' - polyval(data.reg(2).yparas, p');
ceyes = [ceyes, x_new, y_new];
data.reg(2).corrected = compute_correlation(x_new', y_new', p);

ceyes = ceyes';

% function [xme, yme, pme] = trial_average(x, y, p, tr)
% unitr = unique(tr);
% lentr = length(unitr);
% xme = nan(1, lentr);
% yme = nan(1, lentr);
% pme = nan(1, lentr);
% for i = 1:lentr
%     idx = tr == i;
%     xme(i) = mean(x(idx));
%     yme(i) = mean(y(idx));
%     pme(i) = mean(p(idx));
% end

function out = compute_correlation(xme, yme, pme)
[r, pval] = corrcoef(xme, pme);
out.Pearson.r(1) = r(1,2);
out.Pearson.p(1) = pval(1,2);
[r, pval] = corrcoef(yme, pme);
out.Pearson.r(2) = r(1,2);
out.Pearson.p(2) = pval(1,2);
[out.Spearman.r(1), out.Spearman.p(1)]...
    = corr(xme', pme', 'type', 'Spearman');
[out.Spearman.r(2), out.Spearman.p(2)]...
    = corr(yme', pme', 'type', 'Spearman');