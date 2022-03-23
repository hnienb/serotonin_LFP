function [cor, ps, lag] = corrlag(vec1, vec2, maxlag, oneside, type, fig, varargin)
%% compute correlation coefficient with a time-lag (like cross-correlation)
%% without 0-padding
% INPUT: vec1 ... vector 1
%        vec2 ... vector 2
%        maxlag ... lag (same unit as the vector) to take into account (-maxlag to maxlag)
%        oneside ... if 1, compute CC only for left side
%        type ... either 'Pearson' or 'Spearman'
%        fig ... 1, plot
%
% OUTPUT: cor ... correlation coefficient as a function of time
%         ps ... p-values
%         lag ... time-lag
%
% written by Katsuhisa (01.03.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% deal with the input
if nargin < 2
    disp('compute autocorrelation')
    vec2 = vec1;
end
if nargin < 3
    maxlag = 1;
end
if nargin < 4
    oneside = 0;
end
if nargin < 5
    type = 'Pearson';
end
if nargin < 6
    fig = 1;
end

% vec1 and vec2 have to be the same length
if length(vec1) ~= length(vec2)
  error('2 vectors must be the same length.')
end

% forced vector length 
lenv = length(vec1) - maxlag;
% disp(['length of vector for computing correlation: ' num2str(lenv)])

% transpose if necessary
sz1 = size(vec1);
sz2 = size(vec2);
if sz1(1)==1
        vec1 = vec1';
end
if sz2(1)==1
        vec2 = vec2';
end

% initialization
if oneside==1
        lag = 1:maxlag;
else
        lag = -maxlag:maxlag;
end

cor = zeros(1,length(lag));
ps = zeros(1,length(lag));

% loop to compute correlation coefficients with time-lag
for i = 1:maxlag
        if oneside==1
                [cor(i),ps(i)] = corr(vec1(1:lenv),vec2(i:lenv+i-1),'type',type);
        else
                [cor(maxlag-i+1), ps(maxlag-i+1)] = corr(vec1(1:lenv),vec2(i:lenv+i-1),'type',type);
                [cor(maxlag+i), ps(maxlag+i)] = corr(vec2(1:lenv),vec1(i:lenv+i-1),'type',type);
        end
end

if cor(end)==0
        [cor(end), ps(end)] = corr(vec2(1:lenv), vec1(maxlag+1:lenv+maxlag),'type',type);
end

% plot figures
if fig==1
  figure;
  plot(lag, cor, '-k', 'linewidth', 1.5)
  xlabel('lag')
  ylabel([type ' correlation coefficient'])
  yy = get(gca, 'YLim');
  hold on;
  signif = nan(1,length(lag));
  signif(ps < (0.05/length(lag))) = yy(2) + 0.1;
  plot(lag, signif, '*r')
  set(gca, 'box', 'off')
  axis square
end
