function data = fitPM(x,y,n,func,method,fig)
%% 
% fit psychometric function by 'cumulative Weibull' or 'cumulative
% Gaussian'
% INPUT: x,y,n ... unique stimulus range, probability of choice 1, the
% number of choice 1 in each stimulus type
%        func ... 'Weibull','Gaussian'
%        method ... 'LSE'(least squared error), 'MLE' (maximum likelihood)
%        fig ... 0, no plot; 1; plot
%
% OUTPUT: data.params(1): slope (no sign) or bias (with sign); mean
%         data.params(2): psychophysical threshold; SD
%         data.params(3): the lapse rate (only for 'Weibull')
%         data.variance_explained
%         data.fitx, data.fity ... for visualization
%
% NOTE: this function requires 'fminsearchbnd'
% 
% written by Katsuhisa (24.01.2018)
% +++++++++++++++++++++++++++++++++++++++++

% sign of x
if all(x >= 0)
    xsign = 0;
else
    xsign = 1;
end

% initialized parameter
p0 = [median(x) 10*(max(x)-min(x))];

% fitted function
switch func
    case 'Weibull'
        if xsign==0
            % upper bound
            upbnd = [1 inf inf];
            % lower bound
            lwbnd = [0 0 0];
        elseif xsign==1
            % upper bound
            upbnd = [x(end)*2 180 0.45];
            % lower bound
            lwbnd = [-x(end)*2 0 0];
        end
        p0 = [p0 (y(1) + (1 - y(end)))/2];
    case 'Gaussian'
        if xsign==0
            % upper bound
            upbnd = [inf inf];
            % lower bound
            lwbnd = [0 0];
        elseif xsign==1
            % upper bound
            upbnd = [x(end)*2 180];
            % lower bound
            lwbnd = [-x(end)*2 0];
        end
end

%%
% cost function
c = @(p)cost(p, x, y, n, xsign, func, method);

%%
% parameter estimate
options = optimset('MaxFunEvals',10000,'maxiter',10000);
params = fminsearchbnd(c, p0, lwbnd, upbnd, options);

%%
% variance explained
s = unique(x);
fval = nan(1,length(s));
for i = 1:length(s)
    fval(i) = funcval(params, x(i), xsign, func);
end
variance_explained = 1 - (var(abs(fval - y))/var(y));             

%%
% fitted data
fitx = linspace(min(s),max(s),100);
fity = nan(1,100);
for i = 1:100
    fity(i) = funcval(params, fitx(i), xsign, func);
end
data.raw = [x; y; n];
data.func = func;
data.method = method;
data.fitx = fitx;
data.fity = fity;
data.params = params;
data.variance_explained = variance_explained;

%%
% visualization
if fig==1
    close all;
    figure;      
    plot(fitx, fity, '-k')
    hold on;
    plot(x, y,'ok')
    set(gca,'box','off'); set(gca,'TickDir','out')
    axis square
end

%%
% subfunctions
function v = funcval(p, x, xsign, func)
switch func
    case 'Weibull'
        v = cumWeibull(p, x, xsign);
    case 'Gaussian'
        v = cumGauss(p, x, xsign);
end

function cw = cumWeibull(p, x, xsign)
if xsign==0
    cw = 0.5 + (0.5 - p(3))*(1 - exp(-(x/p(1))^p(2)));
elseif xsign==1
%     cw = p(1) + (1 - 2*p(3))*(1/(1 + exp(-(p(1) + p(2)*x))));
    cw = p(3) + (1 - 2*p(3))*logisticfunc(p(1) + p(2)*x);
end

function y = logisticfunc(x)
y = 1/(1 + exp(-x));

function cg = cumGauss(p, x, xsign)
if xsign==0
    cg = 0.5*(1 + erf((x - p(1))/(sqrt(2)*p(2))));
elseif xsign==1
    cg = normcdf(x, p(1), p(2));
end

function f = likelihood(p, x, y, xsign, func)
v = funcval(p, x, xsign, func);
f = (v^y)*(1 - v)^(1 - y);

function c = cost(p, x, y, n, xsign, func, method)
c = 0;
switch method
    case 'MLE'
        for i = 1:length(n)
            f = @(p)likelihood(p, x(i), y(i), xsign, func);
            if f(p) > 0
                c = c - n(i)*log(f(p));
            end
        end
    case 'LSE'
        for i = 1:length(n)
            f = @(p)funcval(p, x(i), xsign, func);
            c = c + n(i)*sqrt((y(i) - f(p))^2);
        end
end
