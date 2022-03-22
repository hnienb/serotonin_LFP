function tcfit = fitOR( mn, sem, or, bootstrp)
%tcfit = fitOR( mn, sem, or )
%
% This function fits a gaussian function to the raw orientation tuning
% curve. 
% The input arguments are vectors:
% - mn:     the spike count or spike rate 
% - sem:    standard error of mn, also a vector
% - or:     stimulus orientation 
% 
% If you add a fourth argument, the tuning fit is bootstrapped and another
% field with these results is added to the output structure.
%
% Fitting is done with restrictions. All parameters (except the preferred
% orientation), have to be positive. To avoid and overestimation of the
% amplitude, this parameters has an upper bound of twice the amplitude of
% the raw tuning curve. Check 'fitHelper' for more information on the
% fitting options.
% 
% 
% The output argument is a struct with the 
%  - fitting results: 
%       parameters mu (preferred OR), sig (sigma, tuning width), a
%       (amplitude), b (offset), the of explained variance of the fit r2,
% - raw tuning curve:
%       val.mn (mean firing rate), val.sem (SEM), val.or (stimulus 
%       samples)
% - well-sampled, fitted gaussian function 
%       x (in units of mn) and y (in degree)
%
%
% 
% @CL

%%% the random number generator seed
%%% do not change 
rng(9123234); 


%%% neglect spontaneous activity in the input arguments
i_noblank = or < 180;
or = or(i_noblank);
mn = mn(i_noblank);
sem = sem(i_noblank);


%%% the gaussian function cannot handle circular data.
%%% thus, center the data around the highest peak and fit the gaussian 
%%% function to it.
pk = or(find(max(mn) == mn, 1, 'first'));
[val.or, val.mn, val.sem] = centerData(or, mn, sem, pk);
[fit_res,~] = fitHelper(pk, val);

%%% repeat the previous two steps to center around the estimated gaussian
%%% mean
[val.or, val.mn, val.sem] = centerData(or, mn, sem, fit_res.mu);
[fit_res,gof2] = fitHelper(fit_res.mu, val);


%%% assign the fitting parameters to a result structure
tcfit.mu = fit_res.mu;        tcfit.sig = fit_res.sig;
tcfit.a = fit_res.a;          tcfit.b = fit_res.b;
tcfit.r2 = gof2.rsquare;


%%% shift the orientations if the preferred orientation is outside 0-180
%%% window
if tcfit.mu > 180 
    tcfit.mu = tcfit.mu -180;
    val.or = val.or-180;
elseif  tcfit.mu < 0
    tcfit.mu = tcfit.mu +180;
    val.or = val.or+180;
end


%%% save the raw and fitted function data, especially helpfull to plot
%%% the results later
tcfit.val = val; % raw tuning curve
tcfit.x = tcfit.mu-100:tcfit.mu+100; % descriptive function, rich of or samples
tcfit.y = gaussian(tcfit.mu, tcfit.sig, tcfit.a, tcfit.b, tcfit.x) ;


%%% bootstrap the fitting results to gain quality estimates of the
%%% parameters
if nargin == 4
    % resample the mean responses and repeat the fitting 
    parfor i = 1:1000
        bootidx(i,:) = randi(length(mn), length(mn), 1);
        boot(i) = fitOR( mn(bootidx(i,:)), sem(bootidx(i,:)), or(bootidx(i,:)));
    end
    tcfit.boot = boot;
end
end


%% Helper functions
function [fit_res,gof] = fitHelper(pk, val)
% fit data to gaussian function using matlab built-in fitting algorithms


%%% define fitting input
amp = max(val.mn)-min(val.mn); % amplitude of the raw tc
x0 = [pk 90 amp 0]; % starting point

fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower', [-inf 0 0 0], 'Upper',[inf inf 2*amp inf],...
        'StartPoint', x0, 'MaxFunEvals', 10^5, 'MaxIter', 10^5);
% fo = fitoptions('Method','NonlinearLeastSquares',...
%     'StartPoint', x0, 'MaxFunEvals', 10^5, 'MaxIter', 10^5); %unrestricted fit

ft = fittype(@(mu, sig, a, b, x) gaussian(mu, sig, a, b, x), 'options', fo);


% fit parameters and derive goodness of fit
[fit_res,gof] = fit(val.or, val.mn, ft);

end

function f = gaussian(mu, sig, a, b, x)
f = (a/(sig*sqrt(2*pi)))*exp(-((x-mu).^2)/(2*sig^2)) + b;
end

function [or, mn, sem] = centerData(or, mn, sem, pk)
%%% center data around peak 
 
or(or > pk+90) = or(or > pk+90) -180;
or(or < pk-90) = or(or < pk-90) +180;

[or, idx] = sort(or);
mn = mn(idx);
sem = sem(idx);

or=or'; sem=sem'; mn=mn';
end