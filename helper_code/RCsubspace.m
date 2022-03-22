function [ res, ratestruc, tcfit ] = RCsubspace( ex, varargin )
% RCsubspace 
% 
% 
% In this batch file, I accumulated most of the analysis for the flashed
% grating experiments. This includes the subspace reverse correlation
% and latency analysis but also the tuning fit.
%
% Optional arguments are:
%   'plot'                  - plots the ML latency estimate plot, see
%                               friedmanpriebe.m for more information
%   'lat_flag', boolean     - computes the ML latency estimate (default) or not.
%   'bootstrap_n', integer  - specifies the number of bootstrap samples
%   'boot_flag', boolean    - performs bootstraping and TC fitting (default) or not (AP 2021)
%

%% define variables and parse input
boot_flag = 1;      % Flag for bootstrapping
nsmpl = 1000;       % number of resampling processes  
j=1; 

while j<= length(varargin)
   switch varargin{j}
       case 'exinfo'
           j=j+1;
       case 'bootstrap_n'
           nsmpl = varargin{j+1};
           j=j+1;
       case 'boot_flag'
           boot_flag = varargin{j+1};
           j=j+1;
   end
   j=j+1;
end

%% get errorbars by resampling
[res, ratestruc] = resampleRC(ex, nsmpl, boot_flag, varargin{:});


%% fit tuning curves

if ratestruc ~= []
    tcfit = fitOR( ratestruc.mn', ratestruc.sem', ratestruc.or');


    % this applies only if there is a second stimulus dimension
    % I hardcoded it for orientation and contrast variation
    if size(res.sdfs.mn_rate(1).mn,1)>1 && size(res.sdfs.mn_rate(1).mn,2)>1

        % I am not quite sure hot to best approach the fitting.
        % For now I fit orientation tuning curve for each contrast and 
        % one contrast tc collapsed across all orientations.
        for i = 1:size(res.sdfs.mn_rate(1).mn,2)
            mn     = res.sdfs.mn_rate(1).mn(:,i);
            sem    = res.sdfs.mn_rate(1).sem(:,i);
            or     = res.sdfs.x(:,1);
            tcfit.others.OR(i) = fitOR( mn, sem, or);
        end

        co_spkrates = mean(res.sdfs.mn_rate(1).mn, 1);
        co_samples = res.sdfs.y(1,:);
        if length(res.sdfs.mn_rate)>1
            co_spkrates = [co_spkrates, res.sdfs.mn_rate(2).mn];
            co_samples = [co_samples, 3];
        end
        tcfit.others.CO = fitCO( co_spkrates, co_samples);
    end
else
    tcfit = [];
end



% %% plot SDFs and latency estiamtes
% if p_flag
% 
%     % this is not yet adapted to handle two stimulus dimensions
%     
%     % the raw and fitted spike tuning curve
%     subplot(2,2,1);
%     plot(tcfit.x, tcfit.y, 'r-'); hold on;
%     errorbar(tcfit.val.or, tcfit.val.mn, tcfit.val.sem, 'rx');
%     xlabel('orientation'); ylabel('spks/frame');
%     
%     % stimulus triggered SDFs
%     subplot(2,2,[3 4]);
%     c = hsv(length(res.sdfs.s));
%     for i = 1:length(res.sdfs.s) % orientation samples
%         plot(res.times/10, res.sdfs.s{i}, 'Color', c(i,:), 'DisplayName', ...
%             sprintf('%1.0f \t n=%1.0f', res.sdfs.x(i), res.sdfs.n(i))); 
%         hold on;
%     end
%     plot(res.times/10, res.sdfs.extras{1}.sdf, 'r:'); ho
%     legend('show', 'Location', 'EastOutside');
%     
%     
%     % deviation across selective SDFs and latency estiamtes
%     latfp = res.latFP;  % latency estimate using ML (parametric)
%     lathmax = res.lat;  % latency estimate using half-maximal response threshold (non-p)
%     dur = res.dur;      % response duration (non-parametric only)
%     
%     ylim_ = get(gca, 'ylim'); 
%     
%     plot([latfp latfp], ylim_, 'k');
%     plot([lathmax lathmax], ylim_, 'k--');
%     plot([lathmax+dur, lathmax+dur], ylim_, 'k');
%     
%     
%     s = horzcat(res.sdfs.s{:}); 
%     meanfr = mean(mean(s(201:400),2)); % mean firing rate across all stimulus conditions
%     
%     title(sprintf('lat: %1.1f (hmax %1.1f), dur: %1.1f, \n average sd: %1.2f, mean fr: %1.2f',...
%         res.latFP, res.lat, res.dur, mean(sqrt(res.vars(201:400))), meanfr));
%     xlim([0 160]); grid on;
%     ylabel('spk/s');
% 
% end

end




%%
function [res, ratestruc] = resampleRC(ex, nsmpl, boot_flag, varargin)
% randomly choose n trials and compute netspikes
% do this nsmpl times
% 
% 
% @CL
% 


rng(9123234);

% result of the 
res = HN_computeLatencyAndNetSpk([], ex, varargin{:});

if boot_flag
    % initialize variable to speed up the bootstrapping
    res_boot = cell(nsmpl,1); 


    parfor i = 1:nsmpl

        % get indices from uniform distribution between 1 and number of Trials
        bootidx = randi(length(ex.Trials), length(ex.Trials), 1); 

        % assign Trials
        ex_boot = ex;
        ex_boot.Trials = ex.Trials(bootidx);

        % use HN function to compute the usual res struct
        res_boot{i} = HN_computeLatencyAndNetSpk([], ex_boot, varargin{:});
    %     latfp_boot(i) = res_boot{i}.latFP;
    %     lathmax_boot(i) = res_boot{i}.lat;
        bootsmpl(:,:,i) = res_boot{i}.netSpikesPerFrame;


        if ~isempty(res.netSpikesPerFrameBlank)
            nspk_blank(i) = res_boot{i}.netSpikesPerFrameBlank(1);
        else

        end

    end


    % get the distribution information
    ratestruc = struct();
    % for orientation data
    for i1=1:size(bootsmpl,1) % contrast
        for i2=1:size(bootsmpl,2) % orientation
            res.sdfs.mn_rate(1).mn(i1,i2) = mean(bootsmpl(i1, i2, :));
            res.sdfs.mn_rate(1).sem(i1,i2) = std(bootsmpl(i1, i2, :));
        end
    end
    ratestruc.or= res.sdfs.x(:,1); 
    [~, maxi]= max(res.sdfs.y(1,:));

    ratestruc.mn  = res.netSpikesPerFrame(:, maxi);
    ratestruc.sem = std( bootsmpl, 0, 3  );


    %for blanks
    if ~isempty(res.netSpikesPerFrameBlank)
        res.sdfs.mn_rate(2).mn = mean(nspk_blank);
        res.sdfs.mn_rate(2).sem = std(nspk_blank);

        ratestruc.or= [ratestruc.or; 1000];
        ratestruc.mn = [ratestruc.mn; res.sdfs.mn_rate(2).mn];
        ratestruc.sem = [ratestruc.sem; res.sdfs.mn_rate(2).sem];
    end


    % res.sdfs.lat_boot = latfp_boot;
    % res.sdfs.lathmax_boot = lathmax_boot;

    ratestruc.bootstrap = bootsmpl;
    
else
    ratestruc = [];
end

end