function res = RCsubspace( ex, varargin )
%RCsubspace 
% 
% 
% In this batch file, I accumulated most of the analysis for the flashed
% grating experiments. This includes the subspace reverse correlation.
%




%% get errorbars by resampling
res = resampleRC(ex, varargin{:});

end




%%
function res = resampleRC(ex, varargin)

rng(9123234);

% result of the 
res = HN_computeLatencyAndNetSpk([], ex, varargin{:});

end