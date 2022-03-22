function params = define_params(fs, ts, trave)
%%
% define parameter sets for chronux toolbox
% INPUT:      fs ... sampling frequency (Hz)
%             ts ... sampling duration (sec)
%             trave ... trial average 1 or 0
%

if nargin < 1; fs = 1000; end
if nargin < 2; ts = 0.2; end
if nargin < 3; trave = 1; end

W = 10/fs;
N = floor(ts*fs); 
K = 2*N*W - 1; % the number of Slepians

params.tapers = [W*N, K]; % was [2, 3]. Chalk et al.,2010 used [3,5]
params.err = 0;
params.Fs = fs;
params.fpass = [0 100];
params.pad = 0; % nearest 2^N
params.trialave = trave; % if 0, the file is too heavy