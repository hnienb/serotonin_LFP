function gpr_estimate = fitGPFA(seq, analysis_channels)
% fit the Gaussian Process Factor Analysis to sequence of (presumably) spike counts
% in order to estimate the HMM parameters
% INPUT: seq ... matrix of spike counts (trials x spike counts x channels) 
%        analysis_channels ... channel number used for parameter estimation
%
% OUTPUT: gpr_estimate ... output structure
% 
% written by Katsuhisa (29.09.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% inputs
if nargin < 2; analysis_channels = 1:size(seq, 3); end

% spike counts across channels
sc = zeros(size(seq, 1), size(seq, 2));
nc = length(analysis_channels);
for n = 1:nc
    sc = sc + seq(:,:,analysis_channels(n));
end

% 1d sequence
sc = sc';
sc = sc(:)';

% avoid 0, just for matlab
sc = sc + 1;

% fit GPFA
t = [1:length(sc)]';
gprMdl = fitrgp(t, sc', 'Basis', 'linear', 'KernelFunction', 'squaredexponential');

% loss of prediction
L = loss(gprMdl, t, sc');

% prediction
predseq = resubPredict(gprMdl);

% structurize
gpr_estimate.Mdl = gprMdl;
gpr_estimate.loss = L;
gpr_estimate.processed_seq = sc - 1;
gpr_estimate.pred = predseq';