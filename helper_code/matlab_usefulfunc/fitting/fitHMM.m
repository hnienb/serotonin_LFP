function hmm_estimate = fitHMM(seq, n_state, analysis_channels)
% fit the Hidden Markov Model (HMM) to sequence of (presumably) spike counts
% in order to estimate the HMM parameters
%
% INPUT: seq ... matrix of spike counts (trials x spike counts x channels) 
%        n_state ... the number of component (2, in default)
%        analysis_channels ... channel number used for parameter estimation
% 
% OUTPUT: hmm_estimate ... output structure
%         n: the number of states, transition: transition matrix
%         emission: emission matrix, fr: firing rate in each state
%         likelihood: model likelihood, likelystates: estimated states
%         duration: duration of each state, variance_explained: variance
%         explained, err_...: confidence intervals of the HMM parameters
%         (this is not implemented...)
%
% % NOTE: - The algorithm works well if seq is binary (0 or 1) as the
%           emission matrix can be smaller (less paramters)
%         - To transform the dimension of likelystates back, use
%         "reshape(likelystates, length(spike counts), length(trials))'"
%
% written by Katsuhisa (29.09.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% inputs
if nargin < 2; n_state = 2; end
if nargin < 3; analysis_channels = 1:size(seq, 3); end

% spike counts across channels
ntr = size(seq, 1);
sc = zeros(ntr, size(seq, 2));
nc = length(analysis_channels);
for n = 1:nc
    sc = sc + seq(:,:,analysis_channels(n));
end

% avoid 0, just for matlab
sc = sc + 1;

% 2-fold cross-validation
alltr = 1:ntr;
idx1 = datasample(alltr, round(ntr/2), 'Replace', false);
idx2 = alltr(~ismember(alltr, idx1));
sc1 = sc(idx1, :);
sc2 = sc(idx2, :);

% 1d sequence
sc1 = sc1';
sc1 = sc1(:)';
sc2 = sc2';
sc2 = sc2(:)';
sccv = {sc1, sc2};
sc = sc';
sc = sc(:)';

% fit HMM 10 times to select sets of parameters yielding the largest likelihood 
% to avoid local optima
repeat = 5;
li = zeros(1, 2*repeat);
ttr_temp = cell(1, 2*repeat); emt_temp = cell(1, 2*repeat);
c = 1;
for r = 1:repeat    
    % cross-validation    
    for t = 1:2
        % initial guess
        [tr_guess, em_guess] = params_initializer(sccv{t}, n_state);
    
        % train HMM
        [ttr_temp{c}, emt_temp{c}] = hmmtrain(sccv{t}, tr_guess, em_guess, ...
            'Algorithm', 'BaumWelch', 'Tolerance', 1e-06, 'Maxiterations', 10000);

        % posterior probability (test)
        [~, ~, frw] = hmmdecode(sccv{-t+3}, ttr_temp{c}, emt_temp{c});
        li(c) = mean(abs(frw(1,:) - 0.5)) + 0.5;
        
        c = c + 1;
    end
end
[li, maxidx] = max(li);
ttr = ttr_temp{maxidx};
emt = emt_temp{maxidx};

% estimate the states
[~, ~, frw] = hmmdecode(sc, ttr, emt);
likelystates = hmmviterbi(sc, ttr, emt);

% fr & duration of each state, variance explained
frvec = sc - 1;
for n = 1:n_state
    hmm_estimate.state(n).fr = mean(sc(likelystates==n) - 1);
    frvec(likelystates==n) = hmm_estimate.state(n).fr;
    hmm_estimate.state(n).duration = [];
end
l = 1;
for i = 2:length(likelystates)
    if likelystates(i)==likelystates(i-1)
        l = l + 1;
    else
        hmm_estimate.state(likelystates(i-1)).duration = ...
            [hmm_estimate.state(likelystates(i-1)).duration, l];
        l = 1;
    end
end
 
% variance explained
sc = sc - 1;
var_tot = sum((sc - mean(sc)).^2);
var_res = abs(var_tot - sum((sc - frvec).^2));
hmm_estimate.variance_explained = 1 - (var_res/var_tot);
 
% structurize
hmm_estimate.n_state = n_state;
hmm_estimate.trnsMat = ttr;
hmm_estimate.emtMat = emt;
hmm_estimate.cvscore = li;
hmm_estimate.likelihood = mean(abs(frw(1,:) - 0.5)) + 0.5;
hmm_estimate.processed_seq = sc;
hmm_estimate.likelystates = likelystates;

function r = drchrnd(a,n)
% random sampling from dirichlet distribution
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

function [tr_guess, em_guess] = params_initializer(sc, n_state)
% initialize parameters
uni = 1:max(sc);
lenuni = length(uni);
tr_guess = zeros(n_state, n_state);
em_guess = zeros(n_state, lenuni);
for n = 1:n_state
    % transition matrix initialized by dirichlet distribution
    tr_guess(n, :) = drchrnd(100*ones(1,n_state)/n_state, 1);
    
    % emission matrix initialized by uniform distribution
    v = rand(1, lenuni);
    em_guess(n, :) = v/sum(v);
end