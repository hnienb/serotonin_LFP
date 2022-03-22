function MP_single(ex, filename)
%%
% MP decomposition preprocessing for a single ex-file with preprocessed LFP data
% 
% INPUT: ex ... ex-file after 'loadCluster.m'
%        filename ... file name of ex (e.g.
%        ma_0014_c1_sortLH_5.37PM.grating.ORxRC_5HT.mat)
%
% NOTE: This function assumes that you are in the '../LFP_project/serotonin_code' path
% 

% preset parameters
Fs = 1000; % sampling frequency
folderName = '../resources/Data/LFPprepro/MP/'; % Needs to be relative because of string length limitations

max_iter = 500;
wrap = 1;
L = 4096; % 2^N
if (2/ex.exp.StimPerTrial + 0.2)*Fs < 1024
    L = 1024;
elseif (2/ex.exp.StimPerTrial + 0.2)*Fs < 2048
    L = 2048;
end
signalRange = [1 L];
T = length(ex.Trials(end).LFP_prepro_time);
ntr = length(ex.Trials);
lfpfield = 'LFP_prepro';


% MP decomposition
freq = 0:Fs/L:100;
lenf = length(freq);

inputSignal = nan(L, ntr, 1);
idxs = ones(1, ntr);
for k = 1:ntr
    % padding and mean subtraction for signal processing
    [z, idxs(k)] = padding(ex.Trials(k).(lfpfield), L);
    inputSignal(:, k, 1) = z - nanmean(z);
end

% MP =========================
% perform Gabor decomposition
importData(inputSignal, folderName, filename, signalRange, Fs);
runGabor(folderName, filename, L, max_iter);

% signal reconstruction
exn.MPtime = ex.Trials(end).LFP_prepro_time;
exn.freq = freq;
MPtrials = struct('signal', [], 'energy', []);

% load data
gaborInfo = getGaborData(folderName, filename, 1);
parfor k = 1:ntr
    % 50 Hz line noise and >100 Hz are removed
    cf = gaborInfo{k}.gaborData(2, :); % atom frequency
    atomList = find((cf > 0 & cf < 49.5) | (cf > 50.5 & cf < 150));    
    
    % signal
    mp_signal = reconstructSignalFromAtomsMPP(...
        gaborInfo{k}.gaborData, L, wrap, atomList);
    mp_signal = mp_signal(idxs(k):idxs(k)+T-1);
    MPtrials(k).signal = mp_signal;

    % energy
    rEnergy = reconstructEnergyFromAtomsMPP(...
        gaborInfo{k}.gaborData, L, wrap, atomList);
    rEnergy = rEnergy(1:lenf, idxs(k):idxs(k)+T-1);
    MPtrials(k).energy = rEnergy;
end
exn.Trials = MPtrials;

% autosave =======================
save([folderName 'Trials/' filename], 'exn', '-v7.3')

disp([filename ' saved after MP!'])

%% subfunction
function [b, pre] = padding(a, L)
% add padding for signal processing 
lena = length(a);
b = zeros(1, L);
pre = floor((L - lena)/2) + 1;
me = nanmean(a);
b(1:pre-1) = me;
b(pre:pre+lena-1) = a;
b(pre+lena:end) = me;