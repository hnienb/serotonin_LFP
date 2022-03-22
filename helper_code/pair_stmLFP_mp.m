function para = pair_stmLFP_mp(ex0, ex2, exn0, exn2, lfpfield, thin, analysis, z)
%%
% LFP analysis by stimulus type with a pair of sessions
% INPUT: ex0, ex2 ... ex file after 'loadCluster.m'
%        exn0, exn2 ... exn file after 'MP_single.m'    
%        lfpfield ... 'LFP_prepro'
%        thin ... 0, no thining, 1, thinning (firing rate correction)
%        analysis ... 'all', 'sta'
%

% inputs =================
if nargin < 5; lfpfield = 'LFP_prepro'; end
if nargin < 6; thin = 1; end
if nargin < 7; analysis = {'all'}; end
if nargin < 8; z = 1; end

% stimulus type ==============
[stimparam, vals] = getStimParam(ex0);
lenv = length(vals);
if lenv==1 && strcmp(stimparam, 'or')
    para.stm.param = 'rc';
else
    para.stm.param = stimparam;
end
para.stm.vals = vals;

% drug & animal ==================
try
    para.ismango = contains(ex2.fileName, 'ma');
    para.is5ht = contains(ex2.fileName, '5HT');
catch
    para.ismango = nan;
    para.is5ht = nan;
end
    
% preset parameters ==================
fs = 1000;
para.ts = ex0.Trials(end).LFP_prepro_time;
para.window = {[0.8 2]};
stmdur = 2;
para.wnd = 0.1;

% LFP & spikes as a function of stimulus
exs = {ex0, ex2};
exns = {exn0, exn2};
for d = 1:2
    % initialization
    para.cond(d).lfpfull = cell(1, lenv);
    para.cond(d).spk = cell(3, lenv);
    para.cond(d).spkc = cell(3, lenv);
    para.cond(d).ntr = zeros(1, lenv);
    para.cond(d).trials = cell(1, lenv);
    para.cond(d).spk_tu = cell(1, 3); 
    
    % stimulus types
    for i = 1:lenv    
        % trials with the unique stimulus
        para.cond(d).trials{i} = exs{d}.Trials([exs{d}.Trials.(stimparam)] == vals(i));
        mptrs = exns{d}.Trials([exs{d}.Trials.(stimparam)] == vals(i));
        para.cond(d).ntr(i) = length(para.cond(d).trials{i});        

        % firing rate per trial in analysis window
        [para.cond(d).spk{1, i}, para.cond(d).spkc{1, i}]  = ...
            getSpks(para.cond(d).trials{i}, [para.window{end}(1), ...
            stmdur - para.window{end}(2)]);    

        % firing rate per trial during stimulus presentation
        [para.cond(d).spk{2, i}, para.cond(d).spkc{2, i}]  = getSpks(para.cond(d).trials{i}, [0 0]);

        % elicited spikes (mean, SD, n)
        dur = para.window{end}(2) - para.window{end}(1);
        para.cond(d).spk_tu{1}(i, :) = [mean(para.cond(d).spkc{1, i})/dur, ...
            std(para.cond(d).spkc{1, i})/dur, sum(para.cond(d).spkc{1, i})];        
        para.cond(d).spk_tu{2}(i, :) = [mean(para.cond(d).spkc{2, i})/stmdur, ...
            std(para.cond(d).spkc{2, i})/stmdur, sum(para.cond(d).spkc{2, i})];    
        
        % swap preprocessed LFP
        for n = 1:para.cond(d).ntr(i)
            % replace with MPs
            para.cond(d).trials{i}(n).(lfpfield) = mptrs(n).signal;
            para.cond(d).trials{i}(n).energy = mptrs(n).energy;
            para.cond(d).lfpfull{i}(n, 1:length(para.cond(d).trials{i}(n).(lfpfield))) = ...
                para.cond(d).trials{i}(n).(lfpfield);
        end
    end
end

% zscoring
if z==1
    lfpall = [];
    for d = 1:2
        for i = 1:lenv    
            for n = 1:para.cond(d).ntr(i)
                lfpall = [lfpall, para.cond(d).trials{i}(n).(lfpfield)];
            end
        end
    end
    me = nanmean(lfpall);
    sd = nanstd(lfpall);
else
    me = 0; sd = 1;
end

for d = 1:2
    lenlfp = length(para.cond(d).trials{end}(end).(lfpfield));
    for i = 1:lenv    
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).trials{i}(n).(lfpfield) = (...
                para.cond(d).trials{i}(n).(lfpfield) - me(1))/sd(1);
            para.cond(d).lfpfull{i}(n, 1:lenlfp) = para.cond(d).trials{i}(n).(lfpfield);
        end
    end
end

% 'thinning' to control the difference in firing rate 
j = 1;
if thin==1
    j = 3;
    a = [1, 2];
    for i = 1:lenv
        % rate ratio
        alpha = para.cond(2).spk_tu{1}(i, 1)/para.cond(1).spk_tu{1}(i, 1);
        if alpha == 1
            para.cond(1).spk_tu{3}(i, :) = para.cond(1).spk_tu{1}(i, :);
            para.cond(2).spk_tu{3}(i, :) = para.cond(2).spk_tu{1}(i, :);
            para.cond(1).spk{3, i} = para.cond(1).spk{1, i};
            para.cond(1).spkc{3, i} = para.cond(1).spkc{1, i};
            para.cond(2).spk{3, i} = para.cond(2).spk{1, i};
            para.cond(2).spkc{3, i} = para.cond(2).spkc{1, i};
            continue
        elseif alpha < 1
            d = 1;
        elseif alpha > 1
            d = 2;
            alpha = 1/alpha;
        end
        para.cond(d).spk{3, i} = para.cond(d).spk{1, i};
        para.cond(d).spkc{3, i} = para.cond(d).spkc{1, i};
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).spk{3, i}{n} = thinning(para.cond(d).spk{1, i}{n}, alpha);
            para.cond(d).spkc{3, i}(n) = length(para.cond(d).spk{3, i}{n});
        end
        b = a(~ismember(a, d));
        para.cond(b).spk{3, i} = para.cond(b).spk{1, i};
        para.cond(b).spkc{3, i} = para.cond(b).spkc{1, i};

        % elicited spikes (mean, SD, n) after thinning
        para.cond(d).spk_tu{3}(i, :) = [mean(para.cond(d).spkc{3, i})/dur, ...
                std(para.cond(d).spkc{3, i})/dur, sum(para.cond(d).spkc{3, i})];     
        para.cond(b).spk_tu{3}(i, :) = para.cond(b).spk_tu{1}(i, :);     
    end
end

% Spectrogram =============================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'spectrogram'))
    for d = 1:2
        for i = 1:lenv
            para.cond(d).spectrogram.f{i} = exns{1}.freq;
            para.cond(d).spectrogram.t{i} = exns{1}.MPtime;
            for m = 1:para.cond(d).ntr(i)
                if m==1
                    para.cond(d).spectrogram.S{i} = para.cond(d).trials{i}(m).energy;
                else
                    if sum(isnan(para.cond(d).trials{i}(m).energy(:)))==0
                        para.cond(d).spectrogram.S{i} = para.cond(d).spectrogram.S{i} + para.cond(d).trials{i}(m).energy;
                    end
                end
            end
            para.cond(d).spectrogram.S{i} = para.cond(d).spectrogram.S{i}/para.cond(d).ntr(i);
        end   
    end
end


% Spike-triggered LFP =====================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'sta'))
    for d = 1:2
        for i = 1:lenv
            sta = [];
            for n = 1:para.cond(d).ntr(i)
                % compute a standard STA
                stlfp = getSTA(para.cond(d).trials{i}(n).(lfpfield), ...
                    para.cond(d).trials{i}(n).LFP_prepro_time, ...
                    para.cond(d).spk{j, i}{n}, para.wnd, fs);        
                   
                if ~isnan(stlfp)
                    % stack stLFP
                    sta = [sta; stlfp];  
                end
            end
            if ~isempty(sta)
                % STA as a function of stimulus type
                para.cond(d).sta.mean(i,:) = nanmean(sta, 1);
                para.cond(d).sta.sd(i,:) = nanstd(sta, [], 1);
                para.cond(d).sta.nspk(i) = size(sta, 1);

                % mean STA spectrogram
                [~, para.cond(d).sta.f{i}, ~, p] = ...
                    spectrogram_frange(para.cond(d).sta.mean(i,:), 95, fs, [0 100]);
                para.cond(d).sta.p{i} = nanmean(p, 3);
                para.cond(d).sta.t{i} = linspace(-para.wnd, para.wnd, length(para.cond(d).sta.mean(1,:)));
            else
                para.cond(d).sta.mean(i,:) = nan(1, length(-para.wnd:1/fs:para.wnd));
                para.cond(d).sta.sd(i,:) = para.cond(d).sta.mean(i,:);
                para.cond(d).sta.nspk(i) = 0;
                para.cond(d).sta.f{i} = nan;
                para.cond(d).sta.t{i} = nan;
                para.cond(d).sta.p{i} = nan;
            end
        end  
    end     
end


% Spike-LFP coherency ========================
if ismember(1, contains(analysis, 'all')) || sum(contains(analysis, 'coherence'))
    ts = 0.3; % 3 tapers
    params = define_params(fs, ts, 0);
    for d = 1:2
        for i = 1:lenv
            % coherence analysis by chronux toolbox
            [~, ~, para.cond(d).coherence.S12{i}, para.cond(d).coherence.S1{i}, ...
                para.cond(d).coherence.S2{i}, para.cond(d).coherence.f{i}] = ...
                spike_lfp_coherence(para.cond(d).lfpfull{i}, para.cond(d).trials{i}(end).LFP_prepro_time, ...
                para.cond(d).spk{j, i}, para.window{end}, params);
            
            % trial average to deal with nans
            para.cond(d).coherence.S12{i} = squeeze(nanmean(para.cond(d).coherence.S12{i}, 2));
            para.cond(d).coherence.S1{i} = squeeze(nanmean(para.cond(d).coherence.S1{i}, 2));
            para.cond(d).coherence.S2{i} = squeeze(nanmean(para.cond(d).coherence.S2{i}, 2));
            para.cond(d).coherence.f{i} = nanmean(para.cond(d).coherence.f{i}, 1)';
            
            % trial average coherence and phase
            C12 = para.cond(d).coherence.S12{i}./sqrt(para.cond(d).coherence.S1{i}.*para.cond(d).coherence.S2{i});
            para.cond(d).coherence.C{i} = abs(C12);
            para.cond(d).coherence.phi{i} = angle(C12);
        end   
    end
end

% remove fields for memory =====================
rms = {'spk', 'spkc', 'trials', 'lfpfull'};
for i = 1:length(rms)
    para.cond = rmfield(para.cond, rms{i});
end
para = rmfield(para, 'ts');


%% subfunction =========================================
function [C, phi, S12, S1, S2, f] = spike_lfp_coherence(lfpv, lfpt, spk, wnd, params)
% spike-LFP coherence analysis

ntr = length(spk);
lfp = lfpv(:, lfpt >= wnd(1) & lfpt <= wnd(2));
s = struct('times', []);

for n = 1:ntr
    % adjust spike timings
    s(n).times = spk{n} - wnd(1);
end
% get coherence by chronux function     
[C, phi, S12, S1, S2, f] = coherencycpt(lfp', s, params, 1); 



function spk = thinning(spk, fr_ratio)
% perform thinning of spikes
% INPUT: spk ... spike times (1D)
%              fr_ratio ... ratio of mean firing rate across two conditions
%                              (u_small/u_large)
%              wnd ... window to apply thinning
%

if isempty(spk)
    return
end

r = rand(1, length(spk));
spk(r <= 1 - fr_ratio) = []; 



function [nsc, nov, nff] = stff_params(L, overlap)
nsc = floor(L/18);
nov = floor(nsc/(100/overlap));
nff = 1024;

function [s, f, t, p] = spectrogram_frange(v, overlap, fs, frange)
[nsc, nov, nff] = stff_params(length(v), overlap);
[s, f, t, p] = spectrogram(v, nsc, nov, nff, fs);
outrange = f >= frange(1) & f <= frange(2);
s = s(outrange, :);
p = p(outrange, :);
f = f(outrange);
