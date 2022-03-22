function label_seq = label4StmPerTr(ex)
% assigns label of where the trial is positioned in the stimulus sequence:
% 0; fixation breaks
% 1, 2, 3...; the stimulus presentation number per trial
% written by Katsuhisa (09.08.18), modified by Aashay (03.15.22)

[~,idx] = sort(arrayfun(@(x) x.times_fpOn, ex.Trials));
oTrials = ex.Trials(idx);
rewards = [oTrials.Reward];
label_seq = rewards;









