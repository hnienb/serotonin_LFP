function [ex0, ex2, sidx1, sidx2] = ex_spliter(ex, es, thre, wnd)
%%
% split ex-file based on firing rate such that the difference between the
% new ex-files is close to the specified effect size 'es'
%

ntr = length(ex.Trials);
nhalf = round(ntr/2);

% get spike counts
[~, spkc] = getSpks(ex.Trials, wnd);

% sort
[~, sortidx] = sort(spkc);

sidx1 = sortidx(1:nhalf);
sidx2 = sortidx(nhalf+1:end);
fr1 = mean(spkc(sidx1));
fr2 = mean(spkc(sidx2));
if fr1*fr2 > 0
    counts = [0, 0, 0, 0];
    disp(['before ...' num2str(fr1/fr2)])
    while fr1/fr2 - es >= thre
        counts(1) = counts(1) + 1;
        if fr1/fr2 - es > 0
            if mod(counts(1), 2)==1
                counts(2) = counts(2) + 1;
                sidx1 = sidx1(1:end-1);
            else
                counts(3) = counts(3) + 1;
                sidx2 = sidx2(1+1:end);
            end
        else
            counts(4) = counts(4) + 1;
            ch1 = sidx1(end - counts(4) + 1);
            ch2 = sidx2(counts(4));
            sidx1(end - counts(4) + 1) = ch2;
            sidx2(counts(4)) = ch1;
        end
        fr1 = mean(spkc(sidx1));
        fr2 = mean(spkc(sidx2));
    end
    disp(['after ...' num2str(fr1/fr2)])
    
end

% median split
ex2 = ex;
ex2.Trials = ex.Trials(sidx1); % Low FR
ex0 = ex;
ex0.Trials = ex.Trials(sidx2); % High FR

