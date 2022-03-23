function [binned] = tcbin(v, binsize)
%%
% bin the vector v (can be matrix) based on the specified binsize.
% note that this binning simply splits 'v' into 'binsize' elements
% from the beginning to the end (different from a 'hist' like method).
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++

lenv = size(v,2);
r = mod(lenv, binsize);
binned = floor(lenv/binsize)*ones(size(v,1), binsize);
i = 1;
while r > 0
    binned(1,i) = binned(1,i) + 1;
    r = r - 1;
    if i <= binsize
        i = i + 1;
    else
        i = 1;
    end
end
begin = 1;
for b = 1:binsize
    step = binned(1,b);
    binned(:,b) = nanmean(v(:,begin:begin+step-1),2);
    begin = begin + step;
end


