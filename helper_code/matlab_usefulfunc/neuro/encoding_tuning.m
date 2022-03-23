function tu = encoding_tuning(stm, res, stmtype)
% compute encoding properties that the tuning curve has
% assuming "rate encoding"
%
% INPUT:
% stm ... stimulus parameters
% res ... neural responses (e.g. firing rate)
% stmtype ... 'or', or 'other'
% 
% OUTPUT:
% - tuning curve (average response to each stimulus)
% - reliability (var average / var all)
% - selectivity (p-val from anova)
% - SNR square (mean / std)^2
% - discriminability (decoding accuracy)
% - metabolic cost of encoding (entropy, conditional entropy, mutual information)
%

ntr = size(stm, 1);
if ~isequal(ntr, size(res, 1))
    error('The size of stm and res must match!')
end
if nargin < 3; stmtype = 'other'; end

if ~strcmp(stmtype, 'or')
    stmtype = 'other';
end

%%
% tuning curve
tu.unistm = unique(stm)';
lenuni = length(tu.unistm);
tu.mean = nan(1, lenuni);
tu.std = nan(1, lenuni);
tu.ntr = nan(1, lenuni);
for u = 1:lenuni
    tu.ntr(u) = sum(stm==tu.unistm(u));
    tu.mean(u) = mean(res(stm==tu.unistm(u)));
    tu.std(u) = std(res(stm==tu.unistm(u)));
end

%%
% reliability
tu.reliability = var(tu.mean)/var(res);

%%
% selectivity
[tu.selectivity, tbl] = anova1(res, stm, 'off');

%%
% SNR2
tu.snr2 = (tu.mean./tu.std).^2;
tu.snr2(isnan(tu.snr2)) = 0;
tu.snr2(isinf(tu.snr2)) = 100;

%% 
% discriminability (eta square --- effect size of anova)
tu.discriminability = tbl{2, 2}/tbl{4,2};
or = tu.unistm;
if strcmp(stmtype, 'or')
    if max(tu.unistm) - min(tu.unistm) > 2*pi
        or = tu.unistm * pi /180;
    end   
    or = mod(or, 2*pi);
end
% acc = zeros(ntr, 1);
% tr = 1:ntr;
% for i = 1:ntr
%     % leave-one-out
%     trs = tr(~ismember(tr, i));
% 
%     % multinomial logistic regression
%     B = mnrfit(res(trs), categ(trs));
%     
%     % model prediction of stimulus type
%     prb = mnrval(B, res(i));
%     acc(i) = ismember(find(tu.unistm==stm(i)), find(prb==max(prb)));
% end
% tu.discriminability = sum(acc)/ntr;

%%
% metabolic cost (entropy, conditional entropy, mutual information)
tu.metabcost = zeros(1, 3);
int_res = round(res);
unires = unique(int_res);
lenr = length(unires);
for r = 1:lenr
    % entropy
    pr = sum(int_res==unires(r))/ntr;
    tu.metabcost(1) = tu.metabcost(1) + pr*log2(1/pr);
    
    % conditional entropy
    for s = 1:lenuni
        ps = sum(stm==tu.unistm(s))/ntr;
        prs = sum(int_res==unires(r) & stm==tu.unistm(s))/sum(stm==tu.unistm(s));
        if prs > 0
            tu.metabcost(2) = tu.metabcost(2) + ps*prs*log2(1/prs);
        end
    end
end
% mutual information
tu.metabcost(3) = tu.metabcost(1) - tu.metabcost(2);

%%
% stimulus specific quantity
switch stmtype
    case 'or'
        % circular variance --- Ringach et al. (2002)
        % compute weighted sum of cos and sin of angles
        mn = tu.mean;
        if sum(mn < 0) > 0
            mn = mn + abs(min(mn));
        end
        r = sum(mn.*exp(1i*or));
        % obtain length 
        r = abs(r)./sum(mn);
        tu.unique.circularvariance = 1 - r;
        
        % direction selectivity
        [rp, irp] = max(tu.mean);
        pdeg = tu.unistm(irp); 
        if pdeg - 180 < 0
            udeg = pdeg + 180;
        else
            udeg = pdeg - 180;
        end
        ru = tu.mean(tu.unistm==udeg);
        tu.unique.directionsel = (rp - ru)/(rp + ru);
    otherwise
        tu.unique = nan;
end
