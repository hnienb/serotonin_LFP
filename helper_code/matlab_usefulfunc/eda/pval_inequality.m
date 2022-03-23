function p = pval_inequality(p)
% simply convert raw p-val to p-val with inequality

if p >= 0.05
    p = round(1000*p)/1000;
elseif p < 0.05 && p >= 0.01
    p = 0.05;
elseif p < 0.01 && p >= 1e-3
    p = 0.01;
else
    for n = 3:9
        if p < 1/(10^n) && p >= 1/(10^(n+1))
            p = 1/(10^n);
            break
        end
    end
end