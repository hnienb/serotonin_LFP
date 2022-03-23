function variance_explained = varexp(yorig, ypred)
% % compute variance explained (R-square, coefficient of determination)

SS_tot = sum((yorig - mean(yorig)).^2);
SS_res = sum((yorig - ypred).^2); 
variance_explained = 1 - (SS_res/SS_tot);