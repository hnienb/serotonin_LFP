function  [FitPa,fval,exitflag,fitx,fity] = FitRatioOfGaussians(x,y, ...
     varargin)
 
 %hn 9/18/12
 %fits ratio of Gaussian model (based on Cavanaugh, Bair, Movshon, 2002)
 %to size tuning data.
 %R(x) = k_c*L_c(x)/(1+k_s*L_s(x)),
 %where 
 %  k_c: gain of center RF
 %  k_s: gain of surround field
 %  L_c = (2/(pi) * integral over 0 to x of exp(-(y/w_c)^2) dy)^2
 %  w_c: width (SD of Gaussian) of center RF
 %  L_s = (2/(pi) * integral over 0 to x of exp(-(y/w_s)^2) dy)^2
 %  w_s: width (SD of Gaussian) of surround field
 % 
 %  constraint: w_c < w_s
 
 
 FitPa.k_c = [];
 FitPa.k_s = [];
 FitPa.w_c = [];
 FitPa.w_s = [];
 fitx = [];
 fity = [];
 

 %options = optimset('MaxFunEvals',1000000,'maxiter',100000);
 options = optimset('MaxFunEvals',10000,'maxiter',10000);
LB = [ 0 0 0 0 0 0 0 0 0];
exitflag = 0;
fval = NaN;
fit_flag=0;


idx = find(y>=0.98*max(y));
if ~isempty(idx)
    w_c_guess = x(idx(1));
else
    w_c_guess = x(end);
end
%guess = [w_c_guess w_c_guess+1 0.05 0.1] ; % w_c, w_s, k_c, k_s --> starting for mouse (HN)
guess = [w_c_guess w_c_guess+1 0.05 0.1] ; % w_c, w_s, k_c, k_s --> starting for monkey (CL)

[fitparams,fval,exitflag,output]= fminsearch(@MyRatioOfGaussians,...
    guess,options,x,y,LB); %,lower_bound,upper_bound,options)

 FitPa.k_c = fitparams(3);
 FitPa.k_s = fitparams(4);
 FitPa.w_c = fitparams(1);
 FitPa.w_s = fitparams(2); 

 
function f = MyRatioOfGaussians(X0,x,y,LB);

f1 = @(t1) exp(-(t1/X0(1)).^2);
f2 = @(t2) exp(-(t2/X0(2)).^2);
for n=1:length(x)
Q(n) = X0(3)*(2/sqrt(pi)*integral(f1,0,x(n)))^2/...
    (1+(2/sqrt(pi)*X0(4)*integral(f2,0,x(n)))^2);
% Q(n) = X0(3)*(2/sqrt(pi)*integral(f1,0,x(n)))^2/...
%     (1+X0(4)*(2/sqrt(pi)*integral(f2,0,x(n)))^2);
end
f = sum((Q-y).^2);

if (X0(1)>X0(2) | X0(1)<0 |X0(2)<0 |X0(3)<0 |X0(4)<0)
    f=inf;
end

 
 