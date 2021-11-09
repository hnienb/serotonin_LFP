function [FitPa,fval,exitflag,type,fitx,fity] = FitGauss_HN(x,y, ...
     varargin)
 % hn 25/08/03
 % fits gaussian to data
 % Options:
 %  'both: Gauss in lin or log- space, whichever gives better fit (default)
 %  'lin': Gauss in lin space
 %  'log':  "     " log  "
 %  'plot': calculates datapoints to plot fitted Gaussian & returns them in: fitx, fity
 
 % returns :
 % fitparams(1) = Amplitude
 %  "       (2) = SD
 %  "       (3) = Mean
 %  "       (4) = y-offset(base)   ... og Gauss
 % fval      residual (sq) at fitparams
 % exitflag
 % type     'lin' or 'log' - space
 % fitx     (optional) dadapoints to plot fit
 % fity
basefix = 'off';
type = 'both';
ploton = 0    ;          
AMP = 1;
SD = 2;
MEAN = 3;
BASE = 4;
guess (1:4) =NaN;
LB = [ 0 0 0 0];
nvar = nargin - 2;
j = 1;
while j  <= nvar
  str = varargin{j};
  if strcmpi('log',str)
    type = varargin{j};    
  elseif strcmpi('lin',str)
    type = varargin{j};
  elseif strcmpi('plot',str)
    ploton = 1;
  elseif strcmpi ('amp', str)
    j = j+1;
    guess (AMP) = varargin{j};    
    elseif strcmpi ('sd', str)
        j = j+1;
        guess (SD) = varargin{j};
    elseif strcmpi ('mean', str)
        j = j+1;
        guess (MEAN) = varargin{j};
    elseif strcmpi ('base',str)
        j = j+1;
        guess (BASE) = varargin {j};
    elseif strcmpi ('LB',str)
        j = j+1;
        LB = varargin {j};
  elseif strcmpi ('basefix',str)
        j = j+1;
        basefix = 'on';
        fixedbase = varargin{j};
        
  end
  j = j+1;
end

if isempty (x) |( size(x,2) ~= size(y,2) ) | (size(x,1) ~=size(x,1))
    
    disp( ' to break')
    
    FitPa.amp = [];
    FitPa.sd = [];
    FitPa.mean = [];
    FitPa.base = [];
    
    exitflag = [];
    fval = [];
    type = [];

    fitx = [];
    fity = [];
    
    return;
end

if size(x,2) < size(x,1)
    x = x';
    y=y';
end



if isnan (guess(AMP)) 
    guess(AMP) = (max(y) - min(y));
end

if  isnan (guess(MEAN))
    [ maximum maxpos] = max(y);
   
    guess(MEAN) = x(maxpos);         
end

if isnan (guess(SD))
    [minimum minpos] = min(y);
    guess(SD) = abs(x(minpos) -x(maxpos))/3;
    guess(SD) = std(x .* y)/mean(y);
end

if isnan (guess(BASE))
    guess(BASE) = min(y);
end

options = optimset('MaxFunEvals',100000,'maxiter',10000);

exitflag1 = 0;
fval1 = NaN;

if strcmpi(basefix,'on')
     LB(4) = fixedbase;
     guess(1) = guess(1)-fixedbase;
    guess(4) = [];
     if strcmpi('lin',type) | strcmpi('both',type)
        [fitparams1,fval1,exitflag1,output1]= fminsearch(@gauss0,...
            guess,options,x,y,LB); %,lower_bound,upper_bound,options)
    end
    fitparams1(4) = fixedbase;
    
    exitflag2 = 0;
    fval2 = NaN;
    if strcmpi('log',type) | strcmpi('both',type)
        [fitparams2,fval2,exitflag2,output2]= fminsearch(@gausslog0,...
        guess,options,x,y,LB); %fit to gaussian in log-space
    end
    
    fitparams2(4) = fixedbase;
    
else 
    if strcmpi('lin',type) | strcmpi('both',type)
        [fitparams1,fval1,exitflag1,output1]= fminsearch(@gauss,...
            guess,options,x,y,LB); %,lower_bound,upper_bound,options)
    end

    exitflag2 = 0;
    fval2 = NaN;
    if strcmpi('log',type) | strcmpi('both',type)
        [fitparams2,fval2,exitflag2,output2]= fminsearch(@gausslog,...
        guess,options,x,y,LB); %fit to gaussian in log-space
    end

end


if strcmpi('lin',type) |((exitflag1==1) & ((fval1<=fval2) | exitflag2 <=0 ))| ...
     (exitflag1<=0 & exitflag2 <=0 & fval1<=fval2 )   %looks for better fit
    fitparams=fitparams1;
    fval=fval1;
    exitflag=exitflag1;
    output=output1; 
    type='lin';  
else
    fitparams=fitparams2;
    fval=fval2;
    exitflag=exitflag2;
    output=output2; 
    type='log';
end

FitPa.amp = fitparams(1);
FitPa.sd = fitparams(2);
FitPa.mean = fitparams(3);
FitPa.base = fitparams(4);

fitx = [];
fity = [];
if ploton
    fitx = min(x): 0.2 : max(x);
    if strcmpi(type ,'lin')
        fity = fitparams(AMP) * ...
            exp(- (fitx-fitparams(MEAN)).^2/(2*fitparams(SD)^2)) + fitparams(BASE);
    else
        fity = fitparams(AMP) * ...
            exp(- (log(fitx)-fitparams(MEAN)).^2/(2*fitparams(SD)^2)) + fitparams(BASE);
    end
end


