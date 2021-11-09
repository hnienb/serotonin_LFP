function [FitPa,fval,exitflag,fitx,fity] = FitPower(x,y, ...
     varargin)
 % hn 25/08/03

ploton = 0    ;     
EXPONENT = 1;
SLOPE = 2;
guess (1:2) =NaN;

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
  elseif strcmpi ('EXPONENT', str)
    j = j+1;
    guess (EXPONENT) = varargin{j};    
    elseif strcmpi ('SLOPE', str)
        j = j+1;
        guess (SLOPE) = varargin{j};
   
  end
  j = j+1;
end
if isempty (x) | size(x) ~= size(y) 
    
    disp( ' to break')
    
    FitPa.exponent = [];
    FitPa.slope = [];
    
    exitflag = [];
    fval = [];
    type = [];

    fitx = [];
    fity = [];
    
    return;
end


if isnan (guess(SLOPE)) 
    
    guess(SLOPE) = 1;
end

if  isnan (guess(EXPONENT))
    guess(EXPONENT) = 1;
   
end


options = optimset('MaxFunEvals',100000,'maxiter',10000);
LB = [ 0 0 ];
exitflag = 0;
fval = NaN;


[fitparams,fval,exitflag,output]= fminsearch(@powerlaw,...
    guess,options,x,y,LB); %,lower_bound,upper_bound,options)
FitPa.exponent = fitparams(1);
    
FitPa.slope = fitparams(2);


fitx = [];
fity = [];

function f=powerlaw(X0,x,y,LB)
f=[];
funct = X0(2)*x.^X0(1);
f=sum((funct-y).^2);
if X0(1)<1
    f=inf;
end
if X0(2)> 3.5
    f=inf;
end


