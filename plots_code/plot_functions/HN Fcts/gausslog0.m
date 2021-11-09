function f = gausslog(X0,x,y,LB);

f=[];
% Test file - try to fit single gaussian to data
% y = X0(1)*exp(-(x-X0(2))^2/X0(3)) +LB(4)
% compared with gauss, gauss0 has the base fixed at LB(4)

for p = 1:size(y,2)
	funct(p) = X0(1) * exp(-(log(x(p))-X0(3)).^2 /(2*(X0(2).^2)))+LB(4);
    if funct(p)<0 funct(p)=0; end;
    f = [f, funct(p) - y(p)];
    
end
f=f.^2;
f=sum(f);

for n=[1]              % check whether X0 is within the limits
    if X0(n) < LB(n)
        f=inf;
    else f=f;
    end
end;

