
function [a,f0] = f1f0Sine (y,time,tf)
% hn 15/08/03
% calculates f1 and f0 via sin- and cos- component.
% corresponds to fitting a Sine-wave at frequency=tf to x.
% units: time and tf must have the same units in time


s = (sin(2*pi*tf*time).*y)/length(y);
c = (cos(2*pi*tf*time).*y)/length(y);
  
a = 2*sqrt((sum(c))^2+(sum(s))^2);          % factor 2 because I only calculate f1 for tf, not for -tf 
f0 = mean(y);

