function [u,c] = unique_count(x)
%% 
% count the occurence of unique numbers
%
% +++++++++++++++++++++++++++++++++++++++++

u = unique(x);
c = zeros(1, length(u));
for i = 1:length(u)
    c(i) = sum(x==u(i));
end