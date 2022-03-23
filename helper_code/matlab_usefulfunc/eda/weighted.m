function [average,sem,ci] = weighted(values,weights,stds, varargin)
%% compute weighted average and weighted SEM
% INPUT: values, weights, stds ... same length of vectors but values can be
%                            a matrix
% OUTPUT: weighted average, weighted SEM and weighted 95% CI
%
% written by Katsuhisa (14.02.17)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% % replace NaN with the average of its column
% for i = 1:matsize(1)
%       rep = find(isnan(values(i,:)));
%       for n = 1:length(rep)
%             values(i,rep(n)) = nanmean(values(:,rep(n)));
%       end
% end

if nargin < 1
        error('At least a matrix has to be given as an input argument.')
end

% transpose if necessary
if nargin>1
        if size(weights,1)==1
                weights = weights';
        end
end
if size(values,1)==1
        values = values';
end
if nargin==3
        if size(stds,1)==1
                stds = stds';
        end
end

% % remove nan
% okrow = find(all(~isnan(values),2));
% values = values(okrow,:);
% if nargin>1
%         weights = weights(okrow);
% end
% if nargin==3
%         stds = stds(okrow,:);
% end

% replace nan in weights with their median
if nargin > 1        
        weights(isnan(weights)) = nanmedian(weights);   

        % matrix for weights
        weimat = nan(size(values));
        lenwei = nan(1,size(values,2));
        for c = 1:size(values,2)        
                okpos = find(~isnan(values(:,c)));
                lenwei(c) = length(okpos);
                if lenwei(c) < size(values,1)
                        new_weight = nan(size(values,1),1);
                        new_weight(okpos) = weights(okpos)/sum(weights(okpos));
                else
                        new_weight = weights/sum(weights);
                end
                weimat(:,c) = new_weight;
        end
end

% nancol = find(all(isnan(values),1));
% for n = 1:length(nancol)
%         v = values(:,n);
%         v(isnan(v)) = nanmedian(v);
%         values(:,n) = v;
% end
% if nargin>1
%         weights(isnan(weights)) = nanmedian(weights);
% end

len_w = size(values,1);

% compute weighted average
switch nargin
        case 1  % standard error of the mean without weighting
                average = mean(values,1);
                sem = std(values,[],1)/sqrt(size(values,1));
                
        case 2  % compute weighted standard error of the mean
                average = nansum(values.*weimat,1);               
                sem = sqrt(nansum(weimat.*(values - ones(size(values,1),1)*average).^2,1)./lenwei);            
                
%         case 2 % bootstrap
%                 mat = zeros(size(values));
%                 for i = 1:len_w
%                         mat(i,:) = values(i,:).*weights(i);
%                 end
%                 average = sum(mat,1)/sum(weights);                
%                 
%                 mat = mat/sum(weights);
%                 repeats = 500;
%                 aves = zeros(repeats, size(values,2));
%                 for r = 1:repeats
%                         res = zeros(size(values));
%                         for k = 1:size(values,1)
%                                 shu = randperm(size(values,1));
%                                 res(k,:) = mat(shu(1),:);
%                         end
%                         aves(r,:) = mean(res,1);
%                 end
%                 sem = std(aves,[],1);               
                
        case 3 % pooled sem (I don't think this is true because this is used for the mean difference in 2 samples!)
                mat = zeros(size(values));
                for i = 1:len_w
                        mat(i,:) = values(i,:).*weights(i);
                end
                average = sum(mat,1)/sum(weights);
                
                variance = zeros(1,size(values,2));
                winv = 0;
                for i = 1:len_w
                        variance = variance + (weights(i)-1)*stds(i,:).^2;
                        winv = winv + 1/weights(i);
                end
                variance = variance/(sum(weights)-len_w);
                sem = sqrt(variance*winv);        
end

% compute weighted 95% confidence interval of the mean
ci = sem*tinv(0.975, len_w-1); 



