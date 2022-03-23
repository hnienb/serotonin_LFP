function dm = dimensionality_reduction(datatable, method)
%%
% perform dimensionality reduction and plot 2d scatter
%
% INPUT:
% datatable ... matrix of data (observations x variables)
% method ... 'pca' or 'tsne' 
%
% OUTPUT:
% dm ... structure containing information for dimensionality reduction
%
% EXAMPLE: dm = dimensioanlity_reduction(rand(100, 6), 'pca');
%

% inputs
if nargin < 2; method = 'pca'; end

% dimensionality reduction
dm.method = method;
switch dm.method
    case 'pca'
        [dm.coef, dm.Y, dm.latent, dm.tsquared, dm.explained]...
            = pca(datatable);
        labs = {'PC1', 'PC2'};
    case 'tsne'
        dm.Y = tsne(datatable);
        labs = {'tsne dim 1', 'tsne dim 2'};
end

% visualization
scatter(dm.Y(:,1), dm.Y(:,2), 30, 'markerfacecolor', 'k', 'markeredgecolor', 'w')
xlabel(labs{1})
ylabel(labs{2})
set(gca, 'box', 'off', 'tickdir', 'out')