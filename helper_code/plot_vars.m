function plot_vars(anaT, varname1, varname2)
%%
% scatter for a pair of variables
%

% Get the 'relative' folder path to get helper functions
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

addpath(genpath([strjoin(parts(1:end-3), '/'), '/external_libraries/matlab_usefulfunc/']))

figure;

% vals
X = anaT.table(:, contains(anaT.varnames, varname1));
Y = anaT.table(:, contains(anaT.varnames, varname2));

mat = [X, Y, anaT.lists(:,1), anaT.lists(:,2)];
ndrug = 1;
drugnames = {'FR'};
if strcmp(anaT.pairnames{1,1}, 'base')
    ndrug = 2;
    drugnames = {'NaCl', '5HT'};
else
    mat(:,3) = 0;
end
for d = 1:ndrug
    subplot(1,ndrug,d)
    if strcmp(varname1(1:end-4), varname2(1:end-4))
        unity_scatter(mat(mat(:,3)==d-1, 1), mat(mat(:,3)==d-1, 2), ...
            mat(mat(:,3)==d-1, 4))
    else
        ls_scatter(mat(mat(:,3)==d-1, 1), mat(mat(:,3)==d-1, 2), ...
            mat(mat(:,3)==d-1, 4))
    end
    title(drugnames{d})
    if length(drugnames) > 1
        xlabel(varname1)
        ylabel(varname2)
    else
        xlabel([varname1(1:end-4) ' (high FR)'])
        ylabel([varname2(1:end-4) ' (low FR)'])
    end
end
set(gcf, 'position', [282   340   560   225])
set(gcf, 'Name', [varname1 ' vs ' varname2], 'NumberTitle', 'off')

% stats    
if ndrug > 1
    kaki_nacl = mat(mat(:,3)==0 & mat(:,4)==0, 1:2);
    kaki_5ht = mat(mat(:,3)==1 & mat(:,4)==0, 1:2);
    mango_nacl = mat(mat(:,3)==0 & mat(:,4)==1, 1:2);
    mango_5ht = mat(mat(:,3)==1 & mat(:,4)==1, 1:2);
    stats_kaki = pair_tests(kaki_nacl, kaki_5ht);
    disp(['Kaki (NaCl; n=' num2str(size(kaki_nacl, 1)) ') ------------------------'])
    disp(stats_kaki.pair(1).table)
    disp(['Kaki (5HT; n=' num2str(size(kaki_5ht, 1)) ') ------------------------'])
    disp(stats_kaki.pair(2).table)
    disp(['Kaki (NaCl vs 5HT; n=' num2str(size(kaki_5ht, 1)+size(kaki_nacl, 1)) ') ---------'])
    disp(stats_kaki.table)
    stats_mango = pair_tests(mango_nacl, mango_5ht);
    disp(['Mango (NaCl; n=' num2str(size(mango_nacl, 1)) ') ------------------------'])
    disp(stats_mango.pair(1).table)
    disp(['Mango (5HT; n=' num2str(size(mango_5ht, 1)) ') ------------------------'])
    disp(stats_mango.pair(2).table)
    disp(['Mango (NaCl vs 5HT; n=' num2str(size(mango_5ht, 1)+size(mango_nacl, 1)) ') ---------'])
    disp(stats_mango.table)
    stats_all = pair_tests([kaki_nacl; mango_nacl], [kaki_5ht; mango_5ht]);
    disp(['Both (NaCl; n=' num2str(size(mango_nacl, 1)+size(kaki_nacl, 1)) ') ------------------------'])
    disp(stats_all.pair(1).table)
    disp(['Both (5HT; n=' num2str(size(mango_5ht, 1)+size(kaki_5ht, 1)) ') ------------------------'])
    disp(stats_all.pair(2).table)
    disp(['Both (NaCl vs 5HT; n=' num2str(size(mango_nacl, 1)+size(mango_5ht, 1)+size(kaki_nacl, 1)+size(kaki_5ht, 1)) ') ---------'])
    disp(stats_all.table)
    
else % Firing Rate Control
    kaki_high = mat(mat(:,3) == 0 & mat(:,4) == 0, 1);
    kaki_low = mat(mat(:,3) == 0 & mat(:,4) == 0, 2);
    mango_high = mat(mat(:,3) == 0 & mat(:,4) == 1, 1);
    mango_low = mat(mat(:,3) == 0 & mat(:,4) == 1, 2);
    
    stats_kaki = pair_tests([kaki_high, kaki_low]);
    disp(['Kaki (n=' num2str(size(kaki_high, 1)) ') ------------------------'])
    disp(stats_kaki.pair(1).table)
    
    stats_mango = pair_tests([mango_high, mango_low]);
    disp(['Mango (n=' num2str(size(mango_high, 1)) ') ------------------------'])
    disp(stats_mango.pair(1).table)
    
    stats_all = pair_tests([[kaki_high;mango_high], [kaki_low; mango_low]]);
    disp(['Both (n=' num2str(size([kaki_high;mango_high], 1)) ') ------------------------'])
    disp(stats_all.pair(1).table)
    
end