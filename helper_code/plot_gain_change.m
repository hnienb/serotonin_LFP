function plot_gain_change(analysis)

% Determine which analysis steps to conduct, default is all of them
% Options:

% All - runs all the sections
% Loading - Only loads in the experiment data for each session
% Subspace - Only performs RC Subspace mapping for each session
% Fitting - Only performs linear regression L2 fitting for each session
% Plotting - Only plots the gain change and additive change data for each session
% Stats - Only performs statistical tests on the data for each session


if nargin < 1; analysis = {'all'}; end

path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

root_path = strjoin(parts(1:end-1), '/');
load_path = [strjoin(parts(1:end-2), '/'), '/resources/Data/LFPprepro/'];
save_path = [strjoin(parts(1:end-2), '/'), '/resources/Data/others/'];

if ~exist(save_path, 'dir'); mkdir(save_path); end

%% Load the ex data from the raw data files
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'Loading'))==1

    Lfps_ = load([strjoin(parts(1:end-2), '/'), '/resources/Data/Lfps_pair/Lfps_rc.mat']);
    Lfps = Lfps_.Lfps;
    
    % Search for only files that contain RC data and load those
    parfor file_index = 1:length(Lfps.lfplist)
       
        ex_base_ = load([load_path, Lfps.lfplist{file_index}{1}]);
        ex_base = ex_base_.ex;

        ex_drug_ = load([load_path, Lfps.lfplist{file_index}{2}]);
        ex_drug = ex_drug_.ex;

        % Add base and drug ex info to the results structure
        ex_results{file_index}.ex_base = ex_base;
        ex_results{file_index}.file_base = Lfps.lfplist{file_index}{1};
        ex_results{file_index}.ex_drug = ex_drug;
        ex_results{file_index}.file_drug = Lfps.lfplist{file_index}{2};
    end

    % Save the resulting structure
    save([save_path, '/ex_results.mat'], 'ex_results', '-v7.3');
    disp('experiment data results saved')
    
else % No flag
    load([save_path, '/ex_results.mat']);
    disp('experiment data results loaded')
end


%% Run each ex file pair through RCsubspace.m

if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'Subspace'))==1
    
    cd(root_path)
    helper_path = strjoin([root_path, "/HN Fcts/"], "");
    addpath(genpath(helper_path))

    % Iterate through the resulting ex structs from above
    subspace_results = cell(length(ex_results), 1);

    for ex_index = 1:length(ex_results)
        session = ex_results{ex_index};

        % Copy over the session names for both conditions
        subspace_results{ex_index}.file_base = session.file_base;
        subspace_results{ex_index}.file_drug = session.file_drug;

        fprintf('Running %s through RCsubspace.m ...\n', session.file_base);
        subspace_results{ex_index}.base_RC_res = RCsubspace(session.ex_base, 'lat_flag', 0, 'boot_flag', 0);

        fprintf('Running %s through RCsubspace.m ...\n', session.file_drug);
        subspace_results{ex_index}.drug_RC_res = RCsubspace(session.ex_drug, 'lat_flag', 0, 'boot_flag', 0);
    end

    % Save the resulting structure
    save([save_path, '/subspace_results.mat'], 'subspace_results', '-v7.3');
    disp('RC subspace maps saved')
else
    load([save_path, '/subspace_results.mat']);
    disp('RC subspace maps loaded')
end

%% Run each RC subspace through linear regression fit

if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'Fitting'))==1

    % Iterate through the resulting subspace structs from above
    fit_results = cell(length(subspace_results), 1);

    for subspace_index = 1:length(subspace_results)
        session = subspace_results{subspace_index};

        % Copy over the session names for both conditions
        fit_results{subspace_index}.file_base = session.file_base;
        fit_results{subspace_index}.file_drug = session.file_drug;

        % Compute the regression fit
        x = session.base_RC_res.netSpikesPerFrame;
        y = session.drug_RC_res.netSpikesPerFrame;
        [fit_results{subspace_index}.b0, fit_results{subspace_index}.b1] = ...
            fit_bothsubj2error(x, y);

        % Determine if the unit achieves criterion R^2 (explained variance)
        f = fit_results{subspace_index}.b0 + fit_results{subspace_index}.b1 * x;
        % Find coefficient of correlation and square it to get R^2
        r = corrcoef([y,f]);
        R2 = r(2,1).^2;

        if R2 >= 0.68
            fit_results{subspace_index}.goodR2 = 1;
        else
            fit_results{subspace_index}.goodR2 = 0;
        end
        fit_results{subspace_index}.R2 = R2;

        % Normalize the additive change to the peak response of the baseline
        % condition of the unit
        fit_results{subspace_index}.b0 = fit_results{subspace_index}.b0 / max(session.base_RC_res.netSpikesPerFrame);

        % Determine the drug condition, 5HT or NaCl
        if contains(session.file_drug, '5HT')
            fit_results{subspace_index}.is5HT = 1;
        else
            fit_results{subspace_index}.is5HT = 0;
        end

        % Determine which monkey performed the session
        if contains(session.file_base, 'ka')
            fit_results{subspace_index}.isKaki = 1;
        else
            fit_results{subspace_index}.isKaki = 0;
        end
    end
    
    % Save the resulting structure
    save([save_path, '/fit_results.mat'], 'fit_results', '-v7.3');
    disp('linear regression fit results saved')
else
    load([save_path, '/fit_results.mat']);
    disp('linear regression fit results loaded')
end


%% Plot the results

if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'Plotting'))==1
    
    % Convert the cell array from above to struct array
    s_fit_results = [fit_results{:}];
    s_fit_results = s_fit_results([s_fit_results.goodR2]' == 1);

    % Obtain the drug condition labels for plotting
    drug_labels = [s_fit_results.is5HT]';
    % Obtain the monkey name for plotting
    monkey_id = [s_fit_results.isKaki]';

    % Colors for the drug conditions (red for 5HT, black for NaCl)
    c = ['k', 'r'];
    % Marker shapes for the monkey (circle for Mango, square for Kaki)
    s = ['o', 's'];

    figure;

    % Histogram for relative gain change (x-axis)
    ax1 = subplot(5, 5, [1,4]);
    plotHist([s_fit_results.b1], logical(drug_labels), ax1, 'log');
    set(gca, 'XScale', 'log', 'XTick', [], 'YTick', [0, 12]);
    ax1.LineWidth = 1;
    title('');

    % Colors of the markers
    ax1.Children(3).Color = 'r';
    ax1.Children(4).Color = 'r';
    ax1.Children(4).MarkerFaceColor = 'r';

    % Fix positions of text markers indicating the median for drug condition
    max_val = max([ax1.Children(5).Values ax1.Children(6).Values]);
    ax1.Children(2).YData = max_val + 1.5;
    ax1.Children(4).YData = max_val + 1.5;

    ax1.Children(1).HorizontalAlignment = 'center';
    ax1.Children(3).HorizontalAlignment = 'center';

    ax1.Children(1).Position(1) = ax1.Children(2).XData;
    ax1.Children(1).Position(2) = ax1.Children(2).YData + 2.5;
    ax1.Children(3).Position(1) = ax1.Children(4).XData;
    ax1.Children(3).Position(2) = ax1.Children(4).YData + 2.5;

    % Fix the total height of the plot to fully show markers
    ax1.YLim(2) = ax1.Children(1).Position(2);
    
    % Manually adjust the xlimits so that they will eventually match the
    % main scatter plot's x-axis 
    xlim([0.125, 8]);

    % Histogram for the additive change (y-axis)
    ax3 = subplot(5, 5, [10, 25]);
    plotHist([s_fit_results.b0], logical(drug_labels), ax3);
    set(gca, 'XLim', [-0.6 0.6], 'XTick', [], 'YLim', [0 max(ylim)], 'YTick', [16]);
    ax3.View = [90 -90]; % Rotate plot to align with the scatterplot
    ax3.LineWidth = 1;
    title('');

    % Colors and orientation of the markers
    ax3.Children(3).Color = 'r';
    ax3.Children(4).Color = 'r';
    ax3.Children(4).MarkerFaceColor = 'r';

    ax3.Children(2).Marker = '<';
    ax3.Children(4).Marker = '<';
    ax3.Children(1).Rotation = -90;
    ax3.Children(3).Rotation = -90;

    % Fix positions of text markers indicating the medians for drug condition
    max_val = max([ax3.Children(5).Values ax3.Children(6).Values]);
    ax3.Children(2).YData = max_val + 1.2;
    ax3.Children(4).YData = max_val + 1.2;

    ax3.Children(1).HorizontalAlignment = 'center';
    ax3.Children(3).HorizontalAlignment = 'center';

    ax3.Children(1).Position(1) = ax3.Children(2).XData + 0.05;
    ax3.Children(1).Position(2) = ax3.Children(2).YData + 1.5;
    ax3.Children(3).Position(1) = ax3.Children(4).XData - 0.05;
    ax3.Children(3).Position(2) = ax3.Children(4).YData + 1.5;

    % Fix the total height of the plot to fully show markers
    ax3.YLim(2) = 16;

    % Scatterplot
    ax2 = subplot(5, 5, [6, 24]);
    ax2.LineWidth = 1;
    % Plot gray (unity) lines
    hold on;
    plot([1 1], [-0.6 0.6], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    plot([0.125 8], [0 0], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

    % Iterate through the drug conditions, NaCl or 5HT
    for cond = 0:1
        % Iterate through the monkeys, Mango or Kaki
        for monkey = 0:1
            % Subset the data according to condition and monkey
            x_vals = [s_fit_results(drug_labels==cond & monkey_id==monkey).b1]';
            y_vals = [s_fit_results(drug_labels==cond & monkey_id==monkey).b0]';

            % Plot the subsetted data
            scatter(x_vals, y_vals, [], c(cond+1), s(monkey+1), 'filled', ...
                'markerfacealpha', 0.4, 'markeredgealpha', 0.6, 'markeredgecolor', c(cond+1));
        end

        % Add the number of units in each drug condition
        text(0.13, 0.50 + cond*0.05, ['n = ' num2str(sum(drug_labels==cond))], ...
            'fontsize', 9, 'color', c(cond+1));
    end


    % Format the axes 
    set(gca, 'XScale', 'log', 'XTick', [0.125, 0.25, 1, 4, 8], 'YTick', [-0.6, 0, 0.6]);
    xlim([0.125, 8]);
    ylim([-0.6, 0.6])
    xlabel('relative gain', 'fontsize', 12)
    ylabel('normalized additive change', 'fontsize', 12)

    % Set the figure background to white
    set(gcf, 'color', 'w');  
end

%% Statistical Tests - Wilcoxon Rank Sum

if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'Stats'))==1
    
    % Compute statistics for each individual monkey and combined data

    % Kaki Individual
    kaki_NaCl_b0 = [s_fit_results(drug_labels==0 & monkey_id==1).b0]';
    kaki_NaCl_b1 = [s_fit_results(drug_labels==0 & monkey_id==1).b1]';
    kaki_5HT_b0 = [s_fit_results(drug_labels==1 & monkey_id==1).b0]';
    kaki_5HT_b1 = [s_fit_results(drug_labels==1 & monkey_id==1).b1]';

    pval_kaki_additive = ranksum(kaki_5HT_b0, kaki_NaCl_b0);
    pval_kaki_gain = ranksum(kaki_5HT_b1, kaki_NaCl_b1);


    % Mango Individual
    mango_NaCl_b0 = [s_fit_results(drug_labels==0 & monkey_id==0).b0]';
    mango_NaCl_b1 = [s_fit_results(drug_labels==0 & monkey_id==0).b1]';
    mango_5HT_b0 = [s_fit_results(drug_labels==1 & monkey_id==0).b0]';
    mango_5HT_b1 = [s_fit_results(drug_labels==1 & monkey_id==0).b1]';

    pval_mango_additive = ranksum(mango_5HT_b0, mango_NaCl_b0);
    pval_mango_gain = ranksum(mango_5HT_b1, mango_NaCl_b1);


    % Both monkeys
    all_NaCl_b0 = [s_fit_results(drug_labels==0).b0]';
    all_NaCl_b1 = [s_fit_results(drug_labels==0).b1]';
    all_5HT_b0 = [s_fit_results(drug_labels==1).b0]';
    all_5HT_b1 = [s_fit_results(drug_labels==1).b1]';

    % NaCl vs 5HT
    pval_all_additive = ranksum(all_5HT_b0, all_NaCl_b0);
    pval_all_gain = ranksum(all_5HT_b1, all_NaCl_b1);
   

    % Print out the results to command line
    fprintf('<strong>================ Kaki Stats ================</strong>\n');
    fprintf('Relative gain: p = %.5f, n_(NaCl) = %d, n_(5HT) = %d\n', pval_kaki_gain, size(kaki_NaCl_b1, 1), size(kaki_5HT_b1, 1));
    fprintf('Additive change: p = %.5f, n_(NaCl) = %d, n_(5HT) = %d\n\n', pval_kaki_additive, size(kaki_NaCl_b0, 1), size(kaki_5HT_b0, 1));

    fprintf('<strong>================ Mango Stats ================</strong>\n');
    fprintf('Relative gain: p = %.5f, n_(NaCl) = %d, n_(5HT) = %d\n', pval_mango_gain, size(mango_NaCl_b1, 1), size(mango_5HT_b1, 1));
    fprintf('Additive change: p = %.5f, n_(NaCl) = %d, n_(5HT) = %d\n\n', pval_mango_additive, size(mango_NaCl_b0, 1), size(mango_5HT_b0, 1));

    fprintf('<strong>================ Both Monkey Stats ================</strong>\n');
    fprintf('Relative gain: p = %.5f, n_(NaCl) = %d, n_(5HT) = %d\n', pval_all_gain, size(all_NaCl_b1, 1), size(all_5HT_b1, 1));
    fprintf('Additive change: p = %.5f, n_(NaCl) = %d, n_(5HT) = %d\n\n', pval_all_additive, size(all_NaCl_b0, 1), size(all_5HT_b0, 1));
end