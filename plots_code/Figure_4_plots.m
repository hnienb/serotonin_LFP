% This creates the subplots for Figure 4.
close all

% Get the 'relative' folder path to get the resource folderpath
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

root_path = strjoin(parts(1:end-1), '/');

% For data
dir_path = strjoin(parts(1:end-2), '/');
data_path = [dir_path, '/resources/Data/'];
addpath(genpath([dir_path, '/external_libraries/']))
addpath([root_path, '/plot_functions']) % Plotting helper functions


%% Create the subplot folder if needed, save figures there
save_path = [strjoin(parts(1:end-2), '/'), '/resources/Figure_4/subplots/'];

if ~exist(save_path, 'dir') % Check if folder exists
    mkdir(save_path)
end
addpath(save_path)


%% Plot spike-triggered mixture model correlation

% Plots each condition's (including FR control) spike-triggered mixture model correlation
disp('Statistics for spike-triggered mixture model:')
plot_c2s();

% Remove p-val significance marker from the plot where necessary
fig = gcf;
delete(fig.Children(1).Children(2))
delete(fig.Children(1).Children(1))
delete(fig.Children(4).Children(2))
delete(fig.Children(4).Children(1))

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/corr_mi.fig'])

%% Plot example LFP, predicted spike train and original spike train

% Plots spike-triggered mixture model's predictions for spike train given
% the example LFP and also plots the original spike train (for baseline and
% drug)
c2s_example()

% Get the figure handles
figHandles = findobj('Type', 'figure');

base_fig = figHandles(2);
drug_fig = figHandles(1);

% Save the figures
savefig(base_fig, [save_path, '/ka173_example200_base.fig'])
savefig(drug_fig, [save_path, '/ka173_example200_drug.fig'])

%% Plot an example output of the spike triggered mixture model

% Plots an example output of the spike triggered mixture model
stm_example()

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/stm_example_ka278.fig'])
