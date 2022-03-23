% This creates the subplots for Figure 1.
close all

% Get the 'relative' folder path to get the resource folderpath
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

root_path = strjoin(parts(1:end-2), '/');

% For data
data_path = [root_path, '/resources/Data/'];

addpath(genpath([root_path, '/helper_code'])) % Plotting helper functions

%% Create the subplot folder if needed, save figures there
save_path = [strjoin(parts(1:end-2), '/'), '/resources/Figure_1/subplots/'];

if ~exist(save_path, 'dir') % Check if folder exists
    mkdir(save_path)
end
addpath(save_path)

%% Plot example baseline LFP and single unit recording
% Plots session ka_258 baseline data, trial 20
visualize_lfp_spk()

% Give the figure window a name to distinguish it
set(gcf, 'name', 'ka_258_20_base');

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path '/ka_258_20_base.fig'])

%% Plot example drug LFP and its single unit recording

% Load the data
load([data_path, '/LFPprepro/filtered/ka_0258_c1_sortLH_17.33.grating.ORxRC_5HT.mat']);
        
% Plots session ka_258 5HT data, trial 20
visualize_lfp_spk(ex, 'r');

% Give the figure window a name to distinguish it
set(gcf, 'name', 'ka_258_20_drug');

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path '/ka_258_20_5HT.fig'])

%% Plot the mean firing rate plot

% Load the LFP data
load([data_path, '/Lfps_pair/Lfps_rc.mat']);

% Process the LFP data into a usable analysis table
anaT = analysis_table(Lfps, 'drug');

% Plot the mean firing rates of baseline condition against drug condition
plot_vars(anaT, 'fr base', 'fr drug');

% Remove p-val significance marker from the plot
fig = gcf;
delete(fig.Children(1).Children(2))
delete(fig.Children(1).Children(1))

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path '/fr.fig'])


%% Plot the Gain Change

% Load, process, and analysis the data to plot gain change
plot_gain_change();

% Give the figure window a name to distinguish it
set(gcf, 'name', 'gain_change');

figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path '/gain_change.fig']);