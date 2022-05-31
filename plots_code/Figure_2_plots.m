% This creates the subplots for Figure 2.
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

addpath(genpath([root_path, '/helper_code/'])) % Plotting helper functions

%% Create the subplot folder if needed, save figures there
save_path = [strjoin(parts(1:end-2), '/'), '/resources/Figure_2/subplots/'];

if ~exist(save_path, 'dir') % Check if folder exists
    mkdir(save_path)
end
addpath(save_path)

%% Load necessary LFP data
load([data_path, '/Lfps_pair/Lfps_rc.mat'])
load([data_path, '/Lfps_pair_nothin_sc/Lfps_rc_sc.mat']);

anaT = analysis_table(Lfps, 'drug');
anaT_sc = analysis_table(Lfps_sc, 'sc');

% Remove any c2 files, since its corresponding LFP is the same as c1
idx = 0;
for i = 1:length(Lfps.lfplist)
    if contains(Lfps.lfplist{i}{1}, 'c2')
        idx = i;
        break
    end
end

anaT.table = anaT.table([1:idx-1,idx+1:end], :);
anaT.lists = anaT.lists([1:idx-1,idx+1:end], :);


%% Plot Gamma Band Scatterplot for drug condition

% Plots scatterplots for the Gamma power between drug condition and baseline
disp('Statistics for Gamma Band Power for drug conditions:')
plot_vars(anaT, 'gamma pow base', 'gamma pow drug')

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Remove p-val significance marker from the plot
fig = gcf;
delete(fig.Children(1).Children(2))
delete(fig.Children(1).Children(1))

% Save the figure
savefig(figHandles(1), [save_path '/gamma_scatter.fig'])

%% Plot Gamma Band Scatterplot for firing rate control

% Plots scatterplots for the Gamma power between low/high firing rate trials
disp('Statistics for Gamma Band Power for firing rate control:')
plot_vars(anaT_sc, 'gamma pow base', 'gamma pow drug')

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path '/gamma_scatter_sc.fig'])


%% Plot Low-Frequency Scatterplot for Drug Condition

% Plots scatterplots for the Gamma power between drug condition and baseline
disp('Statistics for Low-Frequency Band Power for drug conditions:')
plot_vars(anaT, 'low-freq pow base', 'low-freq pow drug')

% Remove p-val significance marker from the plot
fig = gcf;
delete(fig.Children(1).Children(2))
delete(fig.Children(1).Children(1))

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path '/lowfreq_scatter.fig'])


%% Plot Low-Frequency Scatterplot for Firing Rate Control

% Plots scatterplots for the Low-Frequency power between low/high firing rate trials
disp('Statistics for Low-Frequency Band Power for firing rate control:')
plot_vars(anaT_sc, 'low-freq pow base', 'low-freq pow drug')


% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path '/lowfreq_scatter_sc.fig'])

%% Plot the Power Spectra Density and Spectrogram for the Drug Conditions

% Plots the PSD and spectrograms for drug conditions
lfps_pair_visualize(Lfps, 'spectrogram', 'drug')

% Obtain the last two created figure handles
figHandles = findobj('Type', 'figure');

psd_fig = figHandles(2);
spec_fig = figHandles(1);

% Save both figures
savefig(psd_fig, [save_path '/psd.fig'])
savefig(spec_fig, [save_path '/spectrogram.fig'])


%% Plot the Power Spectra Density and Spectrogram for Firing Rate Controls

% Plots the PSD and spectrograms for firing rate controls
lfps_pair_visualize(Lfps_sc, 'spectrogram', 'sc')

% Obtain the last two created figure handles
figHandles = findobj('Type', 'figure');

psd_fig = figHandles(2);
spec_fig = figHandles(1);

% Save both figures
savefig(psd_fig, [save_path '/psd_sc.fig'])
savefig(spec_fig, [save_path '/spectrogram_sc.fig'])
