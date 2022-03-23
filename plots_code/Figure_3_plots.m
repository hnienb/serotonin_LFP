% This creates the subplots for Figure 3.
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

addpath([root_path, '/helper_code/']) % Plotting helper functions

%% Create the subplot folder if needed, save figures there
save_path = [strjoin(parts(1:end-2), '/'), '/resources/Figure_3/subplots/'];

if ~exist(save_path, 'dir') % Check if folder exists
    mkdir(save_path)
end
addpath(save_path)

%% Load necessary LFP data
load([data_path, '/Lfps_pair/Lfps_rc.mat'])
load([data_path, '/Lfps_pair_nothin_sc/Lfps_rc_sc.mat']);

anaT = analysis_table(Lfps, 'drug');
anaT_sc = analysis_table(Lfps_sc, 'sc');

%% Plot stLFP sessions (for drug condition)

% Plots each sessions's stLFP waveform (For drug condition)
lfps_pair_visualize(Lfps, 'sta_batch', 'drug');

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/sta_batch.fig'])


%% Plot STA amplitudes (for drug conditions)

% Plots the amplitudes of each session's STA (for drug condition)
disp('Statistics for stLFP amplitudes for drug conditions:')
plot_vars(anaT, 'sta amp base', 'sta amp drug')

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/sta_amp.fig'])


%% Plot STA amplitudes (for firing rate control)

% These plots aren't used in the actual figure 3, but the p-values
% associated with this information are used

% Plots the amplitudes of each session's STA (for firing rate control)
disp('Statistics for stLFP amplitudes for firing rate control:')
plot_vars(anaT_sc, 'sta amp base', 'sta amp drug')

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/sta_amp_sc.fig'])


%% Plot the spike-LFP coherence for each session (For drug conditions)

% Plots the spike-LFP coherence for each session
disp('Statistics for spike-LFP coherence:')
plot_vars(anaT, 'coh low-freq base', 'coh low-freq drug')

% Remove p-val significance marker from the plot
fig = gcf;
delete(fig.Children(1).Children(2))
delete(fig.Children(1).Children(1))

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/coh_scatter.fig'])


%% Plot the scatterplot for delta STA amp and delta LFP power

% Plots the scatterplot between the changes in STA amplitudes and changes in 
% low-freq LFP power for each session (for drug condition)
disp('Statistics for delta STA amplitudes against delta LFP power:')
plot_vars(anaT, 'd sta amp', 'd low-freq pow')

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/sta_lowfreqpow.fig'])


%% Plots the STA for each drug condition and corresponding change in Power Spectra Density

% Plots the multiplot
lfps_pair_visualize(Lfps, 'sta', 'drug')

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/sta_rc.fig'])


%% Plots the STA for firing rate control and corresponding change in Power Spectra Density

% Plots the multiplot
lfps_pair_visualize(Lfps_sc, 'sta', 'sc')

% Obtain the last created figure handle
figHandles = findobj('Type', 'figure');

% Save the figure
savefig(figHandles(1), [save_path, '/sta_sc.fig'])

