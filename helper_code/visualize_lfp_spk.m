function visualize_lfp_spk(ex, col, varargin)
%%
% display spike train and lfp together in an example trial
%

path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-2), '/');
data_path = [dir_path, '/resources/Data/LFPprepro/filtered/'];
addpath(data_path)

if nargin < 1
    load([data_path, 'ka_0258_c1_sortLH_all.grating.ORxRC.mat']);
end
if nargin < 2; col = 'k'; end

% example trial number
i = 20;

% display
figure;
lfpt = [0.8 2];
trange = ex.Trials(i).LFP_prepro_time >= lfpt(1)...
    & ex.Trials(i).LFP_prepro_time <= lfpt(2);
plot(ex.Trials(i).LFP_prepro_time(trange), ex.Trials(i).LFP_prepro(trange), '-', 'color', col)
yy = get(gca, 'YLim');
nspk = length(ex.Trials(i).Spikes);
for n = 1:nspk
    if ex.Trials(i).Spikes(n) >= lfpt(1)
        [~, idx] = min(abs(ex.Trials(i).LFP_prepro_time - ex.Trials(i).Spikes(n)));
        hold on;
        plot(ex.Trials(i).LFP_prepro_time(idx)*[1 1], [1.1*yy(2)+0.1*(yy(2)-yy(1)) yy(2)], '-', 'color', col)
    end
end
axis off
set(gcf, 'position', [1244         849         341         103])