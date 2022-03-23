function mm = mem2fire(varargin)
%% 
% membrane potential-based model, similar to Seillier et al., 2017
%

%% 
% stimulus type
mm.stm.stmtype = 'or';

% stimulus sequence
mm.stm.ntr = 100;
mm.stm.nframe = 1; 
mm.stm.stm_resolution = 10; % ms

% filter
mm.filter.me = 57; % ms
mm.filter.sd = 6;

% threshold for firing
mm.mpara.c = 15;
mm.mpara.v_thre = 0.05;

% visualize
plot_flag = 0;

j = 1;              
while j <= length(varargin)
    switch varargin{j}
        case 'stmtype'
            mm.stm.stmtype = varargin{j+1};            
        case 'ntr'
            mm.stm.ntr = varargin{j+1};
        case 'nframe'
            mm.stm.nframe = varargin{j+1};
        case 'stm_resolution'
            mm.stm.stm_resolution = varargin{j+1};
        case 'filter_me' 
            mm.filter.me = varargin{j+1};
        case 'filter_sd'
            mm.filter.sd = varargin{j+1};
        case 'c'
            mm.mpara.c = varargin{j+1};       
        case 'v_thre'
            mm.mpara.v_thre = varargin{j+1};    
        case 'plot'
            plot_flag = varargin{j+1};
    end
    j = j + 2;
end

% stimulus type and tuning in membrane potential
switch mm.stm.stmtype
    case 'or'
        mm.vm.SD = 40;
        mm.vm.AMP = 20;
        mm.stm.vals = 22.5:22.5:180;
        mm.vm.v_tu = normpdf(mm.stm.vals, 90, mm.vm.SD);
        mm.vm.v_tu = mm.vm.AMP*mm.vm.v_tu/max(mm.vm.v_tu);
    case 'co'
        mm.vm.N = 2;
        mm.vm.AMP = 20;
        mm.vm.c50 = 15;
        mm.vm.OFFSET = 3;
        mm.stm.vals = linspace(1.56, 100, 7);
        mm.vm.v_tu = mm.vm.AMP*((mm.stm.vals.^mm.vm.N)./...
            (mm.vm.c50^mm.vm.N + mm.stm.vals.^mm.vm.N)) + mm.vm.OFFSET;
end

%%
% stimulus sequence
stms = datasample(mm.stm.vals, mm.stm.ntr*mm.stm.nframe, 'Replace', true)';
mm.stm.stm = zeros(mm.stm.ntr, mm.stm.nframe*mm.stm.stm_resolution);
begin = [1, 1];
for i = 1:mm.stm.nframe
   mm.stm.stm(:, begin(1):begin(1)+mm.stm.stm_resolution-1)= ...
        repmat(stms(begin(2):begin(2) + mm.stm.ntr-1), 1, mm.stm.stm_resolution);
   begin(1) = begin(1) + mm.stm.stm_resolution;
   begin(2) = begin(2) + mm.stm.ntr;
end

% stimulus-driven membrane potential
unistm = unique(mm.stm.vals);
lenuni = length(unistm);
mm.vm.vm = mm.stm.stm;
for i = 1:lenuni
    mm.vm.vm(mm.stm.stm==mm.stm.vals(i)) = mm.vm.v_tu(i);
end

% 1D Gaussian temporal filter
g_x = fspecial('gaussian',[1 10*mm.filter.sd], mm.filter.sd);
nhalf = round(length(g_x)/2);
g_x = [zeros(1, mm.filter.me - nhalf), g_x];

% fluctuating membrane potential
mm.vm.vmf = mm.vm.vm;
for i = 1:mm.stm.ntr
    % membrane potential
    mm.vm.vmf(i,:) = imfilter(mm.vm.vm(i,:), g_x);
end

% spike rates
nv = length(mm.mpara.v_thre);
mm.res.fr = zeros(mm.stm.ntr*nv, mm.stm.nframe*mm.stm.stm_resolution);
begin = 1;
for i = 1:nv
    mm.res.fr(begin:begin+mm.stm.ntr-1, :) = ...
        mm.mpara.c*(mm.vm.vmf - mm.mpara.v_thre(i));
    begin = begin + mm.stm.ntr;
end
mm.res.fr(mm.res.fr <= 0) = 0;

% firing
mm.res.spk = arrayfun(@(x) poissrnd(x), mm.res.fr);

% power
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
range = {[0.5,4], [4, 7], [8, 13], [14, 29], [30, 80]};
lenb = length(bands);
mm.vm.power = zeros(mm.stm.ntr*nv, lenb);
for i = 1:nv
    for b = 1:lenb
        for j = 1:mm.stm.ntr
            mm.vm.power(j + (i-1)*mm.stm.ntr, b) = ...
                bandpower(mm.vm.vmf(j,:) - mm.mpara.v_thre(i), 1000, range{b});
        end
    end
end


%%
% visualization
if ~ismember(plot_flag, 0)
    close all;
    figure;

    % stimulus sequence
    subplot(2,4,1)
    imagesc(mm.stm.stm)
    colorbar('southoutside')
    title('stimulus sequence')

    % stimulus tuning
    subplot(2,4,2)
    plot(mm.stm.vals, mm.vm.v_tu, '-ok')
    title('Vm tuning')

    % temporal filter
    subplot(2,4,3)
    plot(g_x)
    title('temporal filter')

    % fluctuating membrane potential
    subplot(2,4,4)
    imagesc(mm.vm.vmf)
    colorbar('southoutside')
    title('membrane potentital')

    % spike rates
    subplot(2,4,5)
    imagesc(mm.res.fr)
    colorbar('southoutside')
    title('spike rates')

    % firing
    subplot(2,4,6)
    imagesc(mm.res.spk)
    colorbar('southoutside')
    title('firing')

    % tuning curve of firing rate
    subplot(2,4,7)
    cols = lines(nv);
    for j = 1:nv
        tu = zeros(lenuni, 2);
        for i = 1:lenuni
            idx = find(mm.stm.stm(:,1)==unistm(i)) + mm.stm.ntr*(j-1);
            tu(i,1) = mean(mean(mm.res.spk(idx, :), 2));
            v = mm.res.spk(idx, :);
            tu(i,2) = std(v(:));
        end
        hold on;
        errorbar(mm.stm.vals, tu(:,1), tu(:,2), '-o', 'color', cols(j,:),'capsize', 0)
    end
    title('spike tuning')
    
    % power
    subplot(2,4,8)
    pref = mm.stm.vals(mm.vm.v_tu==max(mm.vm.v_tu));
    for j = 1:nv
        idx = find(mm.stm.stm(:,1)==pref) + mm.stm.ntr*(j-1);
        me = mean(mm.vm.power(idx, :), 1);
        sd = std(mm.vm.power(idx, :), [], 1);
        hold on;
        errorbar(1:lenb, me, sd, '-o', 'color', cols(j,:),'capsize', 0)
    end
    set(gca, 'XTick', 1:lenb, 'XTickLabel', bands)
    xtickangle(45)
    title('power at preffered stm.')
end