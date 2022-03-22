function runAlllfp(type)

if nargin < 1; type = 'all'; end

% path =======================================
path = mfilename( 'fullpath' );

if ispc    % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

rootpath = strjoin(parts(1:end-2), '/');

addpath(genpath([rootpath '/helper_code/']))
addpath(genpath([rootpath '/external_libraries/']))

loadpath = [rootpath, '/resources/Data/LFPprepro/'];

    
% LFP filtering =============================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'filter'))   
   
    % load data lists
    a = load([loadpath, '/lfplist.mat']);
    lfplist = a.lfplist;
    N = length(lfplist);
    savepath = [loadpath, 'filtered/'];
    
    if ~exist(savepath, 'file')
        mkdir(savepath)
    end

    parfor i = 1:N
        % baseline 
        try
            loadCluster([loadpath lfplist{i}{1}], 'save', [savepath lfplist{i}{1}]);
        catch
            disp([lfplist{i}{1} ': error in loadCluster for baseline. skip...'])
            continue
        end        
        
        % drug
        try
            loadCluster([loadpath lfplist{i}{2}], 'save', [savepath lfplist{i}{2}]);
        catch
            disp([lfplist{i}{2} ': error in loadCluster for drug. skip...'])
            continue
        end
    end
    disp('all saved!')
end


% MP preprocessing (stimulus response) =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'MP'))
    
    % load data and lists
    a = load([loadpath, '/lfplist.mat']);
    lfplist = a.lfplist;
    N = length(lfplist);
    
    for i = 1:N
        % baseline 
        ex0 = loadCluster([loadpath, lfplist{i}{1}], 'filtering', 0);        
    
        % drug
        ex2 = loadCluster([loadpath, lfplist{i}{2}], 'filtering', 0);
    
        % perform MP & save the ex
        MP_single(ex0, lfplist{i}{1});
        MP_single(ex2, lfplist{i}{2});
    end
end


% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'c2sFormat'))
    
    % Paths
    loadpath_filtered = [loadpath, 'filtered/'];
    savepath = [rootpath, '/resources/Data/c2s/data/'];
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
   
    % load data and lists
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    lfplist = a.lfplist;
    N = length(lfplist);
    
    parfor i = 1:N            
        try
            % load =======================
            d0 = load([loadpath_filtered, lfplist{i}{1}], 'ex');
            d2 = load([loadpath_filtered, lfplist{i}{2}], 'ex');
            
            % analysis
            [~, data] = data4c2s(d0.ex, d2.ex, 'LFP_prepro', 1);

            % save
            c2s_saver(data, lfplist{i}{2}(1:end-4), savepath)
 
            disp(['session ' num2str(i) ' saved!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
end


% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair'))
    % stimulus type
    stmtype = 'rc';
    
    % Paths
    loadpath_filtered = [loadpath, 'filtered/'];
    savepath = [rootpath, '/resources/Data/Lfps_pair/'];
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
   
    % load data and lists
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N);
    is5ht = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N               
        try
            % load =======================
            d0 = load([loadpath_filtered, lfplist{i}{1}], 'ex');
            d2 = load([loadpath_filtered, lfplist{i}{2}], 'ex');
            
            % analysis
            Out1{i} = pair_stmLFP(d0.ex, d2.ex, 'LFP_prepro');

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end

    % unit info
    oks = ~cellfun('isempty', Out1);
    Lfps_all.lfplist = lfplist(oks);
    Lfps_all.animal = animal(oks);
    Lfps_all.is5ht = is5ht(oks);

    % results from analysis    
    Lfps_all.LFP_prepro = Out1(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'animal', 'is5ht'};
        
    for f = 1:length(fields)
        Lfps.(fields{f}) = Lfps_all.(fields{f})(:);
    end
    
    Lfps.LFP_prepro = Lfps_all.LFP_prepro(:);
    
    % autosave
    save([savepath, 'Lfps_' stmtype '.mat'], 'Lfps', '-v7')

    disp([stmtype ': lfps saved!']) 
end

    
% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_nothin_sc'))
    % stimulus types
    stmtype = 'rc';
   
    % Paths
    loadpath_filtered = [loadpath, 'filtered/'];
    savepath = [rootpath, '/resources/Data/Lfps_pair_nothin_sc/'];
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    
    % load data and lists
    a = load([loadpath, 'lfplist.mat'], 'lfplist');  
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    es = 0.6592; % Effect size (approx. ratio between drug and baseline conditions)
    thre = 0.01;        % Threshold to approximate effect size
    wnd = [0.8 0];      % 0.8 sec following stim onset, 0 prior to stim end
    Out1 = cell(1, N);
    is5ht = zeros(1, N); animal = zeros(1, N);

    parfor i = 1:N           
        try
            % load baseline data =======================
            d0 = load([loadpath_filtered, lfplist{i}{1}], 'ex');
            
            % split trials based on spike counts
            [ex0, ex2] = ex_spliter(d0.ex, es, thre, wnd);
            
            % analysis
            Out1{i} = pair_stmLFP(ex0, ex2, 'LFP_prepro', 0);

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end

    % unit info
    outs = cellfun('isempty', Out1);

    lfplist = lfplist(outs==0);
    Lfps_all.lfplist = lfplist;
    Lfps_all.animal = animal(outs==0);
    Lfps_all.is5ht = is5ht(outs==0);

    % results from analysis    
    Lfps_all.LFP_prepro = Out1(outs==0); 
    
    % split by stimulus type
    fields = {'lfplist', 'animal', 'is5ht'};

    % baseline 
    for f = 1:length(fields)
        Lfps.(fields{f}) = Lfps_all.(fields{f})(:);
    end
    
    Lfps.LFP_prepro = Lfps_all.LFP_prepro(:);

    % Rename variables to make it easier to distinguish later down
    Lfps_sc = Lfps;
   
    % autosave
    save([savepath, 'Lfps_' stmtype '_sc.mat'], 'Lfps_sc', '-v7')
    
    disp([stmtype ': sc lfps saved!'])
end


% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp'))
    % stimulus types
    stmtype = 'rc';
    
    % Paths
    loadpath_filtered = [loadpath, 'filtered/'];
    loadpath_mp = [loadpath, 'MP/Trials/'];
    savepath = [rootpath, '/resources/Data/Lfps_pair/'];
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
   
    % Load data and lists
    a = load([loadpath 'lfplist.mat'], 'lfplist');
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); 
    is5ht = zeros(1, N); animal = zeros(1, N);

    parfor i = 1:N            
        try
            % load =======================
            d0 = load([loadpath_filtered, lfplist{i}{1}], 'ex');
            d2 = load([loadpath_filtered, lfplist{i}{2}], 'ex');
            d3 = load([loadpath_mp, lfplist{i}{1}], 'exn');
            d4 = load([loadpath_mp, lfplist{i}{2}], 'exn');
            
            % analysis
            Out1{i} = pair_stmLFP_mp(d0.ex, d2.ex, d3.exn, d4.exn, 'LFP_prepro');

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end

    % unit info
    oks = ~cellfun('isempty', Out1);
    Lfps_all.lfplist = lfplist(oks);
    Lfps_all.animal = animal(oks);
    Lfps_all.is5ht = is5ht(oks);

    % results from analysis    
    Lfps_all.LFP_prepro = Out1(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'animal', 'is5ht'};

    for f = 1:length(fields)
        Lfps.(fields{f}) = Lfps_all.(fields{f})(:);
    end
    Lfps.LFP_prepro = Lfps_all.LFP_prepro(:);
    
    % autosave
    save([savepath, 'Lfps_mp_' stmtype '.mat'], 'Lfps', '-v7')
    disp([stmtype ': MP lfps saved!']) 
end


% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp_nothin_sc'))
    % stimulus types
    stmtype = 'rc';
    
    % Paths
    loadpath_filtered = [loadpath, 'filtered/'];
    loadpath_mp = [loadpath, 'MP/Trials/'];
    savepath = [rootpath, '/resources/Data/Lfps_pair_nothin_sc/'];
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end

    % Load lists
    a = load([loadpath 'lfplist.mat'], 'lfplist');
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    es = 0.6592; % Effect size (approx. ratio between drug and baseline conditions)
    thre = 0.01;        % Threshold to approximate effect size
    wnd = [0.8 0];      % 0.8 sec following stim onset, 0 prior to stim end
    Out1 = cell(1, N); 
    is5ht = zeros(1, N); animal = zeros(1, N);

    parfor i = 1:N            
        try
            % load baseline =======================
            d0 = load([loadpath_filtered, lfplist{i}{1}], 'ex');
            d3 = load([loadpath_mp, lfplist{i}{1}], 'exn');

            % split trials based on sc
            [ex0, ex2, sidx1, sidx2] = ex_spliter(d0.ex, es, thre, wnd);
            exn3_0 = d3.exn;
            exn3_2 = d3.exn;
            exn3_0.Trials = exn3_0.Trials(sidx2);
            exn3_2.Trials = exn3_2.Trials(sidx1);

            % analysis
            Out1{i} = pair_stmLFP_mp(ex0, ex2, exn3_0, exn3_2, 'LFP_prepro', 0);

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end

    % unit info
    outs = cellfun('isempty', Out1);
    lfplist = lfplist(outs==0);
    Lfps_all.lfplist = lfplist;
    Lfps_all.animal = animal(outs==0);
    Lfps_all.is5ht = is5ht(outs==0);

    % results from analysis    
    Lfps_all.LFP_prepro = Out1(outs==0); 
    
    % split by stimulus type
    fields = {'lfplist', 'animal', 'is5ht'};
  
    % baseline 
    for f = 1:length(fields)
        Lfps.(fields{f}) = Lfps_all.(fields{f})(:);
    end   
    Lfps.LFP_prepro = Lfps_all.LFP_prepro(:);
    
    % Rename the variable for later use
    Lfps_sc = Lfps;

    % autosave
    save([savepath, 'Lfps_mp_' stmtype '_sc.mat'], 'Lfps_sc', '-v7')
    disp([stmtype ': MP sc lfps saved!'])
end


function c2s_saver(data, fname, savepath)
if ~exist([savepath, fname], 'dir')
    mkdir([savepath, fname])
end
save([savepath, fname, '/data.mat'], 'data', '-v7.3') % v7.3 due to file size