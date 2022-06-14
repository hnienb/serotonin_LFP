# serotonin_LFP
Custom analysis code for Patel et al. Serotonergic modulation of local network processing in V1 shares signatures with the modulation by spatial attention, 2022

# Contents
## Data
Source data needed to recreate paper figures - located in the /resources/Data/ folder

MP-processed LFP data is not included due file-size constraints, but it can be
replicated using the available code (see below).

#### Format
The spiking information and LFP information are stored in separate files with the following naming convention:

Example:

**ma_0014_c1_sortLH_5.30PM.grating.ORxRC[_drug].mat**

- **ma**: Mango (subject name) ["ka" if Kaki]
- **0014**: Unique session number
- **c1**: Identified single-unit activity
- **5.30PM**: Recording time stamp
- **grating**: Stimulus type (circular sinusoidal luminance grating)
- **ORxRC**: Flashed gratings of random orientations
- **_drug**: Drug condition, if not baseline condition (5HT or NaCl)


Each .mat file includes a struct *ex* with the following information:

- **fileName**: Internally stored file name for the session data
- **exp**: Contains information about the experiment for the session
  - **e1**: First experiment information
    - **type**: Stimulus type (orientation [OR] / contrast [CO])
    - **vals**: Values of the stimulus type used in the session
  - **e2**: (Optional) Second experiment, same format as *e1*
- **Trials**: Includes information about trials (see below)
- **setup.refreshRate**: Refresh rate of the stimulus monitor
- **Header.VisStimVersion**: (Optional) Some sessions used stimulus presentation version
  that requires manual adjustment for times

*ex.Trials* contains the following information:
- **st**: Boolean value if stimulus was presented during trial
- **me**: Monocular eye condition (always set to 0 since only binocular data used)
- **or**: Orientation used for grating stimulus for the trial
- **co**: Contrast used for grating stimulus for the trial
- **TrialEnd**: End time of trial (in sec)
- **TrialStart**: Start time of trial (in sec), i.e. animal starts fixating
- **Start**: Start time of every stimulus video frame [there is typically a small, variable delay between TrialStart (fixation onset) and the first stimulus frame presented]
- **Reward**: Specifies whether the subject response was correct (1) or the trials was discarded (0) [Only trials that are 1 are valid trials]
- **RewardSize**: Size of the reward delivered
- **Spikes**: Spike times (relative to TrialStart), in sec
- **LFP**: Measured LFP voltage for the trial
- **LFP_ts**: Time points corresponding to LFP voltages
- **times_fpOn**: Time (in sec) when fixation point begins presentation
- **st_seq**: Sequence of stimuli presented during the trial
- **or_seq**: Sequence of orientation values presented during the trial
- **co_seq**: Sequence of contrast values presented during the trial



Following preprocessing and filtering, these raw session .mat files are used to
create compiled *Lfps_rc.mat* files. There are corresponding versions for MP-processed
data, spike-count separated data, and combined MP/spike-count data, as well.

Each *Lfps_rc.mat* file includes a struct with the following information:
- **lfp_list**: A list of pairs of baseline and drug condition sessions' filenames
- **animal**: A boolean list containing animal ID (1 if Mango, 0 if Kaki)
- **is5ht**: A boolean list containing drug ID (1 if serotonin, 0 if saline control)
- **LFP_prepro**: A struct list containing the preprocessed LFP data (see below)

*Lfps.LFP_prepro* contains the following information:
- **stm**: Stimulus information, including the modified parameter and its corresponding value
- **ismango**: Boolean to determine animal ID (1 if Mango, 0 if Kaki)
- **is5ht**: Boolean to determine drug ID (1 if serotonin, 0 if saline control)
- **window**: Analysis window used ([0.8, 2] from stimulus onset)
- **wnd**: Window step size used for analysis
- **cond**: 1x2 struct that contains information for both baseline and drug condition
  - **ntr**: Number of trials within the session
  - **spk_tu**: Evoked spike tuning information of format [mean, SD, n]. 1st entry
  corresponds to spikes within analysis window, 2nd entry corresponds to spikes
  within stimulus presentation, 3rd entry corresponds to thinned spike data within
  analysis window
  - **spectrogram**: Spectrogram information (including frequency, time bins)
  - **sta**: Spike-triggered average information, including its corresponding spectrogram
  - **coherence**: Coherency data (including spike-spike, spike-LFP, LFP-LFP)


## Code

### Preprocessing
The central file used for data analysis is *runAlllfp.m*, located in **/analysis_code/**.
This file uses arguments to specify which analysis step to carry out, defaulting to
"all". The following possible values are:
- **"filter"**: Iterates through the session raw files and preprocesses the LFP
  data. **Note**: This step has to be run first for all non-MP and non-c2s data.
- **"MP"**: Iterates through the session raw files and preprocesses the LFP data
  using the [MP algorithm](https://github.com/supratimray/MP). **Note**: This step has to be run first for all MP data and this will expectedly run for multiple hours.
- **"c2sFormat"**: Iterates through the session raw files and organizes the data
  such that the [c2s algorithm](https://github.com/jonasrauber/c2s-docker) (Spike-triggered mixture model) can read it.
- **"lfps_pair"**: Converts the individual preprocessed data into a single struct
  for later analysis.
- **"lfps_pair_nothin_sc"**: Converts the individual preprocessed data into a
  single struct for later analysis, using an algorithm to match the spike rate
  ratio induced by 5HT application.
- **"lfps_pair_mp"**: Converts the MP individual preprocessed data into a single
  struct for later analysis.
- **"lfps_pair_mp_nothin_sc"**: Converts the MP individual preprocessed data into
  a single struct for later analysis, using an algorithm to match the spike rate
  ratio induced by 5HT application.

This code greatly benefits from using the Matlab Parallel Computing Toolbox. Due
to the immense size of the intermediate data files, it is recommended to run this
processing code on a hard drive with at least 170GB and decent read/write speeds.

**Note**: After running "c2sFormat" in *runAlllfp.m*, run *fix_data_dim.m*,
located in /analysis_code. This will split the data into baseline, drug, and FR
control high/low. After this, use a BASH-supported terminal to run the c2s algorithm as follows:
`source run_c2s.sh`

This will additionally require [Docker](https://docs.docker.com/get-docker/) to be installed.

After running *run_c2s.sh*, run *extract_c2s.m* (located in **/analysis_code/**),
which will condense the individual data files into a singular matrix for analysis (*met_cv10.mat*).

### Plotting
Each of the paper's figures are correspondingly generated from the files included
in **/plots_code/** and **/figures_code/**. First run the files in **/plots_code/**
and then run the files in **/figures_code/**.

### External Libraries
This additional libraries will need to be installed for the analysis and plotting
code to run properly. Download these and place them in the **/external_libraries/**
directory (may have to create the directory first in the root folder).

- **[cbrewer](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)** (from Mathworks File Exchange)
- **[MP](https://github.com/supratimray/MP)** (may have to generate the executable file)
- **[Circular Statistics Toolbox](https://github.com/mrkrause/circstat-matlab)**
- **[Chronux Toolbox](http://chronux.org/)**
