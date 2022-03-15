# serotonin_LFP
Custom analysis code for [authorship]. Serotonergic modulation of local network processing in V1 shares signatures with the modulation by spatial attention, [year]

# Contents
## Data
Source data needed to recreate paper figures

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

*ex.Trials* contains the following information:
- **st**: Might remove this
- **me**: Might remove this
- **or**: Orientation used for grating stimulus for the trial
- **co**: Contrast used for grating stimulus for the trial
- **TrialEnd**: End time of trial (in sec)
- **TrialStart**: Start time of trial (in sec), i.e. animal starts fixating
- **Start**: Start time of every stimulus video frame [there is typically a small, variable delay between TrialStart (fixation onset) and the first stimulus frame presented]
- **Reward**: Specifies whether the subject response was correct (1) or the trials was discarded (0) [Only trials that are 1 are valid trials]
- **RewardSize**: Size of the reward delivered
- **Spikes**: Spike times (relative to TrialStart), in sec
-**LFP**: Measured LFP voltage for the trial
-**LFP_ts**: Time points corresponding to LFP voltages
-**times_fpOn**: Time (in sec) when fixation point begins presentation


Following preprocessing and filtering, these raw session .mat files are used to
create compiled Lfps_rc.mat files. There are corresponding versions for MP-processed
data, spike-count separated data, and combined MP/spike-count data, as well.

Each Lfps_rc.mat includes a struct with the following information:
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

## TODO List - Remove when finished
*TODO: Check if trial information is also stored in the LFP files*
