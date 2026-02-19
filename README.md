# TMS-EEG Time Series Analysis Workflow

This MATLAB project provides a complete workflow for analyzing TMS-EEG time series data, including linear and nonlinear measures, statistical comparisons, and visualization.

## Project Structure

The project is organized into modular MATLAB functions:

### Main Script
- `main_workflow.m` - Main orchestration script that runs the complete analysis pipeline

### Data Loading and Preprocessing
- `load_eeg_data.m` - Loads EEG .mat files and channel names
- `select_electrodes.m` - Selects electrodes based on AEM numbers
- `extract_channels.m` - Extracts selected channels from datasets
- `visualize_and_reject_episodes.m` - Interactive visualization for manual artifact rejection
- `select_clean_episodes.m` - Selects N clean episodes per intensity level
- `segment_data.m` - Segments data into preTMS and postTMS periods
- `select_pilot_series.m` - Randomly selects pilot series for measure selection

### Linear Analysis
- `pilot_linear_analysis.m` - Pilot analysis: ACF, PACF, ARIMA modeling, forecasting
- `compute_linear_measures_all.m` - Computes linear measures for all segments

### Nonlinear Analysis
- `pilot_nonlinear_analysis.m` - Pilot analysis: mutual information, FNN, correlation dimension, nonlinear forecasting
- `compute_nonlinear_measures_all.m` - Computes nonlinear measures for all segments
- `compute_mutual_information.m` - Computes mutual information between X(t) and X(t+lag)
- `compute_fnn.m` - Computes False Nearest Neighbors ratio
- `estimate_correlation_dimension.m` - Estimates correlation dimension using Grassberger-Procaccia algorithm
- `compute_nonlinear_forecasting_nrmse.m` - Computes NRMSE for local linear nonlinear forecasting

### Statistical Analysis
- `perform_statistical_analysis.m` - Performs statistical comparisons and creates visualizations

### Data Files
- `Channels.txt` - Channel names file (first 4 channels are non-EEG)

## Setup Instructions

1. **Prepare your data files:**
   - Place your three EEG .mat files in a `data` folder within the Timeseries directory
   - Each .mat file should contain a `cMF` array
   - Update `INTENSITY_LEVELS` in `main_workflow.m` to match your file names

2. **Update configuration in `main_workflow.m`:**
   - `AEM_NUMBERS` - Your team members' AEM numbers (currently [1, 1])
   - `DATA_FOLDER` - Path to folder containing .mat files (default: 'data')
   - `CHANNELS_FILE` - Path to Channels.txt (default: 'Channels.txt')
   - `N_CLEAN_EPISODES` - Number of clean episodes to select (default: 20)

3. **Run the analysis:**
   ```matlab
   main_workflow
   ```

## Workflow Steps

1. **Data Loading**: Loads three EEG .mat files and channel names
2. **Electrode Selection**: Selects two electrodes based on AEM numbers using mod(x,60)+1 rule
3. **Channel Extraction**: Extracts only the selected channels from all episodes
4. **Visualization & Rejection**: Interactive interface to review and reject episodes with artifacts
5. **Clean Episode Selection**: Randomly selects N clean episodes per intensity level
6. **Data Segmentation**: Segments data into preTMS (first 900ms) and postTMS (last 900ms)
7. **Pilot Series Selection**: Randomly selects one preTMS and one postTMS segment from each intensity level
8. **Pilot Linear Analysis**: Computes ACF, PACF, fits ARIMA models, generates forecasts, computes NRMSE
9. **Pilot Nonlinear Analysis**: Computes mutual information, FNN, correlation dimension, nonlinear forecasting
10. **Measure Selection**: Selects one linear and one nonlinear measure based on pilot analysis
11. **Full Analysis**: Computes selected measures for all segments
12. **Statistical Comparisons**: Compares preTMS vs postTMS, across intensity levels, and differences
13. **Visualization**: Creates histograms, boxplots, and scatter plots for all comparisons

## Interactive Episode Review

When reviewing episodes, use the following keyboard controls:
- **[Space]** or **[Right Arrow]** - Next episode
- **[Left Arrow]** - Previous episode
- **[R]** - Toggle reject/accept current episode
- **[Q]** - Quit and save rejections
- **[A]** - Accept all remaining episodes

## Output

The analysis produces:
- ACF/PACF plots for pilot series
- Comparison plots of ACF, PACF, and NRMSE curves
- Nonlinear measure comparison plots
- Statistical comparison plots (histograms, boxplots, scatter plots) for:
  - preTMS vs postTMS within each intensity level
  - Comparison across intensity levels for postTMS
  - Comparison across intensity levels for preTMS
  - Comparison using difference (postTMS - preTMS)

## Notes

- The code automatically handles different data formats ([channels x time x episodes] or [episodes x channels x time])
- Sampling frequency is set to 1450 Hz (can be modified in `main_workflow.m`)
- PreTMS and postTMS durations are 900ms each (can be modified)
- The selected linear and nonlinear measures can be customized in the pilot analysis functions

## Requirements

- MATLAB R2018b or later
- Statistics and Machine Learning Toolbox (for ARIMA, boxplot, etc.)
- Signal Processing Toolbox (for autocorr, parcorr)

## Customization

To change the selected measures, modify:
- `pilot_linear_analysis.m` - Change `selected_measure` and `measure_params`
- `pilot_nonlinear_analysis.m` - Change `selected_measure` and `measure_params`

Available linear measures:
- `NRMSE_step5` - NRMSE at forecast step 5
- `ACF_lag10` - ACF at lag 10
- `PACF_lag5` - PACF at lag 5

Available nonlinear measures:
- `MI_first_min_lag` - Lag of first minimum in mutual information
- `FNN_embedding_dim` - Embedding dimension where FNN drops below 1%
- `Correlation_dimension` - Correlation dimension estimate
- `Nonlinear_NRMSE` - NRMSE of nonlinear forecasting model
