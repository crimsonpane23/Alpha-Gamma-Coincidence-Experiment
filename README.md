# Alpha–Gamma Coincidence Experiment

This repository contains analysis code and supporting material for an undergraduate nuclear-physics experiment on the range of alpha particles in air and alpha–gamma coincidence measurements in the decay of \(^{241}\mathrm{Am}\) to \(^{237}\mathrm{Np}\).

**Repository:** https://github.com/crimsonpane23/Alpha-Gamma-Coincidence-Experiment

## Experiment overview

The work has two parts:

1. **Range of alpha particles:** The residual energy of alpha particles from an \(^{241}\mathrm{Am}\) source was measured at different air pressures while keeping the source–detector separation fixed. The effective path length was calculated for standard atmospheric conditions, and the range was estimated by extrapolating the residual-energy curve.
2. **Alpha–gamma coincidence:** Alpha particles were detected using a silicon surface-barrier detector and gamma rays using a 2-inch × 2-inch NaI(Tl) scintillation detector. Coincidence events were selected using pulse-height calibration and alpha–gamma energy correlations. Several timing definitions were compared, with pulse-height-dependent time-walk correction applied before lifetime fitting.

## Main results

- **Extrapolated alpha-particle range in air at STP:** \((2.61 \pm 0.01)\,\mathrm{cm}\).
- **Coincidence data:** 4,963 events recorded over 96 hours.
- **Reference transition:** the 59.5 keV excited state of \(^{237}\mathrm{Np}\), with a reference half-life of approximately 67.2 ns.
- **Extracted half-lives:** \((63 \pm 2)\,\mathrm{ns}\) using trigger-crossing times and \((82 \pm 3)\,\mathrm{ns}\) using constant-fraction times. The trigger-time estimate is closer to the reference value. The difference between methods indicates that timing-related systematic effects remain relevant.

These values are the results reported in the accompanying experimental report; uncertainties should be interpreted according to the analysis and fit procedures described there.

## Analysis code

The C++ analysis macros use [ROOT](https://root.cern/) for data processing, histograms, fits, and plots.

### 1. `analyse_coincidence_v3.C`

Processes oscilloscope CSV waveforms for the alpha and gamma channels. It estimates baseline RMS, pulse height, peak time, fixed-threshold crossing time, and constant-fraction time; applies event-quality cuts; and writes event-level CSV files and diagnostic plots.

The macro currently defines the input directory and acquisition settings near the top of the source file. **Edit `FOLDER` to point to your local raw-data directory before running it.** Check channel prefixes, trigger thresholds, pulse polarities, RMS cuts, and saturation settings against your own data.

Example command (from the directory containing the macro):

```bash
root -l -b -q 'analyse_coincidence_v3.C+'
```

### 2. `correlation.C`

Builds alpha–gamma energy-correlation heatmaps and selects events within a specified region of the correlation plot. The selection limits and energy-calibration settings are defined in the source. Inspect the heatmaps first and set the selection bounds to the region justified by the decay scheme and measured data.

```bash
root -l -b -q 'correlation.C+'
```

The macro expects the coincidence CSV files and column names described in its source comments. Run it from a working directory containing the expected input files, or update the input paths accordingly.

### 3. `lifetime_v3.C`

Reads the selected correlated events, calculates time differences for peak, fixed-trigger, and constant-fraction timing, applies amplitude-dependent time-walk correction, and fits the selected portions of the time-difference distributions with an exponential function. It produces timing-correction plots, fit plots, and fit-result text output.

```bash
root -l -b -q 'lifetime_v3.C+'
```

The input filename, timing-fit intervals, plot ranges, and correction settings are configurable near the top of the macro. Verify these settings before interpreting the fitted lifetime.

> **Note:** The filenames used in the example commands must match the filenames in your local checkout. Run each macro from the directory where its configured input files are available. ROOT must be installed and accessible as `root` in your shell.

## Typical analysis sequence

1. Place the raw oscilloscope CSV files in the expected directory and update the path and channel settings in `analyse_coincidence_v3.C`.
2. Run the waveform analysis and inspect the generated CSV files and event-check plots.
3. Run `correlation.C` to inspect the linear/log-scale heatmaps, set the energy-selection bounds, and produce the selected-event CSV.
4. Run `lifetime_v3.C` to compare the timing methods, apply time-walk correction, and inspect the exponential-fit regions and fit results.
5. Preserve the configuration values and selected fit intervals with the output plots so that the analysis can be reproduced.

The exact file names and directory layout of raw data and generated outputs depend on the local configuration in the macros.

## Experimental setup

- **Alpha source:** \(^{241}\mathrm{Am}\).
- **Alpha detector:** silicon surface-barrier detector.
- **Gamma detector:** 2-inch × 2-inch NaI(Tl) scintillation detector (SD152F).
- **Gamma-detector high voltage:** 600 V.
- **Waveform acquisition:** Teledyne LeCroy WaveSurfer oscilloscope; the final configuration used 40 kS, 2 GS/s, and a 20 µs acquisition window.
- **Trigger thresholds:** 0.6 V for alpha and 2.3 V for gamma in the reported final configuration.

These are the settings documented in the report and may need to be changed for other runs or setups.

## Dependencies

- ROOT with the C++ interpreter and histogram/fit libraries used by the macros.
- C++17-compatible compiler support for `std::filesystem` if compiling outside ROOT's interpreter/ACLiC workflow.
- Raw oscilloscope CSV files and intermediate CSV files required by each analysis stage.

## Reports and references

The repository is associated with the report **“Range of α-particles and setting up α–γ coincidence”** (October 2026). The report describes the apparatus, data-acquisition procedure, energy calibration, event selection, time-walk correction, and lifetime fits.

A related reference for FPGA-based coincidence counting and nuclear-lifetime measurements is:

M. Vretenar, N. Erceg, and M. Karuza, “Energy-resolved coincidence counting using an FPGA for nuclear lifetime experiments,” *American Journal of Physics* **87**, 997–1004 (2019). https://doi.org/10.1119/1.5122744

## Limitations

The FPGA-based acquisition route was explored but not used for the final reported coincidence data; the oscilloscope was used instead. The extracted lifetime depends on the timing method, and residual timing-related systematic effects require further investigation. The README does not replace the detailed analysis choices and caveats documented in the report.
