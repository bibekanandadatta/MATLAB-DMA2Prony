# DMA master curves and Prony fitting in MATLAB

A compact MATLAB package for processing dynamic mechanical analysis (DMA)
frequency sweeps, constructing time-temperature superposition master curves,
fitting Williams-Landel-Ferry (WLF) shift factors, and identifying generalized
Maxwell/Prony-series parameters.

The package is numerical and script-based. It has no graphical user interface
and keeps each analysis stage directly callable from a MATLAB script or live
script.

## Capabilities

- Import a numerical TA Instruments/TRIOS-style Excel results table.
- Plot storage modulus, loss modulus, and tan delta with publication controls.
- Estimate horizontal shift factors from isothermal frequency sweeps.
- Construct and optionally median-filter a master curve.
- Fit WLF coefficients and generalized-Maxwell/Prony parameters.
- Evaluate smooth frequency- and time-domain model responses.
- Compare user-selected Prony term counts.
- Predict temperature-dependent properties at selected loading frequencies.
- Fit an existing master curve without requiring the original sweep data.

## Requirements

- MATLAB; the package has been verified with MATLAB R2025b.
- Optimization Toolbox for `lsqnonneg`, which is required by the Prony fit.
- The optional `min` relaxation-time mode also requires `fmincon`.

The shift-factor and WLF routines use `lsqcurvefit` when it is available and
otherwise fall back to `fminsearch`.

## Repository layout

```text
.
├── LICENSE.md
├── THIRD_PARTY_NOTICES.md
├── README.md
├── docs/
│   ├── methods.md
│   └── api.md
├── src/
│   └── MATLAB numerical functions
└── test/
    ├── test_dma_pipeline.m
    └── example_data/
```

## Quick start

Clone or download the repository, start MATLAB in its root directory, and run:

```matlab
addpath(fullfile(pwd, 'src'));
run(fullfile(pwd, 'test', 'test_dma_pipeline.m'));
```

The worked example collects the settings intended for user adjustment at the
top of the script. It is an analysis example rather than an assertion-heavy
unit test.

A minimal calculation is:

```matlab
addpath('src');

data = dmaReadTAExcel( ...
    fullfile('test', 'example_data', 'ta_dma_frequency_sweeps.xlsx'), ...
    'Sheet', 'DMA Data');

[shiftFactors, shiftDiagnostics] = ...
    dmaEstimateShiftFactors(data, 25);
masterCurve = dmaBuildMasterCurve(data, shiftFactors);
filteredMaster = dmaFilterMasterCurve(masterCurve, 'MedianWindow', 5);

wlf = dmaFitWLF(shiftFactors, shiftFactors.ReferenceTemperature_C(1));
prony = dmaFitProny(filteredMaster, 15, ...
    'RelaxationTimeMode', 'exact');
response = dmaEvaluateProny(prony, 'PointsPerDecade', 25);
```

For an existing master curve, skip shift-factor estimation and master-curve
construction. Supply `ShiftedFrequency_Hz`, `StorageModulus`, and
`LossModulus` directly to the filtering and Prony functions.

## Documentation

- [Methods and numerical conventions](docs/methods.md)
- [MATLAB API reference](docs/api.md)



## Example data and acknowledgment

This package adapts methods from the
[PyVisco](https://github.com/NatLabRockies/pyvisco) project and adds direct
TA/TRIOS-style Excel import, arbitrary user-defined Prony relaxation-time and
fitting-frequency ranges, and direct temperature-dependent predictions at
selected frequencies with calibration-range flags. The example workbook
contains numerical measurements from
[`freq_user_raw.csv`](https://github.com/NatLabRockies/pyvisco/blob/main/sample_data/freq_user_raw.csv)
in the [PyVisco](https://github.com/NatLabRockies/pyvisco) project:

> Springer, M. (2022). *PYVISCO: A Python library for identifying Prony series
> parameters of linear viscoelastic materials* [Computer software]. Zenodo.
> [https://doi.org/10.5281/zenodo.6384954](https://doi.org/10.5281/zenodo.6384954)

The retained PyVisco copyright and license notice is provided in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## License

This project is distributed under the
[BSD 3-Clause License](LICENSE.md). Third-party notices are listed separately in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
