# MATLAB API reference

Add the source directory before calling the functions:

```matlab
addpath(fullfile(projectDirectory, 'src'));
```

The examples below use canonical table-variable names returned by
`dmaReadTAExcel`. See [methods.md](methods.md) for equations and assumptions.
The method lineage and PyVisco citation are documented in the
[reference section of methods.md](methods.md#reference-implementation-and-citation).

## `dmaReadTAExcel`

Read one numerical TA/TRIOS-style Excel worksheet.

```matlab
data = dmaReadTAExcel(filename)
data = dmaReadTAExcel(filename, 'Sheet', sheet)
```

### Input

- `filename`: Excel workbook path.
- `Sheet`: worksheet name or number; default `1`.

Recognized quantities include temperature, cyclic/angular frequency, storage
and loss modulus, tan delta, oscillation strain/stress, amplitude, time, and
set/step identifiers. The output uses canonical variables such as
`Temperature_C`, `Frequency_Hz`, `StorageModulus`, and `LossModulus`.

The function derives cyclic/angular frequency from the other when only one is
present and derives tan delta when both moduli are present.

## `dmaPlot`

Plot one or more table variables, optionally grouped by set or another column.

```matlab
plotOptions = dmaPlot('defaults');
plotOptions.ModulusType = 'E';
plotOptions.ModulusUnit = 'MPa';

ax = dmaPlot(data, ...
    'XVariable', 'Frequency_Hz', ...
    'YVariables', {'StorageModulus', 'LossModulus'}, ...
    'SecondaryYVariables', {'TanDelta'}, ...
    'GroupVariable', 'Set', ...
    'PlotOptions', plotOptions)
```

Calling `dmaPlot('defaults')` returns the complete default structure. A partial
structure can also be passed directly; every omitted field is filled from the
same internal defaults.

### Name-value inputs

- `XVariable`: table column used for the horizontal axis. If omitted, a
  suitable available variable is selected.
- `YVariables`: one name or a cell/string array of response columns.
- `SecondaryYVariables`: optional response columns plotted on the right y-axis.
- `GroupVariable`: grouping column, `auto`, or an empty string for no grouping.
- `Scale`: `auto`, `linear`, `semilogx`, `semilogy`, or `loglog`. A direct
  `Scale` input overrides `PlotOptions.Scale`.
- `Axes`: existing axes handle; default current axes.
- `PlotOptions`: publication-style option structure.

### `PlotOptions` fields

| Field | Purpose |
|---|---|
| `Scale` | Axis scaling mode |
| `SecondaryYScale` | `linear` or `log` scale for the right y-axis |
| `LineWidth` | General data-line width |
| `Marker`, `MarkerSize` | General marker symbol and size |
| `ShiftedLineWidth`, `ShiftedMarker` | Shifted-master-curve line width and marker used by the example |
| `FilteredMarker`, `FilteredMarkerSize`, `FilteredMarkerLineWidth` | Filtered-master-curve marker settings used by the example |
| `MaxwellLineWidth` | Maxwell-fit line width used by the example |
| `AxesLineWidth` | Axes border/tick width |
| `FontName` | Axes, label, and legend font |
| `AxisLabelFontSize` | Axis-label size |
| `TickLabelFontSize` | Tick-label size |
| `LegendFontSize` | Legend size |
| `TitleFontSize` | Title size if a title is supplied externally |
| `ColorbarFontSize` | Temperature-colorbar tick size |
| `Interpreter` | shared label, tick, legend, title, and colorbar interpreter: `none`, `tex`, or `latex` |
| `Box` | `on` or `off` |
| `Grid` | `on` or `off` |
| `MinorTicks` | `on` or `off` |
| `ShowLegend` | `on` or `off` |
| `LegendLocation` | MATLAB legend location; default `best` |
| `ShiftedLegendLocation` | response-legend location for the shifted plot in the example |
| `FitLegendLocation` | legend location for the filtered/fitted plot in the example |
| `LegendBox` | `on` or `off` |
| `TemperatureColorbar` | `auto`, `on`, or `off`; `auto` uses it for grouped isotherms |
| `ColorMap`, `ColorbarLocation` | temperature colormap and colorbar placement |
| `TemperatureLimits`, `TemperatureTicks` | optional colorbar range and ticks |
| `ModulusType` | `E` for tension/compression or `G` for shear labels |
| `ModulusUnit` | modulus-unit label such as `Pa`, `kPa`, or `MPa` |
| `ResponseLineStyles`, `ResponseColors` | visual distinction between responses |
| `XLim`, `YLim`, `SecondaryYLim` | empty or increasing two-value limits |
| `XTicks`, `YTicks`, `SecondaryYTicks` | empty or increasing tick vectors |

`ModulusType` and `ModulusUnit` change labels only; numerical values are never
rescaled. When the temperature colorbar is active, it replaces the long list
of isotherm legend entries. A compact response legend is retained when several
quantities share one figure.

## `dmaGroupIsotherms`

Assign and summarize isothermal sets.

```matlab
[data, isotherms] = dmaGroupIsotherms(data)
[data, isotherms] = dmaGroupIsotherms(data, ...
    'TemperatureTolerance', 0.5)
```

An existing finite `Set` variable is retained. Otherwise, temperatures are
grouped using `TemperatureTolerance`. `isotherms` contains `Set`, mean
`Temperature_C`, and `NumberOfPoints`.

## `dmaEstimateShiftFactors`

Estimate horizontal shift factors from neighboring storage-modulus isotherms.

```matlab
[shiftFactors, diagnostics] = ...
    dmaEstimateShiftFactors(data, referenceTemperature_C)
```

### Name-value inputs

- `TemperatureTolerance`: grouping tolerance; default `0.5` C.
- `ReferenceTolerance`: maximum distance from a measured isotherm; default
  `1.0` C.
- `NumberOfLevels`: common modulus levels per pair; default `10`.
- `DropFirstPoint`: compare fits with and without the first point; default
  `true`.

### Outputs

`shiftFactors` contains `Set`, `Temperature_C`, `log10_aT`, `Method`, and the
actual measured reference temperature. `diagnostics` is a structure array with
pairwise fitted parameters, standard errors, drop decisions, overlap flags,
modulus levels, and inferred frequencies.

Returned shift factors can be edited directly before calling
`dmaBuildMasterCurve`.

## `dmaBuildMasterCurve`

Apply a shift-factor table without changing the measured response columns.

```matlab
masterCurve = dmaBuildMasterCurve(data, shiftFactors)
masterCurve = dmaBuildMasterCurve(data, shiftFactors, ...
    'TemperatureTolerance', 0.5)
```

Adds `log10_aT`, `ShiftedFrequency_Hz`, and
`ShiftedAngularFrequency_rad_s`, then sorts by shifted frequency.

## `dmaFilterMasterCurve`

Apply explicit moving-median filtering while retaining measured columns.

```matlab
filteredMaster = dmaFilterMasterCurve(masterCurve)
filteredMaster = dmaFilterMasterCurve(masterCurve, ...
    'MedianWindow', 5)
```

Default `MedianWindow` is `5`. The function adds
`StorageModulusFiltered`, `LossModulusFiltered`, and `TanDeltaFiltered`.
Use a window of `1` to disable effective filtering.

`dmaFitProny` does not filter internally. Call this function first when a
filtered fit is desired.

## `dmaWLF`

Evaluate the WLF shift relation.

```matlab
log10aT = dmaWLF(temperature_C, referenceTemperature_C, C1, C2)
```

Inputs may be scalars or compatible numeric arrays.

## `dmaFitWLF`

Fit WLF coefficients to a shift-factor table.

```matlab
wlf = dmaFitWLF(shiftFactors, referenceTemperature_C)
wlf = dmaFitWLF(shiftFactors, referenceTemperature_C, ...
    'TemperatureRange', [-35 70], ...
    'ExcludeTemperatures', [62.5])
```

### Outputs

The returned structure contains `C1`, `C2_C`, reference temperature, RMSE,
\(R^2\), SSE, solver information, covariance, standard errors, approximate 95%
confidence intervals, and a data/residual table.

## `dmaFitProny`

Fit generalized-Maxwell/Prony parameters to a master curve.

```matlab
prony = dmaFitProny(masterCurve, numberOfTerms)
prony = dmaFitProny(masterCurve, numberOfTerms, ...
    'RelaxationTimeMode', 'manual', ...
    'RelaxationTimeRange', [1e-8 1e5], ...
    'FittingFrequencyRange', [1e-5 1e8])
```

### Required master-curve variables

- `ShiftedFrequency_Hz`
- storage-modulus column
- loss-modulus column

If `StorageModulusFiltered` and `LossModulusFiltered` exist, they are selected
automatically. Otherwise, the measured columns are used.

### Name-value inputs

- `StorageVariable`: explicit storage column or `auto`.
- `LossVariable`: explicit loss column or `auto`.
- `RelaxationTimeMode`: `exact` (default), `round`, `manual`, or `min`.
- `RelaxationTimeRange`: increasing two-value range in seconds; required for
  `manual` and ignored with a warning in other modes.
- `FittingFrequencyRange`: optional increasing frequency range in hertz.

Filtering is not an option of this function. Use `dmaFilterMasterCurve`
explicitly before fitting.

### Principal outputs

- `Terms`: table containing `Term`, `tau_i_s`, `f_i_Hz`, `g_i`, and `E_i`.
- `E_0`, `E_inf`: instantaneous and equilibrium modulus.
- `RMSE_StorageModulus`, `RMSE_LossModulus`, `RMSE_TanDelta`: dimensional
  modulus errors and absolute tan-delta error.
- `NRMSE_StorageModulus`, `NRMSE_LossModulus`: range-normalized errors.
- `RMSRelativeError_StorageModulus`, `RMSRelativeError_LossModulus`,
  `RMSRelativeError_TanDelta`: point-relative errors used to judge the fit.
- `FrequencyRange_Hz`, `RelaxationTimeRange_s`: ranges actually fitted.
- `FittedCurve`: measured and fitted responses at fitting frequencies.

## `dmaEvaluateProny`

Evaluate a fitted Prony model without refitting.

```matlab
response = dmaEvaluateProny(prony)
response = dmaEvaluateProny(prony, ...
    'FrequencyRange', [1e-6 1e8], ...
    'PointsPerDecade', 25)
response = dmaEvaluateProny(prony, ...
    'Frequencies', [0.1 1 10])
```

`FrequencyRange` and `Frequencies` are mutually exclusive. Without either, the
fitted Prony frequency range is used. The output table contains frequency,
angular frequency, storage modulus, loss modulus, tan delta, corresponding
`Time_s = 1/f`, and relaxation modulus.

## `dmaComparePronyTerms`

Fit several user-selected term counts using otherwise identical options.

```matlab
[comparison, fits] = dmaComparePronyTerms(masterCurve, ...
    [5 10 15], 'RelaxationTimeMode', 'exact')
```

Additional name-value arguments are passed to `dmaFitProny`. `comparison`
reports dimensional and normalized errors; `fits` contains each complete Prony
result. The function does not automatically choose a term count.

## `dmaPredictTemperatureResponse`

Predict frequency-domain properties at requested physical temperatures and
frequencies.

```matlab
prediction = dmaPredictTemperatureResponse( ...
    prony, wlf, temperatures_C, frequencies_Hz)
prediction = dmaPredictTemperatureResponse( ...
    prony, wlf, temperatures_C, frequencies_Hz, ...
    'AllowExtrapolation', true)
```

The function forms the Cartesian product of temperatures and frequencies. The
output contains physical and reduced frequency, `log10_aT`, `aT`, storage
modulus, loss modulus, tan delta, and coverage flags.

`AllowExtrapolation` defaults to `false`. Unsupported predictions are returned
as `NaN` and marked by `PredictionEvaluated=false`.
