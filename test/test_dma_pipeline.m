%% DMA master-curve and Prony-series example
% This is a worked numerical example, not an assertion-based unit test.

close all;
clear;
clc;

%% Inputs to change
testDirectory                       = fileparts(mfilename('fullpath'));
projectDirectory                    = fileparts(testDirectory);
addpath(fullfile(projectDirectory, 'src'));

dataFile                            = fullfile(testDirectory, 'example_data','ta_dma_frequency_sweeps.xlsx');
sheetName                           = 'DMA Data';
referenceTemperature_C              = 25;
maximumSweepFrequency_Hz            = 25;           % example-data QC cutoff; [] keeps all points
medianFilterWindow                  = 3;
numberOfPronyTerms                  = 33;
relaxationTimeMode                  = 'exact';      % exact, round, manual, or min
relaxationTimeRange_s               = [];           % for example, [1e-8 1e5] in manual mode
fittingFrequencyRange_Hz            = [1e-13 1e12]; % example range before the high-frequency upturn
densePointsPerDecade                = 25;
termCountsToCompare                 = [21 27 33];


predictionTemperatures_C            = [-20 0 25 50 70];
predictionFrequencies_Hz            = [0.1 1 10];
allowPredictionExtrapolation        = false;


% Start from complete defaults and change only the preferred presentation.
plotOptions                         = dmaPlot('defaults');

% Common plot options
plotOptions.AxesLineWidth           = 1.0;
plotOptions.FontName                = 'Latin Modern';
plotOptions.AxisLabelFontSize       = 30;
plotOptions.TickLabelFontSize       = 30;
plotOptions.ColorbarFontSize        = 30;
plotOptions.LegendFontSize          = 24;
plotOptions.Interpreter             = 'latex';      % none, tex, or latex
plotOptions.Box                     = 'on';
plotOptions.Grid                    = 'on';
plotOptions.MinorTicks              = 'on';
plotOptions.LegendBox               = 'off';
plotOptions.ColorMap                = 'turbo';
plotOptions.ColorbarLocation        = 'eastoutside';
plotOptions.ModulusType             = 'E';          % E for tension/compression, G for shear
plotOptions.ModulusUnit             = 'MPa';        % label only; data are not rescaled

% Raw DMA plots: storage modulus, loss modulus, and loss factor
plotOptions.RawScale                = {'loglog', 'loglog', 'semilogx'};
plotOptions.RawLineWidth            = 3;
plotOptions.RawMarker               = 'none';
plotOptions.RawMarkerSize           = 7;
plotOptions.RawTempLimits           = [];
plotOptions.RawTempTicks            = [];
plotOptions.RawXLim                 = {[], [], []};
plotOptions.RawXTicks               = {[], [], []};
plotOptions.RawYLim                 = {[], [], []};
plotOptions.RawYTicks               = {[], [], []};

% Shifted master curve
plotOptions.ShiftedScale            = 'loglog';
plotOptions.ShiftedSecondaryYScale  = 'linear';
plotOptions.ShiftedLineWidth        = 5;
plotOptions.ShiftedMarker           = 'none';
plotOptions.ShiftedMarkerSize       = 7;
plotOptions.ShiftedLegendLocation   = 'east';
plotOptions.ShiftedTempLimits       = [];
plotOptions.ShiftedTempTicks        = [];
plotOptions.ShiftedXLim             = [];
plotOptions.ShiftedXTicks           = [];
plotOptions.ShiftedYLim             = [];
plotOptions.ShiftedYTicks           = [];
plotOptions.ShiftedSecondaryYLim    = [];
plotOptions.ShiftedSecondaryYTicks  = [];

% Filtered master curve and Maxwell fit
plotOptions.FitScale                = 'loglog';
plotOptions.FitSecondaryYScale      = 'linear';
plotOptions.FilteredMarker          = '^';
plotOptions.FilteredMarkerSize      = 7;
plotOptions.FilteredMarkerLineWidth = 2;
plotOptions.MaxwellLineWidth        = 4;
plotOptions.FitLegendLocation       = 'east';
plotOptions.FitXLim                 = [];
plotOptions.FitXTicks               = [];
plotOptions.FitYLim                 = [];
plotOptions.FitYTicks               = [];
plotOptions.FitSecondaryYLim        = [];
plotOptions.FitSecondaryYTicks      = [];

% WLF fit
plotOptions.WLFMarker               = 'o';
plotOptions.WLFMarkerSize           = 8;
plotOptions.WLFMarkerLineWidth      = 3;
plotOptions.WLFLineWidth            = 3;
plotOptions.WLFLegendLocation       = 'northeast';
plotOptions.WLFXLim                 = [-50 100];
plotOptions.WLFXTicks               = -50:25:100;
plotOptions.WLFYLim                 = [-20 20];
plotOptions.WLFYTicks               = -20:10:20;

%% 1. Import and plot the measured frequency sweeps
data                                = dmaReadTAExcel(dataFile, 'Sheet', sheetName);
dmaPlot('raw', data, plotOptions);


% The two highest-frequency points in this example show a repeated upward
% loss-modulus tail. Keep them visible above, but omit them from shifting and
% constitutive fitting when the example-data cutoff is enabled.
analysisData = data;
if ~isempty(maximumSweepFrequency_Hz)
    analysisData    = analysisData(analysisData.Frequency_Hz <= maximumSweepFrequency_Hz, :);
end



%% 2. Calculate horizontal shift factors
[shiftFactors, shiftDiagnostics]    = dmaEstimateShiftFactors(analysisData, referenceTemperature_C);
disp(shiftFactors)

shiftDiagnosticsTable               = struct2table(shiftDiagnostics);

disp(shiftDiagnosticsTable(:, {'ReferenceTemperature_C', ...
    'ShiftedTemperature_C', 'IncrementalLog10aT', 'Extrapolated', ...
    'DroppedFirstReference', 'DroppedFirstShifted'}))


%% 3. Construct and plot the shifted measured data
masterCurve                         = dmaBuildMasterCurve(analysisData, shiftFactors);
dmaPlot('shifted', masterCurve, plotOptions);



%% 4. Filter the master curve without overwriting measured values
filteredMasterCurve         = dmaFilterMasterCurve(masterCurve, 'MedianWindow', medianFilterWindow);


%% 5. Fit the WLF equation
fitReferenceTemperature_C   = shiftFactors.ReferenceTemperature_C(1);
wlf                         = dmaFitWLF(shiftFactors, fitReferenceTemperature_C);

fprintf('\nWLF fit at Tref = %.3g C\n', wlf.ReferenceTemperature_C);
fprintf('C1 = %.6g\n', wlf.C1);
fprintf('C2 = %.6g C\n', wlf.C2_C);
fprintf('R^2 = %.6f\n', wlf.RSquared);



%% 6. Fit the Prony series
prony           = dmaFitProny(filteredMasterCurve, numberOfPronyTerms, ...
                                'RelaxationTimeMode', relaxationTimeMode, ...
                                'RelaxationTimeRange', relaxationTimeRange_s, ...
                                'FittingFrequencyRange', fittingFrequencyRange_Hz);

fprintf('\nProny series\n');
fprintf('Instantaneous modulus E_0 = %.6g\n', prony.E_0);
fprintf('Equilibrium modulus E_inf = %.6g\n', prony.E_inf);
fprintf('Fit method = %s\n', prony.Method);
fprintf('RMS relative error: storage = %.4f, loss = %.4f, tan delta = %.4f\n', ...
    prony.RMSRelativeError_StorageModulus, ...
    prony.RMSRelativeError_LossModulus, ...
    prony.RMSRelativeError_TanDelta);
fprintf('Fitted tau range = [%.3e, %.3e] s\n', ...
    prony.RelaxationTimeRange_s);
fprintf('g_i is dimensionless; E_i = E_0*g_i in the input modulus units.\n');
disp(prony.Terms)



%% 7. Evaluate and plot the filtered data and Maxwell reconstruction
denseResponse                   = dmaEvaluateProny(prony, ...
                                 'FrequencyRange', prony.FrequencyRange_Hz, ...
                                 'PointsPerDecade', densePointsPerDecade);

% Plot only the measurements admitted to the fit. The complete shifted master
% curve, including the excluded high-frequency region, is shown in Section 3.
fitMasterCurve                  = filteredMasterCurve( ...
    filteredMasterCurve.ShiftedFrequency_Hz >= prony.FrequencyRange_Hz(1) & ...
    filteredMasterCurve.ShiftedFrequency_Hz <= prony.FrequencyRange_Hz(2), :);

dmaPlot('fit', fitMasterCurve, denseResponse, plotOptions);



%% 8. Compare user-selected Prony term counts
[pronyComparison, ~]    = dmaComparePronyTerms(filteredMasterCurve, ...
    termCountsToCompare, 'RelaxationTimeMode', relaxationTimeMode, ...
    'RelaxationTimeRange', relaxationTimeRange_s, ...
    'FittingFrequencyRange', fittingFrequencyRange_Hz);
disp(pronyComparison)



%% 9. Predict temperature-dependent response at selected physical frequencies
temperaturePrediction = dmaPredictTemperatureResponse(prony, wlf, ...
    predictionTemperatures_C, predictionFrequencies_Hz, ...
    'AllowExtrapolation', allowPredictionExtrapolation);
disp(temperaturePrediction)



%% 10. Plot the WLF fit
dmaPlot('wlf', wlf, plotOptions);
