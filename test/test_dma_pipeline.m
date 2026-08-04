%% DMA master-curve and Prony-series example
% This is a worked numerical example, not an assertion-based unit test.

close all;
clear;
clc;

%% Inputs to change
testDirectory                    = fileparts(mfilename('fullpath'));
dataFile                         = fullfile(testDirectory, 'example_data', ...
    'ta_dma_frequency_sweeps.xlsx');
sheetName                        = 'DMA Data';
referenceTemperature_C           = 25;

medianFilterWindow               = 5;
numberOfPronyTerms               = 15;
relaxationTimeMode               = 'exact';  % exact, round, manual, or min
relaxationTimeRange_s             = [];       % for example, [1e-8 1e5] in manual mode
fittingFrequencyRange_Hz          = [];       % [] uses every shifted frequency
densePointsPerDecade              = 25;
termCountsToCompare               = [5 10 15];

predictionTemperatures_C          = [-20 0 25 50 70];
predictionFrequencies_Hz          = [0.1 1 10];
allowPredictionExtrapolation      = false;

% Publication-style plotting options used by every dmaPlot call.
plotOptions                       = struct;
plotOptions.Scale                 = 'loglog';
plotOptions.LineWidth             = 2;
plotOptions.AxesLineWidth         = 1.0;
plotOptions.FontName              = 'Latin Modern';
plotOptions.AxisLabelFontSize     = 20;
plotOptions.TickLabelFontSize     = 20;
plotOptions.LegendFontSize        = 16;
plotOptions.LegendInterpreter     = 'latex';
plotOptions.XLabelInterpreter     = 'latex';
plotOptions.YLabelInterpreter     = 'latex';
plotOptions.TickLabelInterpreter  = 'latex';
plotOptions.TitleInterpreter      = 'latex';
plotOptions.Box                   = 'on';
plotOptions.Grid                  = 'on';
plotOptions.MinorTicks            = 'on';
plotOptions.ShowLegend            = 'on';
plotOptions.LegendLocation        = 'best';
plotOptions.LegendBox             = 'off';
plotOptions.LegendNumColumns      = 3;
plotOptions.XLim                  = [];
plotOptions.YLim                  = [];
plotOptions.XTicks                = [];
plotOptions.YTicks                = [];

plotVariables                     = {'StorageModulus', 'LossModulus', 'TanDelta'};
filteredVariables                 = {'StorageModulusFiltered', ...
    'LossModulusFiltered', 'TanDeltaFiltered'};
plotNames                         = {'Storage modulus', 'Loss modulus', 'Tan delta'};
plotScales                        = {'loglog', 'loglog', 'semilogx'};

projectDirectory                  = fileparts(testDirectory);
addpath(fullfile(projectDirectory, 'src'));

%% 1. Import and plot the measured frequency sweeps
data = dmaReadTAExcel(dataFile, 'Sheet', sheetName);

for k = 1:numel(plotVariables)
    figureHandle = figure('Name', ['Raw DMA - ', plotNames{k}], ...
        'NumberTitle', 'off', 'Color', 'w');
    plotOptions.Scale = plotScales{k};
    dmaPlot(data, 'XVariable', 'Frequency_Hz', ...
        'YVariables', plotVariables(k), 'GroupVariable', 'Set', ...
        'Axes', axes('Parent', figureHandle), 'PlotOptions', plotOptions);
end

%% 2. Calculate horizontal shift factors
[shiftFactors, shiftDiagnostics] = dmaEstimateShiftFactors( ...
    data, referenceTemperature_C);
disp(shiftFactors)

shiftDiagnosticsTable = struct2table(shiftDiagnostics);
disp(shiftDiagnosticsTable(:, {'ReferenceTemperature_C', ...
    'ShiftedTemperature_C', 'IncrementalLog10aT', 'Extrapolated', ...
    'DroppedFirstReference', 'DroppedFirstShifted'}))

%% 3. Construct and plot the shifted measured data
masterCurve = dmaBuildMasterCurve(data, shiftFactors);

for k = 1:numel(plotVariables)
    figureHandle = figure('Name', ['Shifted DMA - ', plotNames{k}], ...
        'NumberTitle', 'off', 'Color', 'w');
    plotOptions.Scale = plotScales{k};
    dmaPlot(masterCurve, 'XVariable', 'ShiftedFrequency_Hz', ...
        'YVariables', plotVariables(k), 'GroupVariable', 'Set', ...
        'Axes', axes('Parent', figureHandle), 'PlotOptions', plotOptions);
end

%% 4. Filter the master curve without overwriting measured values
filteredMasterCurve = dmaFilterMasterCurve(masterCurve, ...
    'MedianWindow', medianFilterWindow);

%% 5. Fit the WLF equation
fitReferenceTemperature_C = shiftFactors.ReferenceTemperature_C(1);
wlf = dmaFitWLF(shiftFactors, fitReferenceTemperature_C);

fprintf('\nWLF fit at Tref = %.3g C\n', wlf.ReferenceTemperature_C);
fprintf('C1 = %.6g\n', wlf.C1);
fprintf('C2 = %.6g C\n', wlf.C2_C);
fprintf('R^2 = %.6f\n', wlf.RSquared);

%% 6. Fit the Prony series
prony = dmaFitProny(filteredMasterCurve, numberOfPronyTerms, ...
    'RelaxationTimeMode', relaxationTimeMode, ...
    'RelaxationTimeRange', relaxationTimeRange_s, ...
    'FittingFrequencyRange', fittingFrequencyRange_Hz);

fprintf('\nProny series\n');
fprintf('Instantaneous modulus E_0 = %.6g\n', prony.E_0);
fprintf('Equilibrium modulus E_inf = %.6g\n', prony.E_inf);
fprintf('Fitted tau range = [%.3e, %.3e] s\n', ...
    prony.RelaxationTimeRange_s);
fprintf('g_i is dimensionless; E_i = E_0*g_i in the input modulus units.\n');
disp(prony.Terms)

%% 7. Evaluate and plot the smooth generalized-Maxwell reconstruction
denseResponse = dmaEvaluateProny(prony, ...
    'FrequencyRange', prony.FrequencyRange_Hz, ...
    'PointsPerDecade', densePointsPerDecade);

for k = 1:numel(plotVariables)
    figureHandle = figure('Name', ['Prony fit - ', plotNames{k}], ...
        'NumberTitle', 'off', 'Color', 'w');
    ax = axes('Parent', figureHandle);

    fitPlotOptions = plotOptions;
    fitPlotOptions.Scale = plotScales{k};
    fitPlotOptions.ShowLegend = 'off';
    dmaPlot(denseResponse, 'XVariable', 'Frequency_Hz', ...
        'YVariables', plotVariables(k), 'GroupVariable', '', ...
        'Axes', ax, 'PlotOptions', fitPlotOptions);

    fitLine = findobj(ax, 'Type', 'Line');
    set(fitLine(1), 'DisplayName', 'Generalized-Maxwell reconstruction');
    measuredLine = plot(ax, filteredMasterCurve.ShiftedFrequency_Hz, ...
        filteredMasterCurve.(filteredVariables{k}), 'o', ...
        'LineStyle', 'none', 'MarkerSize', 4, 'Color', [0.25 0.25 0.25], ...
        'DisplayName', 'Filtered master data');
    legend(ax, [measuredLine, fitLine(1)], 'Location', 'best', ...
        'Interpreter', plotOptions.LegendInterpreter, ...
        'FontSize', plotOptions.LegendFontSize, ...
        'Box', plotOptions.LegendBox, 'NumColumns', 1);
end

%% 8. Compare user-selected Prony term counts
[pronyComparison, ~] = dmaComparePronyTerms(filteredMasterCurve, ...
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
figure;
clf;
hold on;
box on;
grid off;

[temperaturePlot, order] = sort(wlf.Data.Temperature_C);
plot(wlf.Data.Temperature_C, wlf.Data.MeasuredLog10aT, 'o', ...
    'MarkerSize', 8, 'LineWidth', 3, 'DisplayName', 'Calculated shift factors');

plot(temperaturePlot, wlf.Data.FittedLog10aT(order), '-', ...
    'LineWidth', 3, 'DisplayName', 'WLF fit');

xlabel('Temperature ($^{\circ}$C)', 'interpreter','latex', 'FontSize',20);
ylabel('$\log_{10}(a_T)$','interpreter','latex', 'FontSize',20);
set(gca, 'xlim', [-50, 100], 'xtick', -50:25:100, 'XMinorTick', 'on', ...
    'ylim', [-20, 20], 'ytick', -20:10:20, 'YMinorTick','on', ...
    'TickLabelInterpreter', 'latex', 'FontSize', 20)
legend('fontsize',16, 'Location', 'northeast', 'Box','off');
hold off;
