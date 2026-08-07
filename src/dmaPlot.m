function plotHandles = dmaPlot(plotType, varargin)
%DMAPLOT Create the standard plots used by the DMA analysis pipeline.
%
% plotOptions = dmaPlot('defaults')
% dmaPlot('raw', data, plotOptions)
% dmaPlot('shifted', masterCurve, plotOptions)
% dmaPlot('fit', filteredMasterCurve, denseResponse, plotOptions)
% dmaPlot('wlf', wlf, plotOptions)
%
% The raw mode creates separate storage-modulus, loss-modulus, and loss-factor
% figures. The other modes each create one figure. PlotOptions.ModulusType and
% PlotOptions.ModulusUnit change labels only; the data are not rescaled.

    if nargin == 1 && isTextScalar(plotType) && strcmpi(plotType, 'defaults')
        plotHandles = defaultPlotOptions;
        return
    end

    if ~isTextScalar(plotType)
        error('dmaPlot:InvalidPlotType', 'The plot type must be text.');
    end

    switch lower(char(plotType))
        case 'raw'
            options     = resolvePlotOptions(varargin{2});
            plotHandles = plotRawData(varargin{1}, options);

        case 'shifted'
            options     = resolvePlotOptions(varargin{2});
            plotHandles = plotShiftedData(varargin{1}, options);

        case {'fit', 'filtered'}
            options     = resolvePlotOptions(varargin{3});
            plotHandles = plotFilteredFit(varargin{1}, varargin{2}, options);

        case 'wlf'
            options     = resolvePlotOptions(varargin{2});
            plotHandles = plotWLF(varargin{1}, options);

        otherwise
            error('dmaPlot:InvalidPlotType', ...
                'Plot type must be raw, shifted, fit, or wlf.');
    end
end

function handles = plotRawData(data, options)
% Plot the three measured DMA quantities in separate figures.

    requireTableVariables(data, {'Frequency_Hz', 'StorageModulus', ...
        'LossModulus', 'TanDelta', 'Temperature_C'});

    [groupIndex, temperatures, colors, colorMap, temperatureLimits] = ...
        temperatureGroups(data, options.RawTempLimits, options);

    variables   = {'StorageModulus', 'LossModulus', 'TanDelta'};
    names       = {'Storage modulus', 'Loss modulus', 'Loss factor'};
    handles     = initializeHandles(3);
    labelArgs   = labelArguments(options);

    for k = 1:3
        handles.Figures(k) = figure('Name', ['Raw DMA - ', names{k}], ...
            'NumberTitle', 'off', 'Color', 'w');
        handles.Axes(k)    = axes('Parent', handles.Figures(k));
        ax                 = handles.Axes(k);

        [xScale, yScale] = splitScale(options.RawScale{k});
        prepareAxes(ax, xScale, yScale, options);

        for g = 1:numel(temperatures)
            rows            = groupIndex == g;
            [frequency, ix] = sort(data.Frequency_Hz(rows));
            response        = data.(variables{k})(rows);
            response        = response(ix);
            valid           = isfinite(frequency) & isfinite(response);
            if strcmp(xScale, 'log'), valid = valid & frequency > 0; end
            if strcmp(yScale, 'log'), valid = valid & response > 0; end

            plot(ax, frequency(valid), response(valid), ...
                'LineStyle', '-', ...
                'Marker', options.RawMarker, ...
                'MarkerSize', options.RawMarkerSize, ...
                'LineWidth', options.RawLineWidth, ...
                'Color', colors(g,:), ...
                'HandleVisibility', 'off');
        end

        xlabel(ax, frequencyLabel('f', 'Frequency', options.Interpreter), ...
            labelArgs{:});
        ylabel(ax, responseAxisLabel(variables{k}, options), ...
            labelArgs{:});
        applyAxisWindow(ax, options.RawXLim{k}, options.RawXTicks{k}, ...
            options.RawYLim{k}, options.RawYTicks{k});

        handles.Colorbars(k) = addTemperatureColorbar(ax, colorMap, ...
            temperatureLimits, options.RawTempTicks, options);
    end
end

function handles = plotShiftedData(data, options)
% Plot shifted storage modulus, loss modulus, and loss factor together.

    requireTableVariables(data, {'ShiftedFrequency_Hz', 'StorageModulus', ...
        'LossModulus', 'TanDelta', 'Temperature_C'});

    [groupIndex, temperatures, colors, colorMap, temperatureLimits] = ...
        temperatureGroups(data, options.ShiftedTempLimits, options);

    handles             = initializeHandles(1);
    handles.Figures(1)  = figure('Name', 'Shifted DMA', ...
        'NumberTitle', 'off', 'Color', 'w');
    handles.Axes(1)     = axes('Parent', handles.Figures(1));
    ax                  = handles.Axes(1);
    [xScale, yScale]    = splitScale(options.ShiftedScale);
    labelArgs           = labelArguments(options);

    yyaxis(ax, 'left');
    prepareAxes(ax, xScale, yScale, options);
    for g = 1:numel(temperatures)
        rows            = groupIndex == g;
        [frequency, ix] = sort(data.ShiftedFrequency_Hz(rows));
        storage         = data.StorageModulus(rows);
        loss            = data.LossModulus(rows);
        storage         = storage(ix);
        loss            = loss(ix);

        plotResponse(ax, frequency, storage, '-', colors(g,:), ...
            options.ShiftedMarker, options.ShiftedMarkerSize, ...
            options.ShiftedLineWidth, xScale, yScale);
        plotResponse(ax, frequency, loss, '--', colors(g,:), ...
            options.ShiftedMarker, options.ShiftedMarkerSize, ...
            options.ShiftedLineWidth, xScale, yScale);
    end
    ylabel(ax, combinedModulusLabel(options), labelArgs{:});
    applyWindow(ax, 'y', options.ShiftedYLim, options.ShiftedYTicks);

    yyaxis(ax, 'right');
    prepareAxes(ax, xScale, options.ShiftedSecondaryYScale, options);
    for g = 1:numel(temperatures)
        rows            = groupIndex == g;
        [frequency, ix] = sort(data.ShiftedFrequency_Hz(rows));
        lossFactor      = data.TanDelta(rows);
        lossFactor      = lossFactor(ix);

        plotResponse(ax, frequency, lossFactor, ':', colors(g,:), ...
            options.ShiftedMarker, options.ShiftedMarkerSize, ...
            options.ShiftedLineWidth, xScale, options.ShiftedSecondaryYScale);
    end
    ylabel(ax, lossFactorLabel(options.Interpreter), labelArgs{:});
    applyWindow(ax, 'y', options.ShiftedSecondaryYLim, ...
        options.ShiftedSecondaryYTicks);

    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(2).Color = [0 0 0];
    yyaxis(ax, 'left');
    xlabel(ax, frequencyLabel('f_s', 'Shifted frequency', ...
        options.Interpreter), labelArgs{:});
    applyWindow(ax, 'x', options.ShiftedXLim, options.ShiftedXTicks);

    handles.Colorbars(1) = addTemperatureColorbar(ax, colorMap, ...
        temperatureLimits, options.ShiftedTempTicks, options);
    handles.Legends(1)   = shiftedLegend(ax, options);
end

function handles = plotFilteredFit(filteredData, denseResponse, options)
% Plot filtered master-curve points and the generalized-Maxwell response.

    requireTableVariables(filteredData, {'ShiftedFrequency_Hz', ...
        'StorageModulusFiltered', 'LossModulusFiltered', 'TanDeltaFiltered'});
    requireTableVariables(denseResponse, {'StorageModulus', ...
        'LossModulus', 'TanDelta'});

    if ismember('ReducedFrequency_Hz', denseResponse.Properties.VariableNames)
        fitFrequency = denseResponse.ReducedFrequency_Hz;
    elseif ismember('Frequency_Hz', denseResponse.Properties.VariableNames)
        fitFrequency = denseResponse.Frequency_Hz;
    else
        error('dmaPlot:MissingVariable', ...
            'The Maxwell response needs Frequency_Hz or ReducedFrequency_Hz.');
    end

    handles             = initializeHandles(1);
    handles.Figures(1)  = figure('Name', ...
        'Filtered master curve and Maxwell fit', ...
        'NumberTitle', 'off', 'Color', 'w');
    handles.Axes(1)     = axes('Parent', handles.Figures(1));
    ax                  = handles.Axes(1);
    colors              = options.ResponseColors;
    [xScale, yScale]    = splitScale(options.FitScale);
    labelArgs           = labelArguments(options);

    yyaxis(ax, 'left');
    prepareAxes(ax, xScale, yScale, options);
    plotResponse(ax, filteredData.ShiftedFrequency_Hz, ...
        filteredData.StorageModulusFiltered, 'none', colors(1,:), ...
        options.FilteredMarker, options.FilteredMarkerSize, ...
        options.FilteredMarkerLineWidth, xScale, yScale);
    plotResponse(ax, filteredData.ShiftedFrequency_Hz, ...
        filteredData.LossModulusFiltered, 'none', colors(2,:), ...
        options.FilteredMarker, options.FilteredMarkerSize, ...
        options.FilteredMarkerLineWidth, xScale, yScale);
    plotResponse(ax, fitFrequency, denseResponse.StorageModulus, '-', ...
        colors(1,:), 'none', 1, options.MaxwellLineWidth, xScale, yScale);
    plotResponse(ax, fitFrequency, denseResponse.LossModulus, '-', ...
        colors(2,:), 'none', 1, options.MaxwellLineWidth, xScale, yScale);
    ylabel(ax, combinedModulusLabel(options), labelArgs{:});
    applyWindow(ax, 'y', options.FitYLim, options.FitYTicks);

    yyaxis(ax, 'right');
    prepareAxes(ax, xScale, options.FitSecondaryYScale, options);
    plotResponse(ax, filteredData.ShiftedFrequency_Hz, ...
        filteredData.TanDeltaFiltered, 'none', colors(3,:), ...
        options.FilteredMarker, options.FilteredMarkerSize, ...
        options.FilteredMarkerLineWidth, xScale, ...
        options.FitSecondaryYScale);
    plotResponse(ax, fitFrequency, denseResponse.TanDelta, '-', ...
        colors(3,:), 'none', 1, options.MaxwellLineWidth, xScale, ...
        options.FitSecondaryYScale);
    ylabel(ax, lossFactorLabel(options.Interpreter), labelArgs{:});
    applyWindow(ax, 'y', options.FitSecondaryYLim, options.FitSecondaryYTicks);

    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(2).Color = [0 0 0];
    yyaxis(ax, 'left');
    xlabel(ax, frequencyLabel('f_r', 'Reduced frequency', ...
        options.Interpreter), labelArgs{:});
    applyWindow(ax, 'x', options.FitXLim, options.FitXTicks);
    handles.Legends(1) = fitLegend(ax, options);
end

function handles = plotWLF(wlf, options)
% Plot measured shift factors and the fitted WLF curve.

    if ~isstruct(wlf) || ~isfield(wlf, 'Data') || ~istable(wlf.Data)
        error('dmaPlot:InvalidWLF', 'The WLF input must contain a Data table.');
    end
    requireTableVariables(wlf.Data, {'Temperature_C', ...
        'MeasuredLog10aT', 'FittedLog10aT'});

    [temperature, order] = sort(wlf.Data.Temperature_C);
    measured             = wlf.Data.MeasuredLog10aT(order);
    fitted               = wlf.Data.FittedLog10aT(order);

    handles             = initializeHandles(1);
    handles.Figures(1)  = figure('Name', 'WLF fit', ...
        'NumberTitle', 'off', 'Color', 'w');
    handles.Axes(1)     = axes('Parent', handles.Figures(1));
    ax                  = handles.Axes(1);
    labelArgs           = labelArguments(options);

    prepareAxes(ax, 'linear', 'linear', options);
    measuredHandle = plot(ax, temperature, measured, ...
        'LineStyle', 'none', ...
        'Marker', options.WLFMarker, ...
        'MarkerSize', options.WLFMarkerSize, ...
        'LineWidth', options.WLFMarkerLineWidth, ...
        'Color', options.ResponseColors(1,:));
    fittedHandle = plot(ax, temperature, fitted, '-', ...
        'LineWidth', options.WLFLineWidth, ...
        'Color', options.ResponseColors(2,:));

    xlabel(ax, temperatureAxisLabel(options.Interpreter), ...
        labelArgs{:});
    ylabel(ax, shiftFactorAxisLabel(options.Interpreter), ...
        labelArgs{:});
    applyAxisWindow(ax, options.WLFXLim, options.WLFXTicks, ...
        options.WLFYLim, options.WLFYTicks);

    handles.Legends(1) = legend(ax, [measuredHandle, fittedHandle], ...
        {'Calculated shift factors', 'WLF fit'}, ...
        'Location', options.WLFLegendLocation, ...
        'Interpreter', options.Interpreter, ...
        'FontName', options.FontName, ...
        'FontSize', options.LegendFontSize, ...
        'Box', options.LegendBox);
end

function options = defaultPlotOptions
% Keep the common settings first, followed by the four plot-specific blocks.

    options                           = struct;

    % Common appearance
    options.AxesLineWidth             = 1.0;
    options.FontName                  = 'Arial';
    options.AxisLabelFontSize         = 30;
    options.TickLabelFontSize         = 30;
    options.LegendFontSize            = 24;
    options.ColorbarFontSize          = 30;
    options.Interpreter               = 'none';
    options.Box                       = 'on';
    options.Grid                      = 'on';
    options.MinorTicks                = 'on';
    options.LegendBox                 = 'off';
    options.ColorMap                  = 'turbo';
    options.ColorbarLocation          = 'eastoutside';
    options.ModulusType               = 'E';
    options.ModulusUnit               = 'MPa';
    options.ResponseColors            = [0.0000 0.4470 0.7410; ...
                                         0.8500 0.3250 0.0980; ...
                                         0.1500 0.1500 0.1500];

    % Raw storage, loss, and loss-factor figures
    options.RawScale                  = {'loglog', 'loglog', 'semilogx'};
    options.RawLineWidth              = 3;
    options.RawMarker                 = 'none';
    options.RawMarkerSize             = 7;
    options.RawTempLimits             = [];
    options.RawTempTicks              = [];
    options.RawXLim                   = {[], [], []};
    options.RawXTicks                 = {[], [], []};
    options.RawYLim                   = {[], [], []};
    options.RawYTicks                 = {[], [], []};

    % Shifted master curve
    options.ShiftedScale              = 'loglog';
    options.ShiftedSecondaryYScale    = 'linear';
    options.ShiftedLineWidth          = 3;
    options.ShiftedMarker             = 'none';
    options.ShiftedMarkerSize         = 7;
    options.ShiftedLegendLocation     = 'northwest';
    options.ShiftedTempLimits         = [];
    options.ShiftedTempTicks          = [];
    options.ShiftedXLim               = [];
    options.ShiftedXTicks             = [];
    options.ShiftedYLim               = [];
    options.ShiftedYTicks             = [];
    options.ShiftedSecondaryYLim      = [];
    options.ShiftedSecondaryYTicks    = [];

    % Filtered master curve and Maxwell fit
    options.FitScale                  = 'loglog';
    options.FitSecondaryYScale        = 'linear';
    options.FilteredMarker            = 'o';
    options.FilteredMarkerSize        = 7;
    options.FilteredMarkerLineWidth   = 3;
    options.MaxwellLineWidth          = 3;
    options.FitLegendLocation         = 'east';
    options.FitXLim                   = [];
    options.FitXTicks                 = [];
    options.FitYLim                   = [];
    options.FitYTicks                 = [];
    options.FitSecondaryYLim         = [];
    options.FitSecondaryYTicks       = [];

    % WLF fit
    options.WLFMarker                 = 'o';
    options.WLFMarkerSize             = 8;
    options.WLFMarkerLineWidth        = 3;
    options.WLFLineWidth              = 3;
    options.WLFLegendLocation         = 'northeast';
    options.WLFXLim                   = [];
    options.WLFXTicks                 = [];
    options.WLFYLim                   = [];
    options.WLFYTicks                 = [];
end

function options = resolvePlotOptions(userOptions)
% Fill omitted fields without introducing another public options structure.

    options  = defaultPlotOptions;
    provided = fieldnames(userOptions);
    unknown  = setdiff(provided, fieldnames(options));

    if ~isempty(unknown)
        error('dmaPlot:UnknownPlotOption', ...
            'Unknown PlotOptions field: %s', strjoin(unknown, ', '));
    end

    for i = 1:numel(provided)
        options.(provided{i}) = userOptions.(provided{i});
    end

end

function [groupIndex, temperatures, colors, colorMap, limits] = ...
    temperatureGroups(data, requestedLimits, options)

    if ismember('Set', data.Properties.VariableNames)
        [groupIndex, groups] = findgroups(data.Set);
    else
        [groupIndex, groups] = findgroups(data.Temperature_C);
    end

    temperatures = nan(numel(groups), 1);
    for g = 1:numel(groups)
        temperatures(g) = mean(data.Temperature_C(groupIndex == g), 'omitnan');
    end

    colorMap = resolveColorMap(options.ColorMap);
    limits   = requestedLimits;
    if isempty(limits)
        limits = [min(temperatures), max(temperatures)];
    end
    if limits(1) == limits(2)
        limits = limits + [-0.5 0.5];
    end

    position = min(max((temperatures - limits(1)) / diff(limits), 0), 1);
    colors   = interp1(linspace(0, 1, size(colorMap,1)), colorMap, ...
        position, 'linear');
end

function colorbarHandle = addTemperatureColorbar(ax, colorMap, limits, ticks, options)
    colormap(ax, colorMap);
    clim(ax, limits);
    colorbarHandle = colorbar(ax, 'Location', options.ColorbarLocation);
    set(colorbarHandle, ...
        'FontName', options.FontName, ...
        'FontSize', options.ColorbarFontSize, ...
        'TickLabelInterpreter', options.Interpreter);

    if ~isempty(ticks)
        colorbarHandle.Ticks = ticks;
    end

    colorbarHandle.Label.String      = temperatureAxisLabel(options.Interpreter);
    colorbarHandle.Label.Interpreter = options.Interpreter;
    colorbarHandle.Label.FontName    = options.FontName;
    colorbarHandle.Label.FontSize    = options.AxisLabelFontSize;
end

function legendHandle = shiftedLegend(ax, options)
    yyaxis(ax, 'left');
    storage = plot(ax, nan, nan, 'k-', ...
        'LineWidth', options.ShiftedLineWidth);
    loss = plot(ax, nan, nan, 'k--', ...
        'LineWidth', options.ShiftedLineWidth);
    lossFactor = plot(ax, nan, nan, 'k:', ...
        'LineWidth', options.ShiftedLineWidth);

    legendHandle = legend(ax, [storage, loss, lossFactor], ...
        {modulusSymbol(options.ModulusType, 1, options.Interpreter), ...
         modulusSymbol(options.ModulusType, 2, options.Interpreter), ...
         tanSymbol(options.Interpreter)}, ...
        'Location', options.ShiftedLegendLocation, ...
        'Interpreter', options.Interpreter, ...
        'FontName', options.FontName, ...
        'FontSize', options.LegendFontSize, ...
        'Box', options.LegendBox);
end

function legendHandle = fitLegend(ax, options)
    yyaxis(ax, 'left');
    filtered = plot(ax, nan, nan, ...
        'LineStyle', 'none', ...
        'Marker', options.FilteredMarker, ...
        'MarkerSize', options.FilteredMarkerSize, ...
        'LineWidth', options.FilteredMarkerLineWidth, ...
        'Color', [0 0 0]);
    maxwell = plot(ax, nan, nan, 'k-', ...
        'LineWidth', options.MaxwellLineWidth);

    legendHandle = legend(ax, [filtered, maxwell], ...
        {'Filtered master curve', 'Maxwell fit'}, ...
        'Location', options.FitLegendLocation, ...
        'Interpreter', options.Interpreter, ...
        'FontName', options.FontName, ...
        'FontSize', options.LegendFontSize, ...
        'Box', options.LegendBox);
end

function plotResponse(ax, x, y, lineStyle, color, marker, markerSize, ...
    lineWidth, xScale, yScale)

    valid = isfinite(x) & isfinite(y);
    if strcmp(xScale, 'log'), valid = valid & x > 0; end
    if strcmp(yScale, 'log'), valid = valid & y > 0; end
    plot(ax, x(valid), y(valid), ...
        'LineStyle', lineStyle, ...
        'Marker', marker, ...
        'MarkerSize', markerSize, ...
        'LineWidth', lineWidth, ...
        'Color', color, ...
        'HandleVisibility', 'off');
end

function prepareAxes(ax, xScale, yScale, options)
    set(ax, ...
        'XScale', xScale, ...
        'YScale', yScale, ...
        'FontName', options.FontName, ...
        'FontSize', options.TickLabelFontSize, ...
        'LineWidth', options.AxesLineWidth, ...
        'TickLabelInterpreter', options.Interpreter, ...
        'XMinorTick', options.MinorTicks, ...
        'YMinorTick', options.MinorTicks);
    hold(ax, 'on');
    grid(ax, options.Grid);
    box(ax, options.Box);
end

function applyAxisWindow(ax, xLimits, xTicks, yLimits, yTicks)
    applyWindow(ax, 'x', xLimits, xTicks);
    applyWindow(ax, 'y', yLimits, yTicks);
end

function applyWindow(ax, axisName, limits, ticks)
    if strcmp(axisName, 'x')
        if ~isempty(limits), xlim(ax, limits); end
        if ~isempty(ticks), ax.XTick = ticks; end
    else
        if ~isempty(limits), ylim(ax, limits); end
        if ~isempty(ticks), ax.YTick = ticks; end
    end
end

function args = labelArguments(options)
    args = {'Interpreter', options.Interpreter, ...
        'FontName', options.FontName, ...
        'FontSize', options.AxisLabelFontSize};
end

function label = responseAxisLabel(name, options)
    switch name
        case 'StorageModulus'
            label = sprintf('Storage modulus, %s (%s)', ...
                modulusSymbol(options.ModulusType, 1, options.Interpreter), ...
                options.ModulusUnit);
        case 'LossModulus'
            label = sprintf('Loss modulus, %s (%s)', ...
                modulusSymbol(options.ModulusType, 2, options.Interpreter), ...
                options.ModulusUnit);
        otherwise
            label = lossFactorLabel(options.Interpreter);
    end
end

function label = combinedModulusLabel(options)
    label = sprintf('Storage and loss moduli, %s, %s (%s)', ...
        modulusSymbol(options.ModulusType, 1, options.Interpreter), ...
        modulusSymbol(options.ModulusType, 2, options.Interpreter), ...
        options.ModulusUnit);
end

function label = lossFactorLabel(interpreter)
    label = sprintf('Loss factor, %s', tanSymbol(interpreter));
end

function label = modulusSymbol(type, derivative, interpreter)
    type = upper(char(type));
    if strcmp(interpreter, 'latex')
        if derivative == 1
            label = ['$', type, '^{\prime}$'];
        else
            label = ['$', type, '^{\prime\prime}$'];
        end
    else
        label = [type, repmat(char(39), 1, derivative)];
    end
end

function label = tanSymbol(interpreter)
    if strcmp(interpreter, 'latex')
        label = '$\tan(\delta)$';
    elseif strcmp(interpreter, 'tex')
        label = 'tan(\delta)';
    else
        label = 'tan(delta)';
    end
end

function label = frequencyLabel(symbol, description, interpreter)
    if strcmp(interpreter, 'latex')
        label = sprintf('%s, $%s$ (Hz)', description, symbol);
    else
        label = sprintf('%s, %s (Hz)', description, symbol);
    end
end

function label = temperatureAxisLabel(interpreter)
    if strcmp(interpreter, 'latex')
        label   = 'Temperature ($^{\circ}$C)';
    elseif strcmp(interpreter, 'tex')
        label   = 'Temperature (^{\circ}C)';
    else
        label   = 'Temperature (deg C)';
    end
end

function label = shiftFactorAxisLabel(interpreter)
    if strcmp(interpreter, 'latex')
        label = '$\log_{10}(a_T)$';
    elseif strcmp(interpreter, 'tex')
        label = 'log_{10}(a_T)';
    else
        label = 'log10(a_T)';
    end
end

function [xScale, yScale] = splitScale(scale)
    switch lower(char(scale))
        case 'linear'
            xScale = 'linear';
            yScale = 'linear';
        case 'semilogx'
            xScale = 'log';
            yScale = 'linear';
        case 'semilogy'
            xScale = 'linear';
            yScale = 'log';
        case 'loglog'
            xScale = 'log';
            yScale = 'log';
    end
end

function map = resolveColorMap(input)
    if isnumeric(input)
        map     = input;
    else
        map     = feval(char(input), 256);
    end
end

function handles = initializeHandles(numberOfPlots)
    handles.Figures     = gobjects(numberOfPlots, 1);
    handles.Axes        = gobjects(numberOfPlots, 1);
    handles.Legends     = gobjects(numberOfPlots, 1);
    handles.Colorbars   = gobjects(numberOfPlots, 1);
end

function requireTableVariables(data, names)
    if ~istable(data)
        error('dmaPlot:InvalidData', 'DMA plot data must be a table.');
    end

    missing = names(~ismember(names, data.Properties.VariableNames));
    if ~isempty(missing)
        error('dmaPlot:MissingVariable', ...
            'Missing table variable: %s', strjoin(missing, ', '));
    end
end

function tf = isTextScalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end
