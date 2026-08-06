function [ax, plotHandles] = dmaPlot(data, varargin)
%DMAPLOT Plot DMA signals with optional grouping and a secondary y-axis.
%
% plotOptions = dmaPlot('defaults')
% dmaPlot(data)
% dmaPlot(data, 'XVariable', 'Frequency_Hz', ...
%     'YVariables', {'StorageModulus','LossModulus'}, ...
%     'SecondaryYVariables', {'TanDelta'}, ...
%     'GroupVariable', 'Set', 'PlotOptions', plotOptions)
%
% Grouped isotherms can be colored continuously by temperature. In that
% mode a colorbar replaces the long temperature legend. ModulusType and
% ModulusUnit in PlotOptions control labels only; values are not converted.
%
% PlotOptions may contain any subset of the default fields. Missing fields
% use their default values.
%
% See also DMAREADTAEXCEL.

    if nargin == 1 && isTextScalar(data) && strcmpi(char(data), 'defaults')
        ax              = defaultPlotOptions;
        plotHandles     = struct;
        return
    end

    p       = inputParser;
    p.addRequired('data', @istable);
    p.addParameter('XVariable', '', @isTextScalar);
    p.addParameter('YVariables', {}, @validVariableList);
    p.addParameter('SecondaryYVariables', {}, @validVariableList);
    p.addParameter('GroupVariable', 'auto', @isTextScalar);
    p.addParameter('Scale', 'auto', @isTextScalar);
    p.addParameter('Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
    p.addParameter('PlotOptions', struct, @isstruct);
    p.parse(data, varargin{:});

    options     = resolvePlotOptions(p.Results.PlotOptions);
    xName       = char(p.Results.XVariable);

    if isempty(xName)
        candidates  = {'Frequency_Hz', 'Temperature_C', ...
                       'OscillationStrain_pct', 'Amplitude', 'Time_s'};
        xName       = firstAvailable(data, candidates);
    end

    if isempty(xName) || ~ismember(xName, data.Properties.VariableNames)
        error('dmaPlot:MissingXVariable', ...
            'The requested x variable is unavailable.');
    end

    primaryNames    = normalizeVariableList(p.Results.YVariables);
    secondaryNames  = normalizeVariableList(p.Results.SecondaryYVariables);

    if isempty(primaryNames) && isempty(secondaryNames)
        candidates      = {'StorageModulus', 'LossModulus', 'TanDelta'};
        primaryNames    = candidates(ismember(candidates, ...
            data.Properties.VariableNames));
    end

    allNames    = [primaryNames, secondaryNames];
    if isempty(allNames)
        error('dmaPlot:MissingYVariables', ...
            'No requested response variable is available.');
    end

    missing     = allNames(~ismember(allNames, data.Properties.VariableNames));
    if ~isempty(missing)
        error('dmaPlot:MissingYVariables', ...
            'Unavailable response variable: %s', strjoin(missing, ', '));
    end

    groupName   = char(p.Results.GroupVariable);
    if strcmpi(groupName, 'auto')
        if ismember('Set', data.Properties.VariableNames)
            groupName   = 'Set';
        else
            groupName   = '';
        end
    end

    [groupIndex, groupValues]             = makeGroups(data, groupName);
    [groupTemperatures, hasTemperature]   = findGroupTemperatures( ...
        data, groupIndex, groupValues, groupName);
    useTemperatureColorbar                = chooseTemperatureColorbar( ...
        options.TemperatureColorbar, hasTemperature, numel(groupValues));

    if useTemperatureColorbar
        [groupColors, colorMap, temperatureLimits] = temperatureColors( ...
            groupTemperatures, options);
    else
        groupColors        = lines(max(numel(groupValues), 1));
        colorMap           = [];
        temperatureLimits  = [];
    end

    if ismember('Scale', p.UsingDefaults)
        scale   = lower(char(options.Scale));
    else
        scale   = lower(char(p.Results.Scale));
    end

    if strcmp(scale, 'auto')
        if any(strcmp(xName, ...
                {'Frequency_Hz', 'AngularFrequency_rad_s', ...
                 'ShiftedFrequency_Hz', 'ReducedFrequency_Hz', 'Time_s'}))
            scale   = 'loglog';
        else
            scale   = 'semilogy';
        end
    end

    [xScale, primaryYScale]    = splitScale(scale);
    secondaryYScale            = options.SecondaryYScale;
    validateLogLimits(options, xScale, primaryYScale, ...
        secondaryYScale, ~isempty(secondaryNames));

    if isempty(p.Results.Axes)
        ax  = gca;
    else
        ax  = p.Results.Axes;
    end

    if ~isempty(secondaryNames)
        yyaxis(ax, 'left');
    end

    prepareAxes(ax, xScale, primaryYScale, options);
    primaryHandles  = plotResponses(ax, data, xName, primaryNames, ...
        groupIndex, groupValues, groupColors, groupName, options, ...
        useTemperatureColorbar, 0, numel(allNames));
    applyYWindow(ax, options.YLim, options.YTicks);
    ylabel(ax, makeAxisLabel(primaryNames, options), ...
        'Interpreter', options.Interpreter, ...
        'FontName', options.FontName, ...
        'FontSize', options.AxisLabelFontSize);

    secondaryHandles = gobjects(0,1);
    if ~isempty(secondaryNames)
        yyaxis(ax, 'right');
        prepareAxes(ax, xScale, secondaryYScale, options);
        secondaryHandles = plotResponses(ax, data, xName, secondaryNames, ...
            groupIndex, groupValues, groupColors, groupName, options, ...
            useTemperatureColorbar, numel(primaryNames), numel(allNames));
        applyYWindow(ax, options.SecondaryYLim, options.SecondaryYTicks);
        ylabel(ax, makeAxisLabel(secondaryNames, options), ...
            'Interpreter', options.Interpreter, ...
            'FontName', options.FontName, ...
            'FontSize', options.AxisLabelFontSize);

        ax.YAxis(1).Color  = [0 0 0];
        ax.YAxis(2).Color  = [0 0 0];
        yyaxis(ax, 'left');
    end

    applyXWindow(ax, options.XLim, options.XTicks);
    xlabel(ax, displayName(xName, options.Interpreter), ...
        'Interpreter', options.Interpreter, ...
        'FontName', options.FontName, ...
        'FontSize', options.AxisLabelFontSize);
    set(ax.Title, ...
        'Interpreter', options.Interpreter, ...
        'FontName', options.FontName, ...
        'FontSize', options.TitleFontSize);

    colorbarHandle = gobjects(0);
    if useTemperatureColorbar
        colormap(ax, colorMap);
        clim(ax, temperatureLimits);
        colorbarHandle = colorbar(ax, 'Location', options.ColorbarLocation);
        set(colorbarHandle, ...
            'FontName', options.FontName, ...
            'FontSize', options.ColorbarFontSize, ...
            'TickLabelInterpreter', options.Interpreter);

        if ~isempty(options.TemperatureTicks)
            colorbarHandle.Ticks = options.TemperatureTicks;
        end

        colorbarHandle.Label.String       = temperatureAxisLabel( ...
            options.Interpreter);
        colorbarHandle.Label.Interpreter  = options.Interpreter;
        colorbarHandle.Label.FontName     = options.FontName;
        colorbarHandle.Label.FontSize     = options.AxisLabelFontSize;
    end

    legendHandle = makeLegend(ax, primaryHandles, secondaryHandles, ...
        primaryNames, secondaryNames, options, useTemperatureColorbar);

    plotHandles            = struct;
    plotHandles.Primary    = primaryHandles;
    plotHandles.Secondary  = secondaryHandles;
    plotHandles.Legend     = legendHandle;
    plotHandles.Colorbar   = colorbarHandle;
end

function options = resolvePlotOptions(userOptions)
% Merge a partial user structure with the public defaults.

    options     = defaultPlotOptions;
    provided    = fieldnames(userOptions);
    unknown     = setdiff(provided, fieldnames(options));

    if ~isempty(unknown)
        error('dmaPlot:UnknownPlotOption', ...
            'Unknown PlotOptions field: %s', strjoin(unknown, ', '));
    end

    for i = 1:numel(provided)
        options.(provided{i}) = userOptions.(provided{i});
    end

    positiveFields  = {'LineWidth', 'MarkerSize', 'ShiftedLineWidth', ...
                       'FilteredMarkerSize', 'FilteredMarkerLineWidth', ...
                       'MaxwellLineWidth', 'AxesLineWidth', ...
                       'AxisLabelFontSize', 'TickLabelFontSize', ...
                       'LegendFontSize', 'TitleFontSize', 'ColorbarFontSize'};
    for i = 1:numel(positiveFields)
        field   = positiveFields{i};
        value   = options.(field);

        if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
            error('dmaPlot:InvalidPlotOption', ...
                'PlotOptions.%s must be a positive scalar.', field);
        end
    end

    textFields  = {'FontName', 'LegendLocation', 'ColorbarLocation', ...
                   'ShiftedLegendLocation', 'FitLegendLocation', ...
                   'Marker', 'ShiftedMarker', 'FilteredMarker', ...
                   'ModulusType', 'ModulusUnit'};
    for i = 1:numel(textFields)
        field = textFields{i};
        if ~isTextScalar(options.(field))
            error('dmaPlot:InvalidPlotOption', ...
                'PlotOptions.%s must be text.', field);
        end
        options.(field) = char(options.(field));
    end

    options.ModulusType = upper(options.ModulusType);
    if ~ismember(options.ModulusType, {'E', 'G'})
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.ModulusType must be E or G.');
    end

    options.Scale               = validateChoice(options.Scale, ...
        {'auto', 'linear', 'semilogx', 'semilogy', 'loglog'}, 'Scale');
    options.SecondaryYScale     = validateChoice(options.SecondaryYScale, ...
        {'linear', 'log'}, 'SecondaryYScale');
    options.TemperatureColorbar = validateChoice(options.TemperatureColorbar, ...
        {'auto', 'on', 'off'}, 'TemperatureColorbar');

    options.Interpreter = validateChoice(options.Interpreter, ...
        {'none', 'tex', 'latex'}, 'Interpreter');

    onOffFields    = {'Box', 'Grid', 'MinorTicks', 'ShowLegend', 'LegendBox'};
    for i = 1:numel(onOffFields)
        field               = onOffFields{i};
        options.(field)     = validateChoice(options.(field), ...
            {'on', 'off'}, field);
    end

    limitFields    = {'TemperatureLimits', 'XLim', 'YLim', 'SecondaryYLim'};
    for i = 1:numel(limitFields)
        field               = limitFields{i};
        options.(field)     = validateLimits(options.(field), field);
    end

    tickFields     = {'TemperatureTicks', 'XTicks', 'YTicks', 'SecondaryYTicks'};
    for i = 1:numel(tickFields)
        field               = tickFields{i};
        options.(field)     = validateTicks(options.(field), field);
    end

    if ~(iscell(options.ResponseLineStyles) && ...
            ~isempty(options.ResponseLineStyles) && ...
            all(cellfun(@isTextScalar, options.ResponseLineStyles)))
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.ResponseLineStyles must be a nonempty cell array of text.');
    end
    options.ResponseLineStyles = cellfun(@char, options.ResponseLineStyles, ...
        'UniformOutput', false);

    if ~(isnumeric(options.ResponseColors) && size(options.ResponseColors,2) == 3 && ...
            ~isempty(options.ResponseColors) && all(isfinite(options.ResponseColors), 'all') && ...
            all(options.ResponseColors >= 0, 'all') && ...
            all(options.ResponseColors <= 1, 'all'))
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.ResponseColors must be an N-by-3 array between zero and one.');
    end

    if ~(isTextScalar(options.ColorMap) || ...
            (isnumeric(options.ColorMap) && size(options.ColorMap,2) == 3 && ...
             ~isempty(options.ColorMap) && all(isfinite(options.ColorMap), 'all') && ...
             all(options.ColorMap >= 0, 'all') && all(options.ColorMap <= 1, 'all')))
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.ColorMap must be a colormap name or an N-by-3 array.');
    end
end

function options = defaultPlotOptions
% Keep every graphics default in one place for both public and internal use.

    options                           = struct;
    
    options.Scale                     = 'auto';
    options.SecondaryYScale           = 'linear';
    options.LineWidth                 = 3;
    options.Marker                    = 'none';
    options.MarkerSize                = 7;
    options.ShiftedLineWidth          = 3;
    options.ShiftedMarker             = 'none';
    options.FilteredMarker            = 'o';
    options.FilteredMarkerSize        = 7;
    options.FilteredMarkerLineWidth   = 3;
    options.MaxwellLineWidth          = 3;
    options.AxesLineWidth             = 1.0;
    options.FontName                  = 'Arial';
    options.AxisLabelFontSize         = 30;
    options.TickLabelFontSize         = 30;
    options.LegendFontSize            = 24;
    options.TitleFontSize             = 30;
    options.ColorbarFontSize          = 30;
    options.Interpreter               = 'none';
    options.Box                       = 'on';
    options.Grid                      = 'on';
    options.MinorTicks                = 'on';
    options.ShowLegend                = 'on';
    options.LegendLocation            = 'best';
    options.ShiftedLegendLocation     = 'northwest';
    options.FitLegendLocation         = 'east';
    options.LegendBox                 = 'off';
    options.TemperatureColorbar       = 'auto';
    options.ColorMap                  = 'turbo';
    options.ColorbarLocation          = 'eastoutside';
    options.ModulusType               = 'E';
    options.ModulusUnit               = 'MPa';
    
    options.ResponseLineStyles        = {'-', '--', ':', '-.'};
    options.ResponseColors            = [0.0000 0.4470 0.7410; ...
                                         0.8500 0.3250 0.0980; ...
                                         0.1500 0.1500 0.1500];

    options.TemperatureLimits         = [];
    options.TemperatureTicks          = [];
    options.XLim                      = [];
    options.YLim                      = [];
    options.SecondaryYLim             = [];
    options.XTicks                    = [];
    options.YTicks                    = [];
    options.SecondaryYTicks           = [];
end

function [groupIndex, groupValues] = makeGroups(data, groupName)
    if isempty(groupName)
        groupIndex     = ones(height(data), 1);
        groupValues    = 1;
        return
    end

    if ~ismember(groupName, data.Properties.VariableNames)
        error('dmaPlot:MissingGroupVariable', ...
            'Group variable "%s" is unavailable.', groupName);
    end

    [groupIndex, groupValues] = findgroups(data.(groupName));
end

function [temperatures, available] = findGroupTemperatures( ...
    data, groupIndex, groupValues, groupName)

    available       = false;
    temperatures    = nan(numel(groupValues), 1);

    if strcmp(groupName, 'Temperature_C') && isnumeric(groupValues)
        temperatures    = groupValues(:);
        available       = all(isfinite(temperatures));
    elseif ismember('Temperature_C', data.Properties.VariableNames)
        for g = 1:numel(groupValues)
            temperatures(g) = mean(data.Temperature_C(groupIndex == g), 'omitnan');
        end
        available       = all(isfinite(temperatures));
    end
end

function useColorbar = chooseTemperatureColorbar(choice, available, numberOfGroups)
    switch choice
        case 'auto'
            useColorbar = available && numberOfGroups > 1;
        case 'on'
            if ~available
                error('dmaPlot:TemperatureUnavailable', ...
                    ['TemperatureColorbar is on, but finite Temperature_C ', ...
                     'values are unavailable for every plotted group.']);
            end
            useColorbar = true;
        otherwise
            useColorbar = false;
    end
end

function [colors, map, limits] = temperatureColors(temperatures, options)
    map     = resolveColorMap(options.ColorMap);
    limits  = options.TemperatureLimits;

    if isempty(limits)
        limits  = [min(temperatures), max(temperatures)];
    end

    if limits(1) == limits(2)
        limits  = limits + [-0.5 0.5];
    end

    position    = (temperatures - limits(1)) / diff(limits);
    position    = min(max(position, 0), 1);
    mapPosition = linspace(0, 1, size(map,1));
    colors      = interp1(mapPosition, map, position, 'linear');
end

function map = resolveColorMap(input)
    if isnumeric(input)
        map = input;
        return
    end

    name        = lower(char(input));
    supported   = {'turbo', 'parula', 'jet', 'hot', 'cool', 'spring', ...
                   'summer', 'autumn', 'winter', 'gray'};
    if ~ismember(name, supported)
        error('dmaPlot:InvalidPlotOption', ...
            'Unsupported PlotOptions.ColorMap value "%s".', name);
    end
    map = feval(name, 256);
end

function [xScale, yScale] = splitScale(scale)
    switch scale
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
        otherwise
            error('dmaPlot:InvalidScale', 'Unknown plot scale "%s".', scale);
    end
end

function validateLogLimits(options, xScale, primaryScale, secondaryScale, hasSecondary)
    if strcmp(xScale, 'log') && ...
            ((~isempty(options.XLim) && any(options.XLim <= 0)) || ...
             (~isempty(options.XTicks) && any(options.XTicks <= 0)))
        error('dmaPlot:InvalidLogXAxis', ...
            'XLim and XTicks must be positive for a logarithmic x-axis.');
    end

    if strcmp(primaryScale, 'log') && ...
            ((~isempty(options.YLim) && any(options.YLim <= 0)) || ...
             (~isempty(options.YTicks) && any(options.YTicks <= 0)))
        error('dmaPlot:InvalidLogYAxis', ...
            'YLim and YTicks must be positive for a logarithmic primary y-axis.');
    end

    if hasSecondary && strcmp(secondaryScale, 'log') && ...
            ((~isempty(options.SecondaryYLim) && any(options.SecondaryYLim <= 0)) || ...
             (~isempty(options.SecondaryYTicks) && any(options.SecondaryYTicks <= 0)))
        error('dmaPlot:InvalidLogYAxis', ...
            ['SecondaryYLim and SecondaryYTicks must be positive for a ', ...
             'logarithmic secondary y-axis.']);
    end
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

function handles = plotResponses(ax, data, xName, names, ...
    groupIndex, groupValues, groupColors, groupName, options, ...
    useTemperatureColorbar, responseOffset, totalResponses)

    handles = gobjects(0,1);
    for g = 1:numel(groupValues)
        rows                    = groupIndex == g;
        [x, order]              = sort(data.(xName)(rows));
        rowIndices              = find(rows);
        rowIndices              = rowIndices(order);

        for j = 1:numel(names)
            responseNumber      = responseOffset + j;
            y                   = data.(names{j})(rowIndices);
            valid               = isfinite(x) & isfinite(y);

            if strcmp(ax.XScale, 'log')
                valid   = valid & x > 0;
            end
            if strcmp(ax.YScale, 'log')
                valid   = valid & y > 0;
            end
            if ~any(valid)
                continue
            end

            if useTemperatureColorbar || numel(groupValues) > 1
                color   = groupColors(g,:);
            else
                colorIndex  = 1 + mod(responseNumber-1, size(options.ResponseColors,1));
                color       = options.ResponseColors(colorIndex,:);
            end

            styleIndex  = 1 + mod(responseNumber-1, numel(options.ResponseLineStyles));
            style       = options.ResponseLineStyles{styleIndex};
            label       = makeDisplayLabel(data, rowIndices, groupName, ...
                groupValues(g), names{j}, options, totalResponses);

            handle = plot(ax, x(valid), y(valid), ...
                'LineStyle', style, ...
                'Marker', options.Marker, ...
                'MarkerSize', options.MarkerSize, ...
                'Color', color, ...
                'LineWidth', options.LineWidth, ...
                'DisplayName', label);

            if useTemperatureColorbar
                handle.HandleVisibility = 'off';
            end
            handles(end+1,1) = handle; %#ok<AGROW>
        end
    end
end

function applyXWindow(ax, limits, ticks)
    if ~isempty(limits)
        xlim(ax, limits);
    end
    if ~isempty(ticks)
        ax.XTick = ticks;
    end
end

function applyYWindow(ax, limits, ticks)
    if ~isempty(limits)
        ylim(ax, limits);
    end
    if ~isempty(ticks)
        ax.YTick = ticks;
    end
end

function legendHandle = makeLegend(ax, primaryHandles, secondaryHandles, ...
    primaryNames, secondaryNames, options, useTemperatureColorbar)

    legendHandle    = gobjects(0);
    if strcmp(options.ShowLegend, 'off')
        legend(ax, 'off');
        return
    end

    names       = [primaryNames, secondaryNames];
    dataHandles = [primaryHandles; secondaryHandles];

    if useTemperatureColorbar
        if isscalar(names)
            legend(ax, 'off');
            return
        end

        legendHandles   = gobjects(numel(names), 1);
        labels          = cell(numel(names), 1);
        for i = 1:numel(names)
            styleIndex      = 1 + mod(i-1, numel(options.ResponseLineStyles));
            legendHandles(i) = plot(ax, nan, nan, ...
                'LineStyle', options.ResponseLineStyles{styleIndex}, ...
                'Color', [0.15 0.15 0.15], ...
                'LineWidth', options.LineWidth, ...
                'DisplayName', responseSymbol(names{i}, options, ...
                    options.Interpreter));
            labels{i}       = legendHandles(i).DisplayName;
        end
    else
        legendHandles   = dataHandles;
        labels          = get(legendHandles, {'DisplayName'});
    end

    if isempty(legendHandles)
        return
    end

    legendHandle = legend(ax, legendHandles, labels, ...
        'Location', options.LegendLocation, ...
        'Interpreter', options.Interpreter);
    set(legendHandle, ...
        'FontName', options.FontName, ...
        'FontSize', options.LegendFontSize, ...
        'Box', options.LegendBox);
end

function label = makeDisplayLabel(data, rows, groupName, groupValue, ...
    yName, options, totalResponses)

    responseLabel   = responseSymbol(yName, options, options.Interpreter);
    if isempty(groupName)
        label   = responseLabel;
        return
    end

    groupLabel  = makeGroupLabel(data, rows, groupName, groupValue, ...
        options.Interpreter);
    if totalResponses > 1
        label   = sprintf('%s, %s', groupLabel, responseLabel);
    else
        label   = groupLabel;
    end
end

function label = makeGroupLabel(data, rows, groupName, groupValue, interpreter)
    if strcmp(groupName, 'Set') && ...
            ismember('Temperature_C', data.Properties.VariableNames)
        temperature = mean(data.Temperature_C(rows), 'omitnan');
        label       = formatTemperature(temperature, interpreter);
    elseif strcmp(groupName, 'Temperature_C')
        label       = formatTemperature(groupValue, interpreter);
    elseif strcmp(groupName, 'Frequency_Hz')
        if strcmp(interpreter, 'latex')
            label   = ['$', num2str(groupValue, '%.3g'), '\,\mathrm{Hz}$'];
        else
            label   = sprintf('%.3g Hz', groupValue);
        end
    elseif isnumeric(groupValue)
        label       = sprintf('%s = %g', groupName, groupValue);
    else
        label       = sprintf('%s = %s', groupName, char(string(groupValue)));
    end
end

function label = formatTemperature(temperature, interpreter)
    switch interpreter
        case 'latex'
            label   = ['$', num2str(temperature, '%.3g'), ...
                       '\,^{\circ}\mathrm{C}$'];
        case 'tex'
            label   = [num2str(temperature, '%.3g'), ' ^{\circ}C'];
        otherwise
            label   = sprintf('%.3g deg C', temperature);
    end
end

function label = makeAxisLabel(names, options)
    if isempty(names)
        label = '';
        return
    end

    hasStorage      = any(ismember(names, ...
        {'StorageModulus', 'StorageModulusFiltered'}));
    hasLoss         = any(ismember(names, ...
        {'LossModulus', 'LossModulusFiltered'}));
    hasRelaxation   = any(strcmp(names, 'RelaxationModulus'));
    hasLossFactor   = any(ismember(names, {'TanDelta', 'TanDeltaFiltered'}));

    labels = {};
    if hasStorage && hasLoss
        labels{end+1} = sprintf('Storage and loss moduli, %s, %s', ...
            modulusSymbol(options.ModulusType, 1, options.Interpreter), ...
            modulusSymbol(options.ModulusType, 2, options.Interpreter));
    elseif hasStorage
        labels{end+1} = sprintf('Storage modulus, %s', ...
            modulusSymbol(options.ModulusType, 1, options.Interpreter));
    elseif hasLoss
        labels{end+1} = sprintf('Loss modulus, %s', ...
            modulusSymbol(options.ModulusType, 2, options.Interpreter));
    elseif hasRelaxation
        labels{end+1} = sprintf('Relaxation modulus, %s', ...
            responseSymbol('RelaxationModulus', options, ...
                options.Interpreter));
    end

    hasModulus = hasStorage || hasLoss || hasRelaxation;
    if hasModulus && ~isempty(options.ModulusUnit)
        labels{end} = sprintf('%s (%s)', labels{end}, options.ModulusUnit);
    end

    if hasLossFactor
        labels{end+1} = sprintf('Loss factor, %s', ...
            responseSymbol('TanDelta', options, options.Interpreter));
    end

    knownNames = {'StorageModulus', 'StorageModulusFiltered', ...
                  'LossModulus', 'LossModulusFiltered', ...
                  'RelaxationModulus', 'TanDelta', 'TanDeltaFiltered'};
    otherNames = names(~ismember(names, knownNames));
    for i = 1:numel(otherNames)
        labels{end+1} = displayName(otherNames{i}, options.Interpreter); %#ok<AGROW>
    end

    label = strjoin(labels, '; ');
end

function label = responseSymbol(name, options, interpreter)
    if any(strcmp(name, {'StorageModulus', 'StorageModulusFiltered'}))
        label   = modulusSymbol(options.ModulusType, 1, interpreter);
    elseif any(strcmp(name, {'LossModulus', 'LossModulusFiltered'}))
        label   = modulusSymbol(options.ModulusType, 2, interpreter);
    elseif strcmp(name, 'RelaxationModulus')
        if strcmp(interpreter, 'latex')
            label   = sprintf('$%s(t)$', options.ModulusType);
        else
            label   = sprintf('%s(t)', options.ModulusType);
        end
    elseif any(strcmp(name, {'TanDelta', 'TanDeltaFiltered'}))
        if strcmp(interpreter, 'latex')
            label   = '$\tan(\delta)$';
        elseif strcmp(interpreter, 'tex')
            label   = 'tan(\delta)';
        else
            label   = 'tan(delta)';
        end
    else
        label   = displayName(name, interpreter);
    end
end

function label = modulusSymbol(type, derivative, interpreter)
    if strcmp(interpreter, 'latex')
        if derivative == 1
            label   = ['$', type, '^{\prime}$'];
        else
            label   = ['$', type, '^{\prime\prime}$'];
        end
    else
        label   = [type, repmat(char(39), 1, derivative)];
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

function name = firstAvailable(data, candidates)
    name = '';
    for i = 1:numel(candidates)
        if ismember(candidates{i}, data.Properties.VariableNames)
            name = candidates{i};
            return
        end
    end
end

function label = displayName(name, interpreter)
% Translate canonical table names into compact axis labels.

    if nargin < 2
        interpreter = 'none';
    end

    switch name
        case 'Frequency_Hz'
            label   = symbolAxisLabel('Frequency', 'f', 'Hz', interpreter);
        case 'AngularFrequency_rad_s'
            label   = symbolAxisLabel('Angular frequency', '\omega', ...
                'rad/s', interpreter);
        case 'ShiftedFrequency_Hz'
            label   = symbolAxisLabel('Shifted frequency', 'f_s', ...
                'Hz', interpreter);
        case 'ReducedFrequency_Hz'
            label   = symbolAxisLabel('Reduced frequency', 'f_r', ...
                'Hz', interpreter);
        case 'Temperature_C'
            label   = temperatureAxisLabel(interpreter);
        case 'OscillationStrain_pct'
            label   = 'Oscillation strain (%)';
        case 'Time_s'
            label   = symbolAxisLabel('Time', 't', 's', interpreter);
        otherwise
            label   = name;
    end
end

function label = symbolAxisLabel(description, symbol, unit, interpreter)
    if strcmp(interpreter, 'latex')
        label = sprintf('%s, $%s$ (%s)', description, symbol, unit);
    elseif strcmp(interpreter, 'none') && strcmp(symbol, '\omega')
        label = sprintf('%s, omega (%s)', description, unit);
    else
        label = sprintf('%s, %s (%s)', description, symbol, unit);
    end
end

function value = validateChoice(value, choices, fieldName)
    if ~isTextScalar(value)
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.%s must be text.', fieldName);
    end

    value = lower(char(value));
    if ~ismember(value, choices)
        error('dmaPlot:InvalidPlotOption', ...
            'Invalid value for PlotOptions.%s.', fieldName);
    end
end

function value = validateLimits(value, fieldName)
    if ~(isempty(value) || (isnumeric(value) && numel(value) == 2 && ...
            all(isfinite(value)) && value(1) < value(2)))
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.%s must be empty or an increasing two-value vector.', ...
            fieldName);
    end
    value = value(:).';
end

function value = validateTicks(value, fieldName)
    if ~(isempty(value) || (isnumeric(value) && isvector(value) && ...
            all(isfinite(value)) && all(diff(value) > 0)))
        error('dmaPlot:InvalidPlotOption', ...
            ['PlotOptions.%s must be empty or a strictly increasing ', ...
             'finite numeric vector.'], fieldName);
    end
    value = value(:).';
end

function names = normalizeVariableList(input)
    if isempty(input)
        names   = {};
    elseif ischar(input)
        names   = {input};
    else
        names   = cellstr(input);
    end
end

function tf = validVariableList(value)
    tf = isempty(value) || ischar(value) || isstring(value) || iscellstr(value);
end

function tf = isTextScalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end
