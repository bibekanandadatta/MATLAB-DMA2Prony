function ax = dmaPlot(data, varargin)
%DMAPLOT Plot one or more DMA signals, optionally grouped by isotherm.
%
% dmaPlot(data)
% dmaPlot(data, 'XVariable', 'Frequency_Hz', ...
%     'YVariables', {'StorageModulus','LossModulus'}, ...
%     'GroupVariable', 'Set', 'Scale', 'loglog')
% dmaPlot(data, ..., 'PlotOptions', plotOptions)
%
% For amplitude sweeps, select OscillationStrain_pct or Amplitude as the
% XVariable and use Scale='semilogy' or 'linear'.
% GroupVariable='auto' uses Set when it is available. Pass an empty value to
% draw the table as one group. An existing axes can be supplied with Axes;
% otherwise the current axes is used. The returned value is the axes handle.
%
% PlotOptions fields:
%   Scale, LineWidth, AxesLineWidth, FontName, AxisLabelFontSize,
%   TickLabelFontSize, LegendFontSize, TitleFontSize,
%   LegendInterpreter, XLabelInterpreter, YLabelInterpreter,
%   TickLabelInterpreter, TitleInterpreter, Box, Grid, MinorTicks,
%   ShowLegend, LegendLocation, LegendBox, LegendNumColumns,
%   XLim, YLim, XTicks, YTicks.
%
% See also DMAREADTAEXCEL.

p = inputParser;
p.addRequired('data', @istable);
p.addParameter('XVariable', '', @(x) ischar(x) || isstring(x));
p.addParameter('YVariables', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
p.addParameter('GroupVariable', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('Scale', 'auto', @(x) any(strcmpi(string(x), ...
    ["auto","linear","semilogx","semilogy","loglog"])));
p.addParameter('Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
p.addParameter('PlotOptions', struct(), @isstruct);
p.parse(data, varargin{:});
plotOptions = resolvePlotOptions(p.Results.PlotOptions);

xName = char(p.Results.XVariable);
if isempty(xName)
    % Choose a useful abscissa for the common sweep types, in priority order.
    candidates = {'Frequency_Hz', 'Temperature_C', ...
        'OscillationStrain_pct', 'Amplitude', 'Time_s'};
    xName = firstAvailable(data, candidates);
end
if isempty(xName) || ~ismember(xName, data.Properties.VariableNames)
    error('dmaPlot:MissingXVariable', 'The requested x variable is unavailable.');
end

yNames = p.Results.YVariables;
if isempty(yNames)
    % With no explicit response list, draw the DMA quantities that are present.
    candidates = {'StorageModulus', 'LossModulus', 'TanDelta'};
    yNames = candidates(ismember(candidates, data.Properties.VariableNames));
elseif ischar(yNames)
    yNames = {yNames};
else
    yNames = cellstr(yNames);
end
if isempty(yNames)
    error('dmaPlot:MissingYVariables', 'No requested response variable is available.');
end
missing = yNames(~ismember(yNames, data.Properties.VariableNames));
if ~isempty(missing)
    error('dmaPlot:MissingYVariables', ...
        'Unavailable response variable: %s', strjoin(missing, ', '));
end

groupName = char(p.Results.GroupVariable);
if strcmpi(groupName, 'auto')
    if ismember('Set', data.Properties.VariableNames)
        groupName = 'Set';
    else
        groupName = '';
    end
end

if isempty(groupName)
    % A single artificial group keeps the plotting loop the same in both cases.
    groupIndex = ones(height(data), 1);
    groupValues = 1;
else
    if ~ismember(groupName, data.Properties.VariableNames)
        error('dmaPlot:MissingGroupVariable', ...
            'Group variable "%s" is unavailable.', groupName);
    end
    [groupIndex, groupValues] = findgroups(data.(groupName));
end

colors = lines(max(groupIndex));
styles = {'-', '--', ':', '-.'};
if ismember('Scale', p.UsingDefaults)
    scale = lower(char(plotOptions.Scale));
else
    % A direct name-value Scale input overrides PlotOptions.Scale.
    scale = lower(char(p.Results.Scale));
end
if strcmp(scale, 'auto')
    % Frequency and time usually span decades; other sweep coordinates are
    % easier to inspect on a linear x-axis.
    if any(strcmp(xName, {'Frequency_Hz', 'AngularFrequency_rad_s', 'Time_s'}))
        scale = 'loglog';
    else
        scale = 'semilogy';
    end
end
if any(strcmp(scale, {'semilogx', 'loglog'})) && ...
        ((~isempty(plotOptions.XLim) && any(plotOptions.XLim <= 0)) || ...
         (~isempty(plotOptions.XTicks) && any(plotOptions.XTicks <= 0)))
    error('dmaPlot:InvalidLogXAxis', ...
        'XLim and XTicks must be positive for a logarithmic x-axis.');
end
if any(strcmp(scale, {'semilogy', 'loglog'})) && ...
        ((~isempty(plotOptions.YLim) && any(plotOptions.YLim <= 0)) || ...
         (~isempty(plotOptions.YTicks) && any(plotOptions.YTicks <= 0)))
    error('dmaPlot:InvalidLogYAxis', ...
        'YLim and YTicks must be positive for a logarithmic y-axis.');
end

if isempty(p.Results.Axes)
    ax = gca;
else
    ax = p.Results.Axes;
end

% Set axis scales before plotting. MATLAB otherwise may discard nonpositive
% points while changing a populated axes from linear to logarithmic.
switch scale
    case 'linear'
        set(ax, 'XScale', 'linear', 'YScale', 'linear');
    case 'semilogx'
        set(ax, 'XScale', 'log', 'YScale', 'linear');
    case 'semilogy'
        set(ax, 'XScale', 'linear', 'YScale', 'log');
    case 'loglog'
        set(ax, 'XScale', 'log', 'YScale', 'log');
end
set(ax, ...
    'FontName', plotOptions.FontName, ...
    'FontSize', plotOptions.TickLabelFontSize, ...
    'LineWidth', plotOptions.AxesLineWidth, ...
    'TickLabelInterpreter', plotOptions.TickLabelInterpreter, ...
    'XMinorTick', plotOptions.MinorTicks, ...
    'YMinorTick', plotOptions.MinorTicks);
hold(ax, 'on');
grid(ax, plotOptions.Grid);
box(ax, plotOptions.Box);

for g = 1:numel(groupValues)
    rows = groupIndex == g;
    [x, order] = sort(data.(xName)(rows));
    rowIndices = find(rows);
    rowIndices = rowIndices(order);
    for j = 1:numel(yNames)
        y = data.(yNames{j})(rowIndices);
        valid = isfinite(x) & isfinite(y);

        % Do not pass invalid log coordinates to plot; leaving them out avoids
        % MATLAB warnings and broken lines at zero frequency.
        if any(strcmp(scale, {'semilogx', 'loglog'}))
            valid = valid & x > 0;
        end
        if any(strcmp(scale, {'semilogy', 'loglog'}))
            valid = valid & y > 0;
        end
        if ~any(valid)
            continue
        end

        label = makeLabel(data, rowIndices, groupName, groupValues(g), ...
            yNames{j}, plotOptions.LegendInterpreter);

        % Color identifies the sweep and line style distinguishes responses
        % when more than one y-variable is requested on the same axes.
        style = styles{1 + mod(j-1, numel(styles))};
        plot(ax, x(valid), y(valid), style, 'Color', colors(g,:), ...
            'LineWidth', plotOptions.LineWidth, ...
            'DisplayName', label);
    end
end

if ~isempty(plotOptions.XLim)
    xlim(ax, plotOptions.XLim);
end
if ~isempty(plotOptions.YLim)
    ylim(ax, plotOptions.YLim);
end
if ~isempty(plotOptions.XTicks)
    set(ax, 'XTick', plotOptions.XTicks);
end
if ~isempty(plotOptions.YTicks)
    set(ax, 'YTick', plotOptions.YTicks);
end

xlabel(ax, displayName(xName), ...
    'Interpreter', plotOptions.XLabelInterpreter, ...
    'FontName', plotOptions.FontName, ...
    'FontSize', plotOptions.AxisLabelFontSize);
ylabel(ax, strjoin(cellfun(@displayName, yNames, 'UniformOutput', false), ', '), ...
    'Interpreter', plotOptions.YLabelInterpreter, ...
    'FontName', plotOptions.FontName, ...
    'FontSize', plotOptions.AxisLabelFontSize);
set(ax.Title, ...
    'Interpreter', plotOptions.TitleInterpreter, ...
    'FontName', plotOptions.FontName, ...
    'FontSize', plotOptions.TitleFontSize);

if strcmp(plotOptions.ShowLegend, 'on')
    lgd = legend(ax, 'Location', plotOptions.LegendLocation, ...
        'Interpreter', plotOptions.LegendInterpreter);
    set(lgd, 'FontName', plotOptions.FontName, ...
        'FontSize', plotOptions.LegendFontSize, ...
        'Box', plotOptions.LegendBox, ...
        'NumColumns', plotOptions.LegendNumColumns);
else
    legend(ax, 'off');
end
end

function options = resolvePlotOptions(userOptions)
% Keep defaults in one structure so a test script can override only the few
% properties needed for a particular figure.
options = struct( ...
    'Scale', 'auto', ...
    'LineWidth', 1.5, ...
    'AxesLineWidth', 1.0, ...
    'FontName', 'Arial', ...
    'AxisLabelFontSize', 12, ...
    'TickLabelFontSize', 11, ...
    'LegendFontSize', 10, ...
    'TitleFontSize', 12, ...
    'LegendInterpreter', 'tex', ...
    'XLabelInterpreter', 'tex', ...
    'YLabelInterpreter', 'tex', ...
    'TickLabelInterpreter', 'tex', ...
    'TitleInterpreter', 'tex', ...
    'Box', 'on', ...
    'Grid', 'on', ...
    'MinorTicks', 'on', ...
    'ShowLegend', 'on', ...
    'LegendLocation', 'best', ...
    'LegendBox', 'on', ...
    'LegendNumColumns', 1, ...
    'XLim', [], ...
    'YLim', [], ...
    'XTicks', [], ...
    'YTicks', []);

provided = fieldnames(userOptions);
unknown = setdiff(provided, fieldnames(options));
if ~isempty(unknown)
    error('dmaPlot:UnknownPlotOption', ...
        'Unknown PlotOptions field: %s', strjoin(unknown, ', '));
end
for i = 1:numel(provided)
    options.(provided{i}) = userOptions.(provided{i});
end

% Validate after merging so defaults and user values follow the same path.
positiveFields = {'LineWidth', 'AxesLineWidth', 'AxisLabelFontSize', ...
    'TickLabelFontSize', 'LegendFontSize', 'TitleFontSize'};
for i = 1:numel(positiveFields)
    value = options.(positiveFields{i});
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.%s must be a positive scalar.', positiveFields{i});
    end
end

if ~(isnumeric(options.LegendNumColumns) && ...
        isscalar(options.LegendNumColumns) && ...
        isfinite(options.LegendNumColumns) && ...
        options.LegendNumColumns >= 1 && ...
        options.LegendNumColumns == round(options.LegendNumColumns))
    error('dmaPlot:InvalidPlotOption', ...
        'PlotOptions.LegendNumColumns must be a positive integer.');
end

if ~(ischar(options.FontName) || ...
        (isstring(options.FontName) && isscalar(options.FontName)))
    error('dmaPlot:InvalidPlotOption', ...
        'PlotOptions.FontName must be text.');
end
options.FontName = char(options.FontName);
if ~(ischar(options.LegendLocation) || ...
        (isstring(options.LegendLocation) && isscalar(options.LegendLocation)))
    error('dmaPlot:InvalidPlotOption', ...
        'PlotOptions.LegendLocation must be text.');
end
options.LegendLocation = char(options.LegendLocation);

limitFields = {'XLim', 'YLim'};
for i = 1:numel(limitFields)
    field = limitFields{i};
    value = options.(field);
    if ~(isempty(value) || (isnumeric(value) && numel(value) == 2 && ...
            all(isfinite(value)) && value(1) < value(2)))
        error('dmaPlot:InvalidPlotOption', ...
            'PlotOptions.%s must be empty or an increasing two-value vector.', ...
            field);
    end
    options.(field) = value(:).';
end

tickFields = {'XTicks', 'YTicks'};
for i = 1:numel(tickFields)
    field = tickFields{i};
    value = options.(field);
    if ~(isempty(value) || (isnumeric(value) && isvector(value) && ...
            all(isfinite(value)) && all(diff(value) > 0)))
        error('dmaPlot:InvalidPlotOption', ...
            ['PlotOptions.%s must be empty or a strictly increasing ', ...
             'finite numeric vector.'], field);
    end
    options.(field) = value(:).';
end

options.Scale = validateChoice(options.Scale, ...
    {'auto','linear','semilogx','semilogy','loglog'}, 'Scale');
interpreterFields = {'LegendInterpreter', 'XLabelInterpreter', ...
    'YLabelInterpreter', 'TickLabelInterpreter', 'TitleInterpreter'};
for i = 1:numel(interpreterFields)
    field = interpreterFields{i};
    options.(field) = validateChoice(options.(field), ...
        {'none','tex','latex'}, field);
end
onOffFields = {'Box', 'Grid', 'MinorTicks', 'ShowLegend', 'LegendBox'};
for i = 1:numel(onOffFields)
    field = onOffFields{i};
    options.(field) = validateChoice(options.(field), {'on','off'}, field);
end
end

function value = validateChoice(value, choices, fieldName)
% Normalize string and character inputs once for the graphics calls below.
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('dmaPlot:InvalidPlotOption', ...
        'PlotOptions.%s must be text.', fieldName);
end
value = lower(char(value));
if ~ismember(value, choices)
    error('dmaPlot:InvalidPlotOption', ...
        'Invalid value for PlotOptions.%s.', fieldName);
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

function label = makeLabel(data, rows, groupName, groupValue, yName, interpreter)
% Isotherm legends show temperature only. The response name is already clear
% from the y-axis when storage, loss, and tan delta use separate figures.
if strcmp(groupName, 'Set') && ...
        ismember('Temperature_C', data.Properties.VariableNames)
    temperature = mean(data.Temperature_C(rows), 'omitnan');
    label = formatTemperature(temperature, interpreter);
elseif strcmp(groupName, 'Temperature_C')
    label = formatTemperature(groupValue, interpreter);
elseif strcmp(groupName, 'Frequency_Hz')
    switch interpreter
        case 'latex'
            label = sprintf('$%.3g\\,\\mathrm{Hz}$', groupValue);
        case 'tex'
            label = sprintf('%.3g Hz', groupValue);
        otherwise
            label = sprintf('%.3g Hz', groupValue);
    end
elseif ~isempty(groupName)
    label = sprintf('%s = %g', groupName, groupValue);
else
    label = displayName(yName);
end
end

function label = formatTemperature(temperature, interpreter)
switch interpreter
    case 'latex'
        label = sprintf('$%.3g\\,^{\\circ}\\mathrm{C}$', temperature);
    case 'tex'
        label = sprintf('%.3g ^{\\circ}C', temperature);
    otherwise
        label = sprintf('%.3g deg C', temperature);
end
end

function label = displayName(name)
% Translate canonical table names into compact axis labels. Unknown user
% columns are left unchanged rather than guessed.
switch name
    case 'Frequency_Hz'
        label = 'Frequency (Hz)';
    case 'AngularFrequency_rad_s'
        label = 'Angular frequency (rad/s)';
    case 'ShiftedFrequency_Hz'
        label = 'Shifted frequency (Hz)';
    case 'ReducedFrequency_Hz'
        label = 'Reduced frequency (Hz)';
    case 'Temperature_C'
        label = 'Temperature (C)';
    case 'StorageModulus'
        label = 'Storage modulus';
    case 'LossModulus'
        label = 'Loss modulus';
    case 'TanDelta'
        label = 'tan(delta)';
    case 'StorageModulusFiltered'
        label = 'Filtered storage modulus';
    case 'LossModulusFiltered'
        label = 'Filtered loss modulus';
    case 'TanDeltaFiltered'
        label = 'Filtered tan(delta)';
    case 'RelaxationModulus'
        label = 'Relaxation modulus';
    case 'OscillationStrain_pct'
        label = 'Oscillation strain (%)';
    case 'Time_s'
        label = 'Time (s)';
    otherwise
        label = name;
end
end
