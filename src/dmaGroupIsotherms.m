function [data, isotherms] = dmaGroupIsotherms(data, varargin)
%DMAGROUPISOTHERMS Assign and summarize isothermal DMA frequency sweeps.
%
% [data, isotherms] = dmaGroupIsotherms(data)
% [data, isotherms] = dmaGroupIsotherms(data, 'TemperatureTolerance', 0.5)
%
% An existing Set column is used when available. Otherwise temperatures are
% rounded to TemperatureTolerance and used to form sets. ISOTHERMS contains
% the set identifier, mean measured temperature, and number of points in each
% sweep.
%
% See also DMAESTIMATESHIFTFACTORS, DMABUILDMASTERCURVE.

    p   = inputParser;
    p.addRequired('data', @istable);
    p.addParameter('TemperatureTolerance', 0.5, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.parse(data, varargin{:});
    
    if ~ismember('Temperature_C', data.Properties.VariableNames)
        error('dmaGroupIsotherms:MissingTemperature', ...
            'Temperature_C is required to identify isothermal sweeps.');
    end
    
    if ismember('Set', data.Properties.VariableNames) && ...
            all(isfinite(data.Set))
        % Respect set identifiers exported by the instrument or supplied by the
        % user. The labels do not have to be consecutive.
        [groupIndex, setValues] = findgroups(data.Set);
    else
        % Temperature readings usually wander slightly about the programmed value.
        % Rounding by the tolerance collects those readings into one isotherm.
        tol                     = p.Results.TemperatureTolerance;
        roundedTemperature      = round(data.Temperature_C / tol) * tol;
        [groupIndex, ~]         = findgroups(roundedTemperature);
        data.Set                = groupIndex;
        setValues               = (1:max(groupIndex)).';
    end
    
    % Keep the measured mean temperature rather than the rounded grouping value.
    temperature                 = splitapply(@(x) mean(x, 'omitnan'), data.Temperature_C, groupIndex);
    numberOfPoints              = splitapply(@numel, data.Temperature_C, groupIndex);
    
    isotherms                   = table(setValues, temperature, numberOfPoints, ...
                                  'VariableNames', {'Set', 'Temperature_C', 'NumberOfPoints'});
    isotherms                   = sortrows(isotherms, 'Temperature_C');
end
