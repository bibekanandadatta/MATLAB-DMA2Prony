function prediction = dmaPredictTemperatureResponse( ...
    prony, wlf, temperatures, frequencies, varargin)
%DMAPREDICTTEMPERATURERESPONSE Predict DMA response using WLF and Prony data.
%
% prediction = dmaPredictTemperatureResponse(prony, wlf, ...
%     temperatures_C, frequencies_Hz)
% prediction = dmaPredictTemperatureResponse(..., ...
%     'AllowExtrapolation', true)
%
% TEMPERATURES and FREQUENCIES are vectors. Every temperature-frequency pair
% is evaluated, and the result is returned as a table with one row per pair.
%
% The package shift convention is:
%   reduced frequency = physical frequency * aT
%
% By default, predictions outside the WLF calibration-temperature range or
% Prony fitting-frequency range are returned as NaN and explicitly flagged.
% Set AllowExtrapolation to true to evaluate those points while retaining the
% range flags in the output table.
%
% See also DMAFITWLF, DMAFITPRONY, DMAEVALUATEPRONY.

    p           = inputParser;
    p.addRequired('prony', @isstruct);
    p.addRequired('wlf', @isstruct);
    p.addRequired('temperatures', ...
        @(x) isnumeric(x) && isvector(x) && ~isempty(x) && all(isfinite(x)));
    p.addRequired('frequencies', ...
        @(x) isnumeric(x) && isvector(x) && ~isempty(x) && ...
        all(isfinite(x)) && all(x > 0));
    p.addParameter('AllowExtrapolation', false, ...
        @(x) islogical(x) && isscalar(x));
    p.parse(prony, wlf, temperatures, frequencies, varargin{:});
    
    requiredWLF     = {'ReferenceTemperature_C', 'C1', 'C2_C', 'Data'};
    if ~all(isfield(wlf, requiredWLF)) || ~istable(wlf.Data) || ...
            ~ismember('Temperature_C', wlf.Data.Properties.VariableNames)
        error('dmaPredictTemperatureResponse:InvalidWLFInput', ...
            'wlf must be an output structure returned by dmaFitWLF.');
    end
    if ~isfield(prony, 'FrequencyRange_Hz')
        error('dmaPredictTemperatureResponse:InvalidPronyInput', ...
            'prony.FrequencyRange_Hz is required to check prediction coverage.');
    end
    
    % ndgrid forms the Cartesian product; callers do not have to assemble the
    % temperature-frequency pairs themselves.
    [temperatureGrid, frequencyGrid]    = ndgrid(temperatures(:), frequencies(:));
    temperature                         = temperatureGrid(:);
    frequency                           = frequencyGrid(:);
    log10aT                             = dmaWLF(temperature, wlf.ReferenceTemperature_C, wlf.C1, wlf.C2_C);
    shiftFactor                         = 10.^log10aT;
    reducedFrequency                    = frequency .* shiftFactor;
    
    wlfRange                            = [min(wlf.Data.Temperature_C), max(wlf.Data.Temperature_C)];
    withinWLFRange                      = temperature >= wlfRange(1) & temperature <= wlfRange(2);
    withinMasterRange                   = reducedFrequency >= prony.FrequencyRange_Hz(1) & ...
        reducedFrequency <= prony.FrequencyRange_Hz(2);
    validWLF = isfinite(log10aT) & isfinite(shiftFactor) & shiftFactor > 0;
    
    % Range checks are based on reduced frequency, because that is the frequency
    % at which the reference-temperature Prony model is evaluated.
    if p.Results.AllowExtrapolation
        evaluate    = validWLF;
    else
        evaluate    = validWLF & withinWLFRange & withinMasterRange;
    end
    
    storage         = nan(size(frequency));
    loss            = nan(size(frequency));
    tanDelta        = nan(size(frequency));

    if any(evaluate)
        % Evaluate only accepted points and leave rejected rows as NaN. The flags
        % below show which calibration boundary was crossed.
        response            = dmaEvaluateProny(prony, 'Frequencies', reducedFrequency(evaluate));
        storage(evaluate)   = response.StorageModulus;
        loss(evaluate)      = response.LossModulus;
        tanDelta(evaluate)  = response.TanDelta;
    end
    
    if any(~evaluate)
        if p.Results.AllowExtrapolation
            reason          = 'nonfinite WLF shift factors';
        else
            reason          = 'the calibrated WLF or fitted master-curve range';
        end
        warning('dmaPredictTemperatureResponse:OutsideCalibrationRange', ...
            '%d of %d requested predictions are outside %s.', ...
            nnz(~evaluate), numel(evaluate), reason);
    end
    
    prediction              = table(temperature, frequency, reducedFrequency, log10aT, ...
                                    shiftFactor, storage, loss, tanDelta, withinWLFRange, ...
                                    withinMasterRange, evaluate, ...
                                    'VariableNames', {'Temperature_C', 'Frequency_Hz', ...
                                    'ReducedFrequency_Hz', 'log10_aT', 'aT', 'StorageModulus', ...
                                    'LossModulus', 'TanDelta', 'WithinWLFFitRange', ...
                                    'WithinMasterCurveRange', 'PredictionEvaluated'});
end