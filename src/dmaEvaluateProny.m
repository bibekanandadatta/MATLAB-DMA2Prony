function response = dmaEvaluateProny(prony, varargin)
%DMAEVALUATEPRONY Evaluate a generalized Maxwell model on a dense grid.
%
% response  = dmaEvaluateProny(prony)
% response  = dmaEvaluateProny(prony, 'FrequencyRange', [1e-6 1e8], 'PointsPerDecade', 25)
% response  = dmaEvaluateProny(prony, 'Frequencies', [0.1 1 10])
%
% PRONY is the structure returned by dmaFitProny. Supply either an explicit
% frequency vector or a range and number of points per decade. If neither is
% given, the frequency range used for fitting is reused.
%
% The returned table contains storage modulus, loss modulus, tan delta, and
% the relaxation modulus evaluated at the companion time t = 1/f.
%
% See also DMAFITPRONY, DMAPREDICTTEMPERATURERESPONSE.

    p       = inputParser;
    p.addRequired('prony', @isstruct);
    p.addParameter('FrequencyRange', [], @validPositiveRangeOrEmpty);
    p.addParameter('Frequencies', [], ...
        @(x) isempty(x) || (isnumeric(x) && isvector(x) && ...
        all(isfinite(x)) && all(x > 0)));
    p.addParameter('PointsPerDecade', 25, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2);
    p.parse(prony, varargin{:});
    
    validateProny(prony);
    if ~isempty(p.Results.FrequencyRange) && ~isempty(p.Results.Frequencies)
        error('dmaEvaluateProny:ConflictingFrequencyInputs', ...
            'Specify FrequencyRange or Frequencies, not both.');
    end
    
    if ~isempty(p.Results.Frequencies)
        % Preserve the requested order when the caller supplies individual points.
        frequency       = p.Results.Frequencies(:);
    else
        frequencyRange  = p.Results.FrequencyRange;
        if isempty(frequencyRange)
            if ~isfield(prony, 'FrequencyRange_Hz')
                error('dmaEvaluateProny:MissingFrequencyRange', ...
                    'Supply FrequencyRange because prony.FrequencyRange_Hz is absent.');
            end
            frequencyRange  = prony.FrequencyRange_Hz;
        end
        decades             = log10(frequencyRange(2)/frequencyRange(1));
        numberOfPoints      = max(2, ceil(p.Results.PointsPerDecade*decades) + 1);
        frequency           = logspace(log10(frequencyRange(1)), ...
                                       log10(frequencyRange(2)), numberOfPoints).';
    end
    
    % Evaluate all Maxwell branches at once. Each column of omegaTau corresponds
    % to one Prony term.
    omega           = 2*pi*frequency;
    tau             = prony.Terms.tau_i_s(:).';
    branchModulus   = prony.Terms.E_i(:);
    omegaTau        = omega*tau;
    denominator     = 1 + omegaTau.^2;
    storage         = prony.E_inf + (omegaTau.^2 ./ denominator)*branchModulus;
    loss            = (omegaTau ./ denominator)*branchModulus;
    tanDelta        = loss ./ storage;
    
    % The paired time axis is a plotting convention, not a frequency-to-time
    % transform of the measured data.
    time            = 1 ./ frequency;
    relaxation      = prony.E_inf + exp(-time ./ tau)*branchModulus;
    
    response = table(frequency, omega, storage, loss, tanDelta, time, relaxation, ...
        'VariableNames', {'Frequency_Hz', 'AngularFrequency_rad_s', ...
        'StorageModulus', 'LossModulus', 'TanDelta', 'Time_s', ...
        'RelaxationModulus'});
end


function tf = validPositiveRangeOrEmpty(x)
    tf =     isempty(x) || (isnumeric(x) && numel(x) == 2 && ...
        all(isfinite(x)) && all(x > 0) && x(1) < x(2));
end


function validateProny(prony)
% Fail early here so matrix operations below do not produce cryptic errors.
    requiredFields  = {'Terms', 'E_inf'};
    
    if ~all(isfield(prony, requiredFields)) || ~istable(prony.Terms) || ...
            ~all(ismember({'tau_i_s', 'E_i'}, ...
            prony.Terms.Properties.VariableNames))
        error('dmaEvaluateProny:InvalidPronyInput', ...
            'prony must contain E_inf and a Terms table with tau_i_s and E_i.');
    end
    
    if any(~isfinite(prony.Terms.tau_i_s) | prony.Terms.tau_i_s <= 0) || ...
            any(~isfinite(prony.Terms.E_i) | prony.Terms.E_i < 0) || ...
            ~isfinite(prony.E_inf) || prony.E_inf < 0
        error('dmaEvaluateProny:InvalidPronyParameters', ...
            'Relaxation times and moduli must be finite and nonnegative.');
    end

end
