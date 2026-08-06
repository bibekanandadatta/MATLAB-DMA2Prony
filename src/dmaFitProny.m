function prony = dmaFitProny(masterCurve, numberOfTerms, varargin)
%DMAFITPRONY Fit a generalized Maxwell Prony series in frequency space.
%
% prony = dmaFitProny(masterCurve, numberOfTerms)
% prony = dmaFitProny(..., 'RelaxationTimeMode', 'manual', ...
%     'RelaxationTimeRange', [1e-8 1e5], ...
%     'FittingFrequencyRange', [1e-5 1e8])
%
% RelaxationTimeMode:
%   exact  - span tau = 1/(2*pi*f) for the selected data (default)
%   round  - use interior power-of-ten endpoints
%   manual - use RelaxationTimeRange exactly
%   min    - jointly optimize g_i and tau_i after the initial fit
%
% MASTERCURVE must contain ShiftedFrequency_Hz and storage/loss modulus
% columns. Filtered modulus columns are preferred when present, unless their
% names are supplied explicitly with StorageVariable and LossVariable. This
% function does not smooth the data.
%
% NUMBEROFTERMS is chosen by the user. FittingFrequencyRange limits the data
% used by the fit, while RelaxationTimeRange controls the allowable mechanism
% times in manual mode. The returned structure contains the Prony table,
% fitted curves, modulus limits, and error measures.
%
% See also DMAFILTERMASTERCURVE, DMACOMPAREPRONYTERMS, DMAEVALUATEPRONY.

    p       = inputParser;
    p.addRequired('masterCurve', @istable);
    p.addRequired('numberOfTerms', ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && ...
        x >= 1 && x == round(x));
    p.addParameter('StorageVariable', 'auto', @isTextScalar);
    p.addParameter('LossVariable', 'auto', @isTextScalar);
    p.addParameter('RelaxationTimeMode', 'exact', @isTextScalar);
    p.addParameter('RelaxationTimeRange', [], @validPositiveRangeOrEmpty);
    p.addParameter('FittingFrequencyRange', [], @validPositiveRangeOrEmpty);
    p.parse(masterCurve, numberOfTerms, varargin{:});
    
    if exist('lsqnonneg', 'file') ~= 2
        error('dmaFitProny:OptimizationToolboxRequired', ...
            'dmaFitProny requires lsqnonneg from MATLAB Optimization Toolbox.');
    end
    
    [storageName, lossName] = selectModulusVariables(masterCurve, ...
        p.Results.StorageVariable, p.Results.LossVariable);
    required                = {'ShiftedFrequency_Hz', storageName, lossName};
    missing                 = required(~ismember(required, masterCurve.Properties.VariableNames));
    if ~isempty(missing)
        error('dmaFitProny:MissingVariables', ...
            'Missing required variables: %s', strjoin(missing, ', '));
    end
    
    frequency       = masterCurve.ShiftedFrequency_Hz(:);
    storage         = masterCurve.(storageName)(:);
    loss            = masterCurve.(lossName)(:);
    valid           = isfinite(frequency) & frequency > 0 & ...
        isfinite(storage) & storage > 0 & isfinite(loss) & loss >= 0;
    
    % Apply the requested data window before choosing relaxation times. This keeps
    % the automatic time range tied to the points that actually enter the fit.
    fitRange        = p.Results.FittingFrequencyRange;
    
    if ~isempty(fitRange)
        valid = valid & frequency >= fitRange(1) & frequency <= fitRange(2);
    end
    
    frequency           = frequency(valid);
    storage             = storage(valid);
    loss                = loss(valid);
    [frequency, order]  = sort(frequency);
    storage             = storage(order);
    loss                = loss(order);
    
    if numel(frequency) < max(4, numberOfTerms)
        error('dmaFitProny:TooFewPoints', ...
            ['At least max(4, numberOfTerms) valid points are required in the ', ...
             'selected fitting-frequency range.']);
    end
    
    if frequency(1) == frequency(end)
        error('dmaFitProny:FrequencyWindowTooNarrow', ...
            'The selected fitting data must span more than one frequency.');
    end
    
    % The Maxwell equations use angular frequency. The corresponding mechanism
    % time is tau = 1/omega = 1/(2*pi*f), even though input and plots use hertz.
    omega               = 2*pi*frequency;
    actualTauRange      = [1/omega(end), 1/omega(1)];

    mode                = validatestring(lower(char(p.Results.RelaxationTimeMode)), ...
        {'exact', 'round', 'manual', 'min'});
    [tauDescending, selectedTauRange]   = selectRelaxationTimes( ...
        actualTauRange, numberOfTerms, mode, p.Results.RelaxationTimeRange);
    if ~strcmp(mode, 'manual') && ~isempty(p.Results.RelaxationTimeRange)
        warning('dmaFitProny:UnusedRelaxationTimeRange', ...
            ['RelaxationTimeRange is used only in manual mode and was ignored ', ...
             'for RelaxationTimeMode="%s".'], mode);
    end
    
    mechanismFrequencyRange     = [1/(2*pi*selectedTauRange(2)), 1/(2*pi*selectedTauRange(1))];
    
    % A manual window is allowed to be narrower than the data, but warn because
    % the missing end mechanisms can bias the fitted plateau behavior.
    coverageFactor = 1 + 1e-12;
    if mechanismFrequencyRange(1) > frequency(1)*coverageFactor || ...
            mechanismFrequencyRange(2) < frequency(end)/coverageFactor
        warning('dmaFitProny:RelaxationWindowCoverage', ...
            ['The selected relaxation-time window does not cover the complete ', ...
             'fitting-frequency window under f_i = 1/(2*pi*tau_i).']);
    end
    
    % Solve for E_inf and every branch modulus together. Relative row scaling
    % prevents the much larger storage modulus from dominating the loss fit.
    [E_0, gDescending, residualNorm] = relativeAllDataFit( ...
        omega, storage, loss, tauDescending);
    optimizer                           = 'relative-error all-data nonnegative least squares';
    
    % sum(g_i) above one would imply a negative equilibrium modulus. Rescaling is
    % the same safeguard used by the reference implementation.
    if sum(gDescending) >= 1
        gDescending = gDescending/sum(gDescending);
    end
    if strcmp(mode, 'min')
        if exist('fmincon', 'file') ~= 2
            error('dmaFitProny:FminconRequired', ...
                ['RelaxationTimeMode="min" requires fmincon from MATLAB ', ...
                 'Optimization Toolbox.']);
        end
        [gDescending, tauDescending, residualNorm]  = minimizeTerms( ...
            frequency, storage, loss, E_0, gDescending, tauDescending, ...
            actualTauRange);
        optimizer       = 'relative-error joint g_i and tau_i minimization';
    end
    
    % The nonlinear solver also respects sum(g_i) <= 1, but keep this check here
    % for small numerical violations at the constraint boundary.
    if sum(gDescending) >= 1
        gDescending     = gDescending/sum(gDescending);
    end
    
    % Present terms from shortest to longest relaxation time. Internally they are
    % descending because that order is convenient for the response matrix.
    tau_i               = flipud(tauDescending(:));
    g_i                 = flipud(gDescending(:));
    [tau_i, termOrder]  = sort(tau_i);
    g_i                 = g_i(termOrder);
    frequency_i         = 1 ./ (2*pi*tau_i);
    E_i                 = E_0*g_i;
    E_inf               = E_0*(1 - sum(g_i));
    N                   = numberOfTerms;
    terms               = table((1:N).', tau_i, frequency_i, g_i, E_i, ...
                        'VariableNames', {'Term', 'tau_i_s', 'f_i_Hz', 'g_i', 'E_i'});
    
    omegaTau            = omega*tau_i.';
    denominator         = 1 + omegaTau.^2;
    
    % Each matrix column is one Maxwell branch evaluated at all fitted frequencies.
    storageFit          = E_inf + (omegaTau.^2 ./ denominator)*E_i;
    lossFit             = (omegaTau ./ denominator)*E_i;
    tanDeltaMeasured    = loss ./ storage;
    tanDeltaFit         = lossFit ./ storageFit;
    
    fittedCurve = table(frequency, omega, storage, loss, tanDeltaMeasured, ...
        storageFit, lossFit, tanDeltaFit, ...
        'VariableNames', {'Frequency_Hz', 'AngularFrequency_rad_s', ...
        'MeasuredStorageModulus', 'MeasuredLossModulus', 'MeasuredTanDelta', ...
        'FittedStorageModulus', 'FittedLossModulus', 'FittedTanDelta'});
    
    rmseStorage                             = sqrt(mean((storageFit - storage).^2));
    rmseLoss                                = sqrt(mean((lossFit - loss).^2));
    rmseTanDelta                            = sqrt(mean((tanDeltaFit - tanDeltaMeasured).^2));
    relativeFloor                           = sqrt(eps)*max([storage; loss]);
    rmsRelativeStorage                      = sqrt(mean(((storageFit - storage) ./ ...
        max(storage, relativeFloor)).^2));
    rmsRelativeLoss                         = sqrt(mean(((lossFit - loss) ./ ...
        max(loss, relativeFloor)).^2));
    tanDeltaFloor                           = sqrt(eps)*max([tanDeltaMeasured; 1]);
    rmsRelativeTanDelta                     = sqrt(mean(((tanDeltaFit - tanDeltaMeasured) ./ ...
        max(tanDeltaMeasured, tanDeltaFloor)).^2));
    prony = struct;
    prony.Terms = terms;
    prony.E_0 = E_0;
    prony.E_inf = E_inf;
    prony.NumberOfTerms                     = N;
    prony.ResidualNorm                      = residualNorm;
    prony.RMSE_StorageModulus               = rmseStorage;
    prony.RMSE_LossModulus                  = rmseLoss;
    prony.RMSE_TanDelta                     = rmseTanDelta;
    prony.NRMSE_StorageModulus              = normalizedRMSE(rmseStorage, storage);
    prony.NRMSE_LossModulus                 = normalizedRMSE(rmseLoss, loss);
    prony.RMSRelativeError_StorageModulus   = rmsRelativeStorage;
    prony.RMSRelativeError_LossModulus      = rmsRelativeLoss;
    prony.RMSRelativeError_TanDelta         = rmsRelativeTanDelta;
    prony.FrequencyRange_Hz                 = [frequency(1), frequency(end)];
    prony.RelaxationTimeRange_s             = [min(tau_i), max(tau_i)];
    prony.RelaxationTimeMode                = mode;
    prony.FittingFrequencyRangeRequested_Hz = fitRange;
    prony.StorageVariable                   = storageName;
    prony.LossVariable                      = lossName;
    prony.RelaxationTimeDefinition          = 'tau_i = 1/(2*pi*f_i)';
    prony.Method                            = optimizer;
    prony.Optimizer                         = optimizer;
    prony.FittedCurve                       = fittedCurve;
end

function [E_0, g, residualNorm] = relativeAllDataFit( ...
    omega, storage, loss, tau)
% Fit the dimensional equilibrium and branch moduli at every measured point.
    omegaTau        = omega(:)*tau(:).';
    denominator     = 1 + omegaTau.^2;
    storageBasis    = omegaTau.^2 ./ denominator;
    lossBasis       = omegaTau ./ denominator;
    
    numberOfPoints  = numel(omega);
    designMatrix    = [ones(numberOfPoints, 1), storageBasis; ...
                    zeros(numberOfPoints, 1), lossBasis];
    measured        = [storage(:); loss(:)];
    
    % Scale each equation by its measured magnitude. The small common floor only
    % matters if a measured loss value is zero or extremely close to zero.
    relativeFloor                   = sqrt(eps)*max(measured);
    scale                           = max(measured, relativeFloor);
    [coefficients, residualNorm]    = lsqnonneg( ...
        designMatrix ./ scale, measured ./ scale);
    
    E_inf   = coefficients(1);
    E_i     = coefficients(2:end);
    E_0     = E_inf + sum(E_i);
    if ~isfinite(E_0) || E_0 <= 0
        error('dmaFitProny:InvalidInstantaneousModulus', ...
            'The fitted instantaneous modulus is not positive.');
    end
    g   = E_i/E_0;
end

function [storageName, lossName] = selectModulusVariables(data, storageInput, lossInput)
% Prefer explicitly filtered data when available, but let the caller override
% either column independently.

    storageName     = char(storageInput);
    lossName        = char(lossInput);
    
    if strcmpi(storageName, 'auto')
        if ismember('StorageModulusFiltered', data.Properties.VariableNames)
            storageName = 'StorageModulusFiltered';
        else
            storageName = 'StorageModulus';
        end
    end
    if strcmpi(lossName, 'auto')
        if ismember('LossModulusFiltered', data.Properties.VariableNames)
            lossName = 'LossModulusFiltered';
        else
            lossName = 'LossModulus';
        end
    end

end

function [tau, range] = selectRelaxationTimes(actualRange, N, mode, manualRange)
% Return times in descending order for the response matrix.

    switch mode
        case 'exact'
            range   = actualRange;
            tau     = logspace(log10(range(2)), log10(range(1)), N).';

        case 'round'
            % Keep only complete decades inside the measured time range.
            range   = [10^ceil(log10(actualRange(1))), ...
                       10^floor(log10(actualRange(2)))];

            if range(1) > range(2)
                error('dmaFitProny:FrequencyWindowTooNarrow', ...
                    ['The selected frequency window does not contain a complete ', ...
                     'rounded relaxation-time decade. Use exact or manual mode.']);
            end

            tau     = logspace(log10(range(2)), log10(range(1)), N).';

        case 'manual'
            if isempty(manualRange)
                error('dmaFitProny:MissingRelaxationTimeRange', ...
                    ['RelaxationTimeRange must be supplied when ', ...
                     'RelaxationTimeMode is manual.']);
            end

            range   = manualRange(:).';
            tau     = logspace(log10(range(2)), log10(range(1)), N).';

        case 'min'
            range   = actualRange;

            % Interior starting points leave room for the optimizer to move each
            % mechanism toward either frequency boundary.
            allTau  = logspace(log10(range(2)), log10(range(1)), N + 2).';
            tau     = allTau(2:end-1);
    end
end

function [g, tau, residual] = minimizeTerms(frequency, storage, loss, E_0, g0, tau0, tauBounds)

    N       = numel(g0);
    
    % Optimize log10(tau) so relaxation times remain well scaled across decades.
    x0      = [min(max(g0(:), 0), 1); log10(tau0(:))];
    lower   = [zeros(N,1); repmat(log10(tauBounds(1)), N, 1)];
    upper   = [ones(N,1); repmat(log10(tauBounds(2)), N, 1)];
    A       = [ones(1,N), zeros(1,N)];
    b       = 1;
    
    % Storage and loss are normalized by E_0 so neither residual block is weighted
    % merely by the modulus unit used in the workbook.
    measured            = [storage(:); loss(:)]/E_0;

    relativeFloor       = sqrt(eps)*max(measured);
    scale               = max(measured, relativeFloor);
    objective           = @(x) sum(((normalizedFrequencyResponse( ...
        2*pi*frequency, x(1:N), 10.^x(N+1:end)) - measured) ./ scale).^2);
    options             = optimoptions('fmincon', 'Display', 'off', ...
                                'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 2e4, ...
                                'MaxIterations', 1000);
    [x, residual]       = fmincon(objective, x0, A, b, [], [], lower, upper, [], options);
    g                   = x(1:N);
    tau                 = 10.^x(N+1:end);
end

function response = normalizedFrequencyResponse(omega, g, tau)
    omegaTau        = omega(:)*tau(:).';
    denominator     = 1 + omegaTau.^2;
    storage         = 1 - sum(g) + (omegaTau.^2 ./ denominator)*g(:);
    loss            = (omegaTau ./ denominator)*g(:);
    response        = [storage; loss];
end

function value = normalizedRMSE(rmse, measured)
    scale       = max(measured) - min(measured);
    if scale > 0
        value   = rmse/scale;
    else
        value   = NaN;
    end
end

function tf = validPositiveRangeOrEmpty(x)
    tf      = isempty(x) || (isnumeric(x) && numel(x) == 2 && ...
              all(isfinite(x)) && all(x > 0) && x(1) < x(2));
end

function tf = isTextScalar(x)
    tf      = ischar(x) || (isstring(x) && isscalar(x));
end
