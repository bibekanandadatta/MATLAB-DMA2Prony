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
%   round  - use the interior power-of-ten endpoints used by PyVisco
%   manual - use RelaxationTimeRange exactly
%   min    - PyVisco-style joint optimization of g_i and tau_i
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

p = inputParser;
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
required = {'ShiftedFrequency_Hz', storageName, lossName};
missing = required(~ismember(required, masterCurve.Properties.VariableNames));
if ~isempty(missing)
    error('dmaFitProny:MissingVariables', ...
        'Missing required variables: %s', strjoin(missing, ', '));
end

frequency = masterCurve.ShiftedFrequency_Hz(:);
storage = masterCurve.(storageName)(:);
loss = masterCurve.(lossName)(:);
valid = isfinite(frequency) & frequency > 0 & ...
    isfinite(storage) & storage > 0 & isfinite(loss) & loss >= 0;

% Apply the requested data window before choosing relaxation times. This keeps
% the automatic time range tied to the points that actually enter the fit.
fitRange = p.Results.FittingFrequencyRange;
if ~isempty(fitRange)
    valid = valid & frequency >= fitRange(1) & frequency <= fitRange(2);
end
frequency = frequency(valid);
storage = storage(valid);
loss = loss(valid);
[frequency, order] = sort(frequency);
storage = storage(order);
loss = loss(order);

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
omega = 2*pi*frequency;
actualTauRange = [1/omega(end), 1/omega(1)];
mode = validatestring(lower(char(p.Results.RelaxationTimeMode)), ...
    {'exact', 'round', 'manual', 'min'});
[tauDescending, selectedTauRange] = selectRelaxationTimes( ...
    actualTauRange, numberOfTerms, mode, p.Results.RelaxationTimeRange);
if ~strcmp(mode, 'manual') && ~isempty(p.Results.RelaxationTimeRange)
    warning('dmaFitProny:UnusedRelaxationTimeRange', ...
        ['RelaxationTimeRange is used only in manual mode and was ignored ', ...
         'for RelaxationTimeMode="%s".'], mode);
end

mechanismFrequencyRange = [1/(2*pi*selectedTauRange(2)), ...
    1/(2*pi*selectedTauRange(1))];

% A manual window is allowed to be narrower than the data, but warn because
% the missing end mechanisms can bias the fitted plateau behavior.
coverageFactor = 1 + 1e-12;
if mechanismFrequencyRange(1) > frequency(1)*coverageFactor || ...
        mechanismFrequencyRange(2) < frequency(end)/coverageFactor
    warning('dmaFitProny:RelaxationWindowCoverage', ...
        ['The selected relaxation-time window does not cover the complete ', ...
         'fitting-frequency window under f_i = 1/(2*pi*tau_i).']);
end

% Use the measured high-frequency storage modulus as the instantaneous scale.
% The low-to-high change sets the scale for the collocation equations.
E_0 = storage(end);
E_inf_estimate = storage(1);
modulusRange = E_0 - E_inf_estimate;
if ~isfinite(modulusRange) || modulusRange <= 0
    error('dmaFitProny:InvalidModulusRange', ...
        'Storage modulus must increase from low to high fitting frequency.');
end

[gDescending, residualNorm] = generalizedCollocation( ...
    frequency, storage, loss, tauDescending, modulusRange);

% sum(g_i) above one would imply a negative equilibrium modulus. Rescaling is
% the same safeguard used by the reference implementation.
if sum(gDescending) >= 1
    gDescending = gDescending/sum(gDescending);
end
optimizer = 'nonnegative generalized collocation';

if strcmp(mode, 'min')
    if exist('fmincon', 'file') ~= 2
        error('dmaFitProny:FminconRequired', ...
            ['RelaxationTimeMode="min" requires fmincon from MATLAB ', ...
             'Optimization Toolbox.']);
    end
    [gDescending, tauDescending, residualNorm] = minimizeTerms( ...
        frequency, storage, loss, E_0, gDescending, tauDescending, ...
        actualTauRange);
    optimizer = 'PyVisco-style joint g_i and tau_i minimization';
end

% The nonlinear solver also respects sum(g_i) <= 1, but keep this check here
% for small numerical violations at the constraint boundary.
if sum(gDescending) >= 1
    gDescending = gDescending/sum(gDescending);
end

% Present terms from shortest to longest relaxation time. Internally they are
% descending because that order is convenient for the collocation matrix.
tau_i = flipud(tauDescending(:));
g_i = flipud(gDescending(:));
[tau_i, termOrder] = sort(tau_i);
g_i = g_i(termOrder);
frequency_i = 1 ./ (2*pi*tau_i);
E_i = E_0*g_i;
E_inf = E_0*(1 - sum(g_i));
N = numberOfTerms;
terms = table((1:N).', tau_i, frequency_i, g_i, E_i, ...
    'VariableNames', {'Term', 'tau_i_s', 'f_i_Hz', 'g_i', 'E_i'});

omegaTau = omega*tau_i.';
denominator = 1 + omegaTau.^2;

% Each matrix column is one Maxwell branch evaluated at all fitted frequencies.
storageFit = E_inf + (omegaTau.^2 ./ denominator)*E_i;
lossFit = (omegaTau ./ denominator)*E_i;
tanDeltaMeasured = loss ./ storage;
tanDeltaFit = lossFit ./ storageFit;

fittedCurve = table(frequency, omega, storage, loss, tanDeltaMeasured, ...
    storageFit, lossFit, tanDeltaFit, ...
    'VariableNames', {'Frequency_Hz', 'AngularFrequency_rad_s', ...
    'MeasuredStorageModulus', 'MeasuredLossModulus', 'MeasuredTanDelta', ...
    'FittedStorageModulus', 'FittedLossModulus', 'FittedTanDelta'});

rmseStorage = sqrt(mean((storageFit - storage).^2));
rmseLoss = sqrt(mean((lossFit - loss).^2));
prony = struct;
prony.Terms = terms;
prony.E_0 = E_0;
prony.E_inf = E_inf;
prony.NumberOfTerms = N;
prony.ResidualNorm = residualNorm;
prony.RMSE_StorageModulus = rmseStorage;
prony.RMSE_LossModulus = rmseLoss;
prony.NRMSE_StorageModulus = normalizedRMSE(rmseStorage, storage);
prony.NRMSE_LossModulus = normalizedRMSE(rmseLoss, loss);
prony.FrequencyRange_Hz = [frequency(1), frequency(end)];
prony.RelaxationTimeRange_s = [min(tau_i), max(tau_i)];
prony.RelaxationTimeMode = mode;
prony.FittingFrequencyRangeRequested_Hz = fitRange;
prony.StorageVariable = storageName;
prony.LossVariable = lossName;
prony.RelaxationTimeDefinition = 'tau_i = 1/(2*pi*f_i)';
prony.Method = optimizer;
prony.Optimizer = optimizer;
prony.FittedCurve = fittedCurve;
end

function [storageName, lossName] = selectModulusVariables(data, storageInput, lossInput)
% Prefer explicitly filtered data when available, but let the caller override
% either column independently.
storageName = char(storageInput);
lossName = char(lossInput);
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
% Return times in descending order for generalizedCollocation.
switch mode
    case 'exact'
        range = actualRange;
        tau = logspace(log10(range(2)), log10(range(1)), N).';
    case 'round'
        % Keep only complete decades inside the measured time range.
        range = [10^ceil(log10(actualRange(1))), ...
            10^floor(log10(actualRange(2)))];
        if range(1) > range(2)
            error('dmaFitProny:FrequencyWindowTooNarrow', ...
                ['The selected frequency window does not contain a complete ', ...
                 'rounded relaxation-time decade. Use exact or manual mode.']);
        end
        tau = logspace(log10(range(2)), log10(range(1)), N).';
    case 'manual'
        if isempty(manualRange)
            error('dmaFitProny:MissingRelaxationTimeRange', ...
                ['RelaxationTimeRange must be supplied when ', ...
                 'RelaxationTimeMode is manual.']);
        end
        range = manualRange(:).';
        tau = logspace(log10(range(2)), log10(range(1)), N).';
    case 'min'
        range = actualRange;
        % Interior starting points leave room for the optimizer to move each
        % mechanism toward either frequency boundary.
        allTau = logspace(log10(range(2)), log10(range(1)), N + 2).';
        tau = allTau(2:end-1);
end
end

function [g, residualNorm] = generalizedCollocation( ...
    frequency, storage, loss, tau, modulusRange)
% Sample the master curve at each mechanism's characteristic frequency.
frequency_i = 1 ./ (2*pi*tau);
storage_i = interpolateWithEndpointValues(frequency, storage, frequency_i);
loss_i = interpolateWithEndpointValues(frequency, loss, frequency_i);

N = numel(tau);

% At its own characteristic frequency a Maxwell branch contributes half of
% its stiffness to storage and reaches half of its peak loss contribution.
% The small neighboring loss bands account for nearby mechanisms without
% making the initial nonnegative solve overly dense.
K_storage = tril(ones(N), -1) + 0.5*eye(N);
K_loss = 0.5*eye(N);
bandWeights = [0.1, 0.01, 0.001];
for offset = 1:min(3, N-1)
    band = diag(ones(N-offset, 1), offset);
    K_loss = K_loss + bandWeights(offset)*(band + band.');
end
K = [K_storage; K_loss; ones(1, N)];

% The last row encourages the fitted branches to span the measured modulus
% change. lsqnonneg keeps every branch modulus nonnegative.
target = [storage_i/modulusRange; loss_i/modulusRange; 1];
[g, residualNorm] = lsqnonneg(K, target);
end

function values = interpolateWithEndpointValues(x, y, query)
% numpy.interp, used by PyVisco, holds endpoint values outside the data.
[xUnique, ~, group] = unique(x, 'sorted');
if numel(xUnique) < numel(x)
    % Shifted isotherms may contribute several measurements at one frequency.
    % Average duplicates before calling interp1.
    y = accumarray(group, y, [], @mean);
    x = xUnique;
end
queryInside = min(max(query, x(1)), x(end));
values = interp1(x, y, queryInside, 'linear');
end

function [g, tau, residual] = minimizeTerms( ...
    frequency, storage, loss, E_0, g0, tau0, tauBounds)
N = numel(g0);

% Optimize log10(tau) so relaxation times remain well scaled across decades.
x0 = [min(max(g0(:), 0), 1); log10(tau0(:))];
lower = [zeros(N,1); repmat(log10(tauBounds(1)), N, 1)];
upper = [ones(N,1); repmat(log10(tauBounds(2)), N, 1)];
A = [ones(1,N), zeros(1,N)];
b = 1;

% Storage and loss are normalized by E_0 so neither residual block is weighted
% merely by the modulus unit used in the workbook.
measured = [storage(:); loss(:)]/E_0;
objective = @(x) sum((normalizedFrequencyResponse( ...
    2*pi*frequency, x(1:N), 10.^x(N+1:end)) - measured).^2);
options = optimoptions('fmincon', 'Display', 'off', ...
    'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 2e4, ...
    'MaxIterations', 1000);
[x, residual] = fmincon(objective, x0, A, b, [], [], ...
    lower, upper, [], options);
g = x(1:N);
tau = 10.^x(N+1:end);
end

function response = normalizedFrequencyResponse(omega, g, tau)
omegaTau = omega(:)*tau(:).';
denominator = 1 + omegaTau.^2;
storage = 1 - sum(g) + (omegaTau.^2 ./ denominator)*g(:);
loss = (omegaTau ./ denominator)*g(:);
response = [storage; loss];
end

function value = normalizedRMSE(rmse, measured)
scale = max(measured) - min(measured);
if scale > 0
    value = rmse/scale;
else
    value = NaN;
end
end

function tf = validPositiveRangeOrEmpty(x)
tf = isempty(x) || (isnumeric(x) && numel(x) == 2 && ...
    all(isfinite(x)) && all(x > 0) && x(1) < x(2));
end

function tf = isTextScalar(x)
tf = ischar(x) || (isstring(x) && isscalar(x));
end
