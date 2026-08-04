function [shiftFactors, diagnostics] = dmaEstimateShiftFactors(data, referenceTemperature, varargin)
%DMAESTIMATESHIFTFACTORS Estimate horizontal DMA shift factors.
%
% [shiftFactors, diagnostics] = dmaEstimateShiftFactors(data, RefT)
%
% This is a MATLAB implementation of the pairwise power-law method used by
% PyVisco. Neighboring storage-modulus isotherms are fitted with
% y = a*x^b + e. Their horizontal separation is averaged at common modulus
% levels, and incremental shifts are accumulated away from RefT.
%
% Required variables: Frequency_Hz, StorageModulus, Temperature_C.
% An existing Set variable is used; otherwise dmaGroupIsotherms creates one.
% ReferenceTolerance controls how close RefT must be to a measured isotherm.
% NumberOfLevels sets the number of common modulus levels used for each pair.
% When DropFirstPoint is true, the first point is kept or removed according to
% which power-law fit gives the smaller exponent standard error.
%
% SHIFTFACTORS is ready for dmaBuildMasterCurve. DIAGNOSTICS records every
% pairwise fit, including overlap and first-point decisions.
%
% See also DMAGROUPISOTHERMS, DMABUILDMASTERCURVE.

p = inputParser;
p.addRequired('data', @istable);
p.addRequired('referenceTemperature', ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('TemperatureTolerance', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ReferenceTolerance', 1.0, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('NumberOfLevels', 10, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 3 && x == round(x));
p.addParameter('DropFirstPoint', true, ...
    @(x) islogical(x) && isscalar(x));
p.parse(data, referenceTemperature, varargin{:});

required = {'Frequency_Hz', 'StorageModulus', 'Temperature_C'};
missing = required(~ismember(required, data.Properties.VariableNames));
if ~isempty(missing)
    error('dmaEstimateShiftFactors:MissingVariables', ...
        'Missing required variables: %s', strjoin(missing, ', '));
end

[data, isotherms] = dmaGroupIsotherms(data, ...
    'TemperatureTolerance', p.Results.TemperatureTolerance);

[referenceError, referenceIndex] = min(abs(isotherms.Temperature_C - referenceTemperature));
if referenceError > p.Results.ReferenceTolerance
    error('dmaEstimateShiftFactors:ReferenceNotMeasured', ...
        ['No measured isotherm is within %.3g C of the requested reference ', ...
         'temperature %.3g C.'], p.Results.ReferenceTolerance, referenceTemperature);
end

nSets = height(isotherms);
log10aT = nan(nSets, 1);
log10aT(referenceIndex) = 0;
diagnostics = repmat(emptyDiagnostic(), max(nSets-1, 0), 1);
d = 0;

% Work away from the reference in both temperature directions. Each local
% separation is added to the accumulated shift of the neighboring isotherm.
for i = referenceIndex:-1:2
    d = d + 1;
    refData = getSet(data, isotherms.Set(i));
    shiftData = getSet(data, isotherms.Set(i-1));
    [increment, diagnostics(d)] = fitShiftPair(refData, shiftData, ...
        p.Results.NumberOfLevels, p.Results.DropFirstPoint);
    diagnostics(d).ReferenceTemperature_C = isotherms.Temperature_C(i);
    diagnostics(d).ShiftedTemperature_C = isotherms.Temperature_C(i-1);
    log10aT(i-1) = log10aT(i) + increment;
end

for i = referenceIndex:nSets-1
    d = d + 1;
    refData = getSet(data, isotherms.Set(i));
    shiftData = getSet(data, isotherms.Set(i+1));
    [increment, diagnostics(d)] = fitShiftPair(refData, shiftData, ...
        p.Results.NumberOfLevels, p.Results.DropFirstPoint);
    diagnostics(d).ReferenceTemperature_C = isotherms.Temperature_C(i);
    diagnostics(d).ShiftedTemperature_C = isotherms.Temperature_C(i+1);
    log10aT(i+1) = log10aT(i) + increment;
end

shiftFactors = isotherms(:, {'Set', 'Temperature_C'});
shiftFactors.log10_aT = log10aT;
shiftFactors.Method = repmat("automatic power-law", nSets, 1);
shiftFactors.ReferenceTemperature_C = repmat( ...
    isotherms.Temperature_C(referenceIndex), nSets, 1);
end

function setData = getSet(data, setValue)
% Isolate one sweep and remove values that cannot be used in a log-frequency
% power-law fit.
rows = data.Set == setValue;
setData = data(rows, :);
valid = isfinite(setData.Frequency_Hz) & setData.Frequency_Hz > 0 & ...
    isfinite(setData.StorageModulus) & setData.StorageModulus > 0;
setData = setData(valid, :);
if height(setData) < 4
    error('dmaEstimateShiftFactors:TooFewPoints', ...
        'Every isotherm needs at least four valid frequency/modulus points.');
end
end

function [log10aT, diagnostic] = fitShiftPair(refData, shiftData, numberOfLevels, dropFirst)
xRef = refData.Frequency_Hz(:);
yRef = refData.StorageModulus(:);
xShift = shiftData.Frequency_Hz(:);
yShift = shiftData.StorageModulus(:);

[pRef, seRef] = fitPowerLaw(xRef, yRef);
[pShift, seShift] = fitPowerLaw(xShift, yShift);
droppedRef = false;
droppedShift = false;

% The first frequency point is often less reliable after a temperature step.
% Remove it only when the fitted exponent becomes more precise.
if dropFirst && numel(xRef) >= 5
    [pCandidate, seCandidate] = fitPowerLaw(xRef(2:end), yRef(2:end));
    if seCandidate(2) < seRef(2)
        pRef = pCandidate;
        seRef = seCandidate;
        xRef = xRef(2:end);
        yRef = yRef(2:end);
        droppedRef = true;
    end
end
if dropFirst && numel(xShift) >= 5
    [pCandidate, seCandidate] = fitPowerLaw(xShift(2:end), yShift(2:end));
    if seCandidate(2) < seShift(2)
        pShift = pCandidate;
        seShift = seCandidate;
        xShift = xShift(2:end);
        yShift = yShift(2:end);
        droppedShift = true;
    end
end

yRefFit = powerLaw(pRef, xRef);
yShiftFit = powerLaw(pShift, xShift);

% Label the curves by vertical position so a positive frequency ratio always
% has the same meaning. Direction restores the original reference-to-shifted
% sign convention afterward.
if max(yRefFit) > max(yShiftFit) && min(yRefFit) > min(yShiftFit)
    topY = yRefFit; topP = pRef;
    bottomY = yShiftFit; bottomP = pShift;
    direction = 1;
elseif max(yRefFit) < max(yShiftFit) && min(yRefFit) < min(yShiftFit)
    topY = yShiftFit; topP = pShift;
    bottomY = yRefFit; bottomP = pRef;
    direction = -1;
elseif max(yRefFit) > max(yShiftFit)
    topY = yShiftFit; topP = pShift;
    bottomY = yRefFit; bottomP = pRef;
    direction = -1;
else
    topY = yRefFit; topP = pRef;
    bottomY = yShiftFit; bottomP = pShift;
    direction = 1;
end

hasOverlap = min(topY) < max(bottomY);
if hasOverlap
    % Compare the fits only over their shared modulus interval.
    yMin = min(topY);
    yMax = max(bottomY);
else
    % There is no shared interval, so bridge the gap by extrapolating both
    % fitted power laws. The diagnostic marks this less reliable case.
    yMin = max(bottomY);
    yMax = min(topY);
end

levels = linspace(yMin, yMax, numberOfLevels).';
xTop = inversePowerLaw(levels, topP);
xBottom = inversePowerLaw(levels, bottomP);
ratio = xTop ./ xBottom;
valid = isfinite(ratio) & ratio > 0;
if nnz(valid) < 3
    error('dmaEstimateShiftFactors:PowerLawIntersection', ...
        'The fitted neighboring isotherms did not produce valid shift levels.');
end
log10aT = direction * mean(log10(ratio(valid)));

diagnostic = emptyDiagnostic();
diagnostic.IncrementalLog10aT = log10aT;
diagnostic.Extrapolated = ~hasOverlap;
diagnostic.DroppedFirstReference = droppedRef;
diagnostic.DroppedFirstShifted = droppedShift;
diagnostic.ReferencePowerLaw = pRef;
diagnostic.ShiftedPowerLaw = pShift;
diagnostic.ReferenceParameterSE = seRef;
diagnostic.ShiftedParameterSE = seShift;
diagnostic.ModulusLevels = levels;
diagnostic.ReferenceFrequencyAtLevel_Hz = ...
    inversePowerLaw(levels, pRef);
diagnostic.ShiftedFrequencyAtLevel_Hz = ...
    inversePowerLaw(levels, pShift);
end

function [parameters, standardError] = fitPowerLaw(x, y)
% A log-log straight-line fit gives a stable starting point for a*x^b + e.
% The offset starts at zero and is bounded to the scale of the measured data.
logFit = polyfit(log(x), log(y), 1);
p0 = [exp(logFit(2)), logFit(1), 0];
lower = [-Inf, -Inf, -max(y)];
upper = [ Inf,  Inf,  max(y)];

if exist('lsqcurvefit', 'file') == 2
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    [parameters, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit( ...
        @(p,xdata) powerLaw(p,xdata), p0, x, y, lower, upper, options);
else
    % Base MATLAB fallback: discourage bound violations with a large penalty.
    scale = max(std(y), eps);
    objective = @(p) sum(((powerLaw(p,x) - y) / scale).^2) + ...
        boundPenalty(p, lower, upper);
    options = optimset('Display', 'off', 'MaxFunEvals', 5000, ...
        'MaxIter', 2000, 'TolX', 1e-10, 'TolFun', 1e-10);
    parameters = fminsearch(objective, p0, options);
    parameters = min(max(parameters, lower), upper);
    residual = powerLaw(parameters, x) - y;
    jacobian = numericalJacobian(@(p) powerLaw(p,x), parameters);
end

dof = max(numel(y) - numel(parameters), 1);

% pinv uses SVD internally, so force a full normal matrix even when the
% optimizer happens to return a sparse Jacobian.
normalMatrix = full(jacobian.' * jacobian);
covariance = pinv(normalMatrix) * sum(residual.^2) / dof;
standardError = sqrt(abs(diag(covariance))).';
end

function value = boundPenalty(p, lower, upper)
lowViolation = max(lower - p, 0);
highViolation = max(p - upper, 0);
value = 1e12 * sum(lowViolation(isfinite(lower)).^2) + ...
    1e12 * sum(highViolation(isfinite(upper)).^2);
end

function y = powerLaw(p, x)
y = p(1) .* x.^p(2) + p(3);
end

function x = inversePowerLaw(y, p)
% sign/abs keeps real arithmetic for small extrapolation excursions caused by
% the fitted offset. Invalid inversions are removed by the caller.
base = (y - p(3)) ./ p(1);
x = sign(base) .* abs(base).^(1 ./ p(2));
end

function jacobian = numericalJacobian(fun, p)
% Used only by the fminsearch fallback to estimate parameter uncertainty.
base = fun(p);
jacobian = zeros(numel(base), numel(p));
for j = 1:numel(p)
    step = sqrt(eps) * max(abs(p(j)), 1);
    trial = p;
    trial(j) = trial(j) + step;
    jacobian(:,j) = (fun(trial) - base) / step;
end
end

function diagnostic = emptyDiagnostic()
% A fixed template makes the pairwise diagnostics straightforward to append.
diagnostic = struct( ...
    'ReferenceTemperature_C', NaN, ...
    'ShiftedTemperature_C', NaN, ...
    'IncrementalLog10aT', NaN, ...
    'Extrapolated', false, ...
    'DroppedFirstReference', false, ...
    'DroppedFirstShifted', false, ...
    'ReferencePowerLaw', nan(1,3), ...
    'ShiftedPowerLaw', nan(1,3), ...
    'ReferenceParameterSE', nan(1,3), ...
    'ShiftedParameterSE', nan(1,3), ...
    'ModulusLevels', [], ...
    'ReferenceFrequencyAtLevel_Hz', [], ...
    'ShiftedFrequencyAtLevel_Hz', []);
end
