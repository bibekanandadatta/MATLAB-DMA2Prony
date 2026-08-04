function filteredCurve = dmaFilterMasterCurve(masterCurve, varargin)
%DMAFILTERMASTERCURVE Apply PyVisco-style moving-median filtering.
%
% filteredCurve = dmaFilterMasterCurve(masterCurve)
% filteredCurve = dmaFilterMasterCurve(masterCurve, 'MedianWindow', 5)
%
% MASTERCURVE must contain ShiftedFrequency_Hz, StorageModulus, and
% LossModulus. MedianWindow is the number of neighboring master-curve points
% used by the centered median; its default is 5.
%
% The measured columns are preserved. The function adds:
%   StorageModulusFiltered
%   LossModulusFiltered
%   TanDeltaFiltered
%
% TanDeltaFiltered is calculated from the two filtered moduli. It is not
% filtered independently.
%
% See also DMAFITPRONY.

p = inputParser;
p.addRequired('masterCurve', @istable);
p.addParameter('MedianWindow', 5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == round(x));
p.parse(masterCurve, varargin{:});

required = {'ShiftedFrequency_Hz', 'StorageModulus', 'LossModulus'};
missing = required(~ismember(required, masterCurve.Properties.VariableNames));
if ~isempty(missing)
    error('dmaFilterMasterCurve:MissingVariables', ...
        'Missing required variables: %s', strjoin(missing, ', '));
end

filteredCurve = sortrows(masterCurve, 'ShiftedFrequency_Hz');
validFrequency = isfinite(filteredCurve.ShiftedFrequency_Hz) & ...
    filteredCurve.ShiftedFrequency_Hz > 0;

% Filter storage and loss separately, then form their ratio. Filtering tan
% delta on its own would make the three output columns inconsistent.
storageFiltered = nan(height(filteredCurve), 1);
lossFiltered = nan(height(filteredCurve), 1);
storageFiltered(validFrequency) = movingMedianFinite( ...
    filteredCurve.StorageModulus(validFrequency), p.Results.MedianWindow);
lossFiltered(validFrequency) = movingMedianFinite( ...
    filteredCurve.LossModulus(validFrequency), p.Results.MedianWindow);

tanDeltaFiltered = lossFiltered ./ storageFiltered;
tanDeltaFiltered(~isfinite(tanDeltaFiltered) | storageFiltered <= 0) = NaN;

filteredCurve.StorageModulusFiltered = storageFiltered;
filteredCurve.LossModulusFiltered = lossFiltered;
filteredCurve.TanDeltaFiltered = tanDeltaFiltered;
end

function filtered = movingMedianFinite(values, window)
% Invalid measurements keep their original positions and are left as NaN.
% They are excluded from the rolling sequence so they do not contaminate
% otherwise valid neighboring points.
values = values(:);
filtered = nan(size(values));
valid = isfinite(values);
if any(valid)
    filtered(valid) = movmedian(values(valid), window, 'Endpoints', 'shrink');
end
end
