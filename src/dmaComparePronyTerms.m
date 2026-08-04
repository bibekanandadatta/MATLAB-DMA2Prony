function [comparison, fits] = dmaComparePronyTerms( ...
    masterCurve, termCounts, varargin)
%DMACOMPAREPRONYTERMS Fit several requested Prony-series term counts.
%
% [comparison, fits] = dmaComparePronyTerms(masterCurve, [5 10 15], ...)
%
% Additional name-value arguments are passed directly to dmaFitProny. The
% function reports fit errors but deliberately leaves the engineering choice
% of term count to the user. COMPARISON is a summary table; FITS is a cell
% array containing the complete result structure for each term count.
%
% See also DMAFITPRONY.

if ~istable(masterCurve)
    error('dmaComparePronyTerms:InvalidMasterCurve', ...
        'masterCurve must be a table.');
end
if ~(isnumeric(termCounts) && isvector(termCounts) && ~isempty(termCounts) && ...
        all(isfinite(termCounts)) && all(termCounts >= 1) && ...
        all(termCounts == round(termCounts)))
    error('dmaComparePronyTerms:InvalidTermCounts', ...
        'termCounts must contain positive integers.');
end

termCounts = unique(termCounts(:), 'sorted');
numberOfFits = numel(termCounts);

% Preallocate the summary columns so every fit is stored in the same order as
% the reported term counts.
fits = cell(numberOfFits, 1);
rmseStorage = nan(numberOfFits, 1);
rmseLoss = nan(numberOfFits, 1);
nrmseStorage = nan(numberOfFits, 1);
nrmseLoss = nan(numberOfFits, 1);
residualNorm = nan(numberOfFits, 1);

for i = 1:numberOfFits
    fits{i} = dmaFitProny(masterCurve, termCounts(i), varargin{:});
    rmseStorage(i) = fits{i}.RMSE_StorageModulus;
    rmseLoss(i) = fits{i}.RMSE_LossModulus;
    nrmseStorage(i) = fits{i}.NRMSE_StorageModulus;
    nrmseLoss(i) = fits{i}.NRMSE_LossModulus;
    residualNorm(i) = fits{i}.ResidualNorm;
end

% This single measure is convenient for comparison, but the storage and loss
% errors are kept separately because they need not favor the same model.
combinedNRMSE = sqrt((nrmseStorage.^2 + nrmseLoss.^2)/2);
comparison = table(termCounts, rmseStorage, rmseLoss, nrmseStorage, ...
    nrmseLoss, combinedNRMSE, residualNorm, ...
    'VariableNames', {'NumberOfTerms', 'RMSE_StorageModulus', ...
    'RMSE_LossModulus', 'NRMSE_StorageModulus', ...
    'NRMSE_LossModulus', 'CombinedNRMSE', 'ResidualNorm'});
end
