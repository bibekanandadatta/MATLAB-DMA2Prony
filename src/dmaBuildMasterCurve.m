function masterCurve = dmaBuildMasterCurve(data, shiftFactors, varargin)
%DMABUILDMASTERCURVE Apply horizontal shift factors to DMA frequency sweeps.
%
% masterCurve = dmaBuildMasterCurve(data, shiftFactors)
% masterCurve = dmaBuildMasterCurve(..., 'TemperatureTolerance', 0.5)
%
% DATA must contain Frequency_Hz and Temperature_C. SHIFTFACTORS must contain
% one log10_aT value for every isothermal Set. The measured modulus columns
% are copied without modification; only the frequency axis is shifted:
%   f_shifted = f * 10^(log10_aT)
%
% To manually correct an automatic result, edit shiftFactors.log10_aT before
% calling this function.
%
% See also DMAESTIMATESHIFTFACTORS, DMAGROUPISOTHERMS.

    p       = inputParser;
    
    p.addRequired('data', @istable);
    p.addRequired('shiftFactors', @istable);
    p.addParameter('TemperatureTolerance', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.parse(data, shiftFactors, varargin{:});
    
    requiredData        = {'Frequency_Hz', 'Temperature_C'};
    requiredShift       = {'Set', 'log10_aT'};
    missingData         = requiredData(~ismember(requiredData, data.Properties.VariableNames));
    missingShift        = requiredShift(~ismember(requiredShift, shiftFactors.Properties.VariableNames));
    if ~isempty(missingData) || ~isempty(missingShift)
        error('dmaBuildMasterCurve:MissingVariables', ...
            'The data or shift-factor table is missing required variables.');
    end
    
    [data, ~]           = dmaGroupIsotherms(data, ...
        'TemperatureTolerance', p.Results.TemperatureTolerance);
    
    % Match on Set rather than table row number. This also allows the user to
    % reorder the shift-factor table after making manual corrections.
    [found, location] = ismember(data.Set, shiftFactors.Set);
    if ~all(found)
        error('dmaBuildMasterCurve:MissingShiftFactor', ...
            'A shift factor is not available for every measurement set.');
    end
    
    masterCurve                                 = data;
    masterCurve.log10_aT                        = shiftFactors.log10_aT(location);
    masterCurve.ShiftedFrequency_Hz             = masterCurve.Frequency_Hz .* 10.^masterCurve.log10_aT;
    masterCurve.ShiftedAngularFrequency_rad_s   = 2*pi*masterCurve.ShiftedFrequency_Hz;
    
    % A monotonic frequency column makes later filtering and interpolation
    % predictable, especially where shifted isotherms overlap.
    masterCurve                                 = sortrows(masterCurve, 'ShiftedFrequency_Hz');

end
