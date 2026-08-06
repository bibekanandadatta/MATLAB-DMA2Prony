function data = dmaReadTAExcel(filename, varargin)
%DMAREADTAEXCEL Read a TA/TRIOS Excel data table into canonical DMA columns.
%
% data = dmaReadTAExcel(filename)
% data = dmaReadTAExcel(filename, 'Sheet', 'DMA Data')
%
% The function expects one tabular worksheet with one row per measurement.
% It recognizes common TA/TRIOS header spellings for temperature, frequency,
% angular frequency, storage modulus, loss modulus, tan(delta), oscillation
% strain/stress, amplitude, time, and step/set number. Numerical values are
% not rescaled; the units selected during export are retained.
%
% DATA is a table with consistent variable names used by the remaining DMA
% functions. SourceRow and workbook metadata are retained so questionable
% measurements can be traced back to the original worksheet.
%
% See also DMAPLOT, DMAGROUPISOTHERMS.

    p       = inputParser;
    p.addRequired('filename', @(x) ischar(x) || isstring(x));
    p.addParameter('Sheet', 1, @(x) isnumeric(x) || ischar(x) || isstring(x));
    p.parse(filename, varargin{:});
    
    filename = char(p.Results.filename);

    if ~isfile(filename)
        error('dmaReadTAExcel:FileNotFound', 'File not found: %s', filename);
    end
    
    raw = readtable(filename, 'Sheet', p.Results.Sheet,'VariableNamingRule', 'preserve');
    
    if isempty(raw)
        error('dmaReadTAExcel:EmptySheet', 'The selected worksheet is empty.');
    end
    
    names   = string(raw.Properties.VariableNames);
    
    % Keep all header matching in one place. Exact normalized names are preferred;
    % findColumn uses partial matching only when an exact alias is unavailable.
    idx.temperature     = findColumn(names, ...
        ["temperature", "temperaturec", "temp", "tempc"]);
    idx.frequency       = findColumn(names, ...
        ["frequencyhz", "frequency", "freqhz", "freq", "f"], "angular");
    idx.omega           = findColumn(names, ...
        ["angularfrequency", "angularfrequencyrads", "omegarads", "omega"]);
    idx.storage         = findColumn(names, ...
        ["storagemodulus", "estor", "gstor", "storage"]);
    idx.loss            = findColumn(names, ...
        ["lossmodulus", "eloss", "gloss", "loss"]);
    idx.tanDelta        = findColumn(names, ...
        ["tandelta", "tandel", "lossfactor", "dampingfactor"]);
    idx.strain          = findColumn(names, ...
        ["oscillationstrain", "dynamicstrain", "strainpercent", "strain"], "rate");
    idx.stress          = findColumn(names, ...
        ["oscillationstress", "dynamicstress", "stress"], "rate");
    idx.amplitude       = findColumn(names, ...
        ["oscillationamplitude", "dynamicamplitude", "amplitude"]);
    idx.time            = findColumn(names, ...
        ["steptimes", "times", "steptime", "time"]);
    idx.set             = findColumn(names, ...
        ["set", "sweep", "segment", "stepnumber", "step"]);
    
    if isempty(idx.storage)
        error('dmaReadTAExcel:MissingStorageModulus', ...
            ['No storage-modulus column was recognized. Export a column named ', ...
             '"Storage Modulus" (or E''/G'').']);
    end

    sourceRow   = (1:height(raw)).';
    data        = table(sourceRow, 'VariableNames', {'SourceRow'});
    
    % Build a clean table rather than carrying instrument-specific headers through
    % the numerical routines. Missing optional quantities are simply omitted.
    data    = addNumeric(data, raw, idx.temperature, 'Temperature_C');
    data    = addNumeric(data, raw, idx.frequency, 'Frequency_Hz');
    data    = addNumeric(data, raw, idx.omega, 'AngularFrequency_rad_s');
    data    = addNumeric(data, raw, idx.storage, 'StorageModulus');
    data    = addNumeric(data, raw, idx.loss, 'LossModulus');
    data    = addNumeric(data, raw, idx.tanDelta, 'TanDelta');
    data    = addNumeric(data, raw, idx.strain, 'OscillationStrain_pct');
    data    = addNumeric(data, raw, idx.stress, 'OscillationStress');
    data    = addNumeric(data, raw, idx.amplitude, 'Amplitude');
    data    = addNumeric(data, raw, idx.time, 'Time_s');
    
    if ~isempty(idx.set)
    
        setColumn   = raw{:, idx.set};
    
        if isnumeric(setColumn) || islogical(setColumn)
            data.Set = double(setColumn);
        else
            % Text labels are converted to stable integer groups for later joins.
            [data.Set, ~] = findgroups(string(setColumn));
        end
    end
    
    % Remove unit rows, blank rows, and report/footer rows that do not contain
    % a storage-modulus result.
    data    = data(isfinite(data.StorageModulus), :);
    if isempty(data)
        error('dmaReadTAExcel:NoNumericData', ...
            'No numerical storage-modulus rows were found in the selected sheet.');
    end
    
    if ~ismember('Frequency_Hz', data.Properties.VariableNames) && ...
            ismember('AngularFrequency_rad_s', data.Properties.VariableNames)
        % Frequency in hertz is the package's primary frequency coordinate.
        data.Frequency_Hz = data.AngularFrequency_rad_s / (2*pi);
    end
    if ~ismember('AngularFrequency_rad_s', data.Properties.VariableNames) && ...
            ismember('Frequency_Hz', data.Properties.VariableNames)
        data.AngularFrequency_rad_s = 2*pi*data.Frequency_Hz;
    end
    if ~ismember('TanDelta', data.Properties.VariableNames) && ...
            all(ismember({'StorageModulus', 'LossModulus'}, data.Properties.VariableNames))
        % Derive tan(delta) only when it was not exported by the instrument.
        data.TanDelta = data.LossModulus ./ data.StorageModulus;
    end
    
    data.Properties.Description     = sprintf('DMA data imported from %s', filename);
    metadata                        = struct;
    metadata.SourceFile             = filename;
    metadata.SourceSheet            = p.Results.Sheet;
    metadata.SourceVariableNames    = names;
    data.Properties.UserData        = metadata;
end

function data = addNumeric(data, raw, idx, outputName)
    if isempty(idx)
        return
    end
    data.(outputName)           = toNumeric(raw{:, idx}, raw.Properties.VariableNames{idx});
end


function values = toNumeric(column, columnName)
% readtable may return a numeric, duration, or text column depending on the
% workbook. Convert all three cases without guessing or changing units.
    if isnumeric(column) || islogical(column)
        values  = double(column);
    elseif isduration(column)
        values  = seconds(column);
    else
        text    = strip(string(column));
        text    = replace(text, ",", "");
        values  = str2double(text);
    end
    
    values      = values(:);
    
    if all(~isfinite(values))
        error('dmaReadTAExcel:NonNumericColumn', ...
            'Recognized column "%s", but it contains no numerical values.', columnName);
    end
end


function idx = findColumn(names, aliases, excludedText)
    if nargin < 3
        excludedText = "";
    end
    
    normalized      = normalize(names);
    aliases         = normalize(aliases);
    excludedText    = normalize(string(excludedText));
    
    idx             = [];
    
    % Exact matches avoid confusing short aliases such as "f" or "time" with a
    % longer, unrelated column name.
    for i = 1:numel(aliases)
        hit = find(normalized == aliases(i), 1);
        if ~isempty(hit)
            idx = hit;
            return
        end
    end
    
    % Some TRIOS exports append units or specimen notes to the header. A partial
    % match catches those variants after the exact pass has failed.
    for i = 1:numel(aliases)
        hit = find(contains(normalized, aliases(i)), 1);
        if ~isempty(hit)
            if strlength(excludedText) == 0 || ~contains(normalized(hit), excludedText)
                idx = hit;
                return
            end
        end
    end
end

function text = normalize(text)
    text    = lower(string(text));
    text    = regexprep(text, '[^a-z0-9]', '');
end
