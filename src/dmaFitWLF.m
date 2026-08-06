function fitResult = dmaFitWLF(shiftFactors, referenceTemperature, varargin)
%DMAFITWLF Fit the WLF relation to temperature shift factors.
%
% fitResult = dmaFitWLF(shiftFactors, RefT)
% fitResult = dmaFitWLF(..., 'TemperatureRange', [-20 80])
%
% SHIFTFACTORS must contain Temperature_C and log10_aT. TemperatureRange
% limits the fit without changing the input table; individual temperatures
% can also be removed with ExcludeTemperatures.
%
% The bounded fit follows the PyVisco implementation: 0 <= C1,C2 <= 5000.
% lsqcurvefit is used when Optimization Toolbox is available; otherwise a
% bounded fminsearch parameterization is used. The returned structure includes
% the coefficients, residual statistics, uncertainty estimates, and the data
% actually included in the fit.
%
% See also DMAWLF, DMAPREDICTTEMPERATURERESPONSE.

    p           = inputParser;
    p.addRequired('shiftFactors', @istable);
    p.addRequired('referenceTemperature', ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x));
    p.addParameter('TemperatureRange', [-Inf Inf], ...
        @(x) isnumeric(x) && numel(x) == 2 && x(1) <= x(2));
    p.addParameter('ExcludeTemperatures', [], @isnumeric);
    p.parse(shiftFactors, referenceTemperature, varargin{:});
    
    required    = {'Temperature_C', 'log10_aT'};
    if ~all(ismember(required, shiftFactors.Properties.VariableNames))
        error('dmaFitWLF:MissingVariables', ...
            'shiftFactors must contain Temperature_C and log10_aT.');
    end
    
    T       = shiftFactors.Temperature_C(:);
    y       = shiftFactors.log10_aT(:);
    range   = p.Results.TemperatureRange;
    use     = isfinite(T) & isfinite(y) & T >= range(1) & T <= range(2);
    
    % Compare with a small floating-point tolerance so a value typed as 25 is
    % treated the same as an imported value numerically equal to 25.
    for value = p.Results.ExcludeTemperatures(:).'
        use     = use & abs(T - value) > 100*eps(max(abs(value),1));
    end

    T       = T(use);
    y       = y(use);

    if numel(T) < 3
        error('dmaFitWLF:TooFewPoints', 'At least three shift factors are required.');
    end
    
    model   = @(parameters, temperature) dmaWLF(temperature, referenceTemperature, parameters(1), parameters(2));
    lower   = [0 0];
    upper   = [5000 5000];
    initial = [1000 4999];
    
    if exist('lsqcurvefit', 'file') == 2
        options     = optimoptions('lsqcurvefit', 'Display', 'off');
        [parameters, sse, residual, exitFlag, output, ~, jacobian] = ...
            lsqcurvefit(model, initial, T, y, lower, upper, options);
        solver      = 'lsqcurvefit';

    else
        % Map unconstrained fminsearch variables through a logistic curve. This
        % keeps both WLF parameters inside the same bounds used by lsqcurvefit.
        decode      = @(u) upper ./ (1 + exp(-u));
        encode      = @(v) log(v ./ (upper - v));
        objective   = @(u) sum((model(decode(u), T) - y).^2);
        options     = optimset('Display', 'off', 'MaxFunEvals', 5000, ...
                               'MaxIter', 2000, 'TolX', 1e-10, 'TolFun', 1e-12);
        [u, sse, exitFlag, output] = fminsearch(objective, encode(initial), options);
        parameters  = decode(u);
        residual    = model(parameters, T) - y;
        jacobian    = numericalJacobian(@(v) model(v,T), parameters);
        solver      = 'fminsearch';
    end
    
    prediction = model(parameters, T);
    sst     = sum((y - mean(y)).^2);
    if sst > 0
        rSquared    = 1 - sum(residual.^2) / sst;
    else
        rSquared    = NaN;
    end
    rmse            = sqrt(mean(residual.^2));
    
    dof             = max(numel(y) - 2, 1);
    
    % The covariance estimate is the usual local linear approximation. Convert
    % the normal matrix to full storage before pinv because MATLAB's SVD does not
    % accept sparse matrices.
    normalMatrix        = full(jacobian.' * jacobian);
    covariance          = pinv(normalMatrix) * sum(residual.^2) / dof;
    standardError       = sqrt(abs(diag(covariance)));
    
    if exist('tinv', 'file') == 2
        multiplier      = tinv(0.975, dof);
    else
        multiplier      = 1.96;
    end
    
    confidenceInterval95    = [parameters(:) - multiplier*standardError, ...
                               parameters(:) + multiplier*standardError];
    
    fitResult                           = struct;
    fitResult.ReferenceTemperature_C    = referenceTemperature;
    fitResult.C1                        = parameters(1);
    fitResult.C2_C                      = parameters(2);
    fitResult.RMSE_log10aT              = rmse;
    fitResult.RSquared                  = rSquared;
    fitResult.SSE                       = sse;
    fitResult.ExitFlag                  = exitFlag;
    fitResult.Solver                    = solver;
    fitResult.Output                    = output;
    fitResult.Covariance                = covariance;
    fitResult.StandardError             = standardError;
    fitResult.ConfidenceInterval95      = confidenceInterval95;
    fitResult.Data = table(T, y, prediction, residual, ...
                    'VariableNames', {'Temperature_C', 'MeasuredLog10aT', ...
                    'FittedLog10aT', 'Residual'});

end

function jacobian = numericalJacobian(fun, parameters)
    % Forward differences are used only when lsqcurvefit is unavailable.
    base        = fun(parameters);
    jacobian    = zeros(numel(base), numel(parameters));
    for j = 1:numel(parameters)
        step            = sqrt(eps) * max(abs(parameters(j)), 1);
        trial           = parameters;
        trial(j)        = trial(j) + step;
        jacobian(:,j)   = (fun(trial) - base) / step;
    end
end
