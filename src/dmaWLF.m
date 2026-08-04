function log10aT = dmaWLF(temperature, referenceTemperature, C1, C2)
%DMAWLF Evaluate the Williams-Landel-Ferry shift relation.
%
% log10aT = dmaWLF(T, RefT, C1, C2)
%
% T and RefT are temperatures in the same units, normally degrees Celsius.
% C1 is dimensionless and C2 has temperature units. Array inputs are handled
% element by element, and the result is log10(aT).
%
% See also DMAFITWLF, DMAPREDICTTEMPERATURERESPONSE.

deltaT = temperature - referenceTemperature;
% Writing the equation in terms of deltaT keeps the reference condition
% explicit: deltaT = 0 gives log10(aT) = 0.
log10aT = -C1 .* deltaT ./ (C2 + deltaT);
end
