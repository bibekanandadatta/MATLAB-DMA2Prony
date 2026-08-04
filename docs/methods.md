# Methods and numerical conventions

This document describes the calculations implemented by the MATLAB functions.
The [API reference](api.md) lists their calling syntax and outputs.

## 1. Input data

The raw frequency-sweep workflow expects one row per measurement with:

- temperature in degrees Celsius;
- cyclic frequency in hertz or angular frequency in radians per second;
- storage modulus;
- loss modulus; and
- preferably a set, sweep, segment, or step identifier.

Tan delta and oscillation strain/amplitude are useful for plotting and quality
review but are not independent inputs to the Prony fit. If tan delta is absent,
the importer calculates

$$
\tan\delta = \frac{E''}{E'}.
$$

Modulus values are not rescaled. A consistent unit must be used throughout the
input. The resulting $E_0$, $E_\infty$, and $E_i$ retain that unit, while
$g_i$ is dimensionless.

An existing master curve can enter the workflow directly. It must contain
`ShiftedFrequency_Hz`, `StorageModulus`, and `LossModulus`. Shift-factor and WLF
fitting require the original temperature sweeps or a separate table of measured
shift factors.

## 2. Isotherm grouping

An exported `Set` identifier is retained when it is present. Otherwise,
temperatures are grouped using a user-selected tolerance. The representative
temperature of a set is the mean of its measurement temperatures.

The reference temperature requested for automatic shifting must be close to an
actually measured isotherm. The nearest measured set is used, and its shift
factor is defined as

$$
a_T = 1, \qquad \log_{10}(a_T)=0.
$$

## 3. Automatic horizontal shifting

The automatic shift procedure follows the neighboring-isotherm power-law
approach used by PyVisco. Storage modulus is used for alignment; loss modulus
and tan delta receive the same shift and are never aligned independently.

For each pair of neighboring isotherms, storage modulus is fitted to

$$
E'(f) = a f^b + e.
$$

The fit is performed with all measurements and, optionally, after dropping the
first measurement. The first point is discarded when its removal reduces the
estimated standard error of exponent $b$. The horizontal separation of the
chosen power laws is then evaluated at common modulus levels and averaged in
log-frequency space.

Incremental shifts are accumulated outward from the reference set. Therefore,
an uncertain pairwise shift affects all subsequent temperatures in that
direction.

### Solver sensitivity

MATLAB `lsqcurvefit` and SciPy `curve_fit` can calculate slightly different
parameter covariance estimates. Consequently, the first-point discard decision
is not guaranteed to match PyVisco exactly even when the fitted power-law curves
are nearly identical. On the included dataset, forcing identical discard
decisions makes the pairwise shift increments agree closely, while independent
discard decisions cause larger accumulated differences on the cold side.

The returned diagnostics expose:

- whether the two fitted modulus ranges overlap;
- whether either first point was dropped;
- fitted power-law parameters and standard errors; and
- the modulus levels and inferred frequencies used for shifting.

These diagnostics should be reviewed rather than treating automatic shifts as
experimental ground truth. Automatic first-point removal can be disabled with
`'DropFirstPoint', false`, and the returned shift-factor table can be edited
before constructing the master curve.

## 4. Master-curve convention

This package uses

$$
f_{\mathrm{reduced}} = f_{\mathrm{measured}}a_T
                      = f_{\mathrm{measured}}10^{\log_{10}(a_T)}.
$$

The master curve is sorted by reduced frequency. Its extreme frequency values
are shifted coordinates, not frequencies applied directly by the instrument.

Changing reference temperature translates the master curve horizontally. Under
ideal horizontal time-temperature superposition, the modulus fractions remain
unchanged and the relaxation times scale by the shift factor between reference
temperatures. Prony times and WLF coefficients must always use the same
reference temperature.

Choose a measured reference isotherm with reliable data, preferably near both
the center of the characterization range and the temperature range of practical
interest. Avoid an extreme reference unless that temperature is specifically
required.

## 5. Master-curve filtering

The optional filter is a centered moving median, corresponding to the approach
used by PyVisco:

$$
\widetilde{E'}_j = \operatorname{median}\!\left(E'\text{ in the local window}\right),
$$

with the same operation applied to $E''$. Filtered tan delta is recomputed as

$$
\widetilde{\tan\delta}=\frac{\widetilde{E''}}{\widetilde{E'}}.
$$

The measured columns are preserved. A window of one performs no effective
filtering.

Filtering reduces isolated spikes and local jaggedness. It does not repair
incorrect shifts, systematic endpoint artifacts, or failure of
time-temperature superposition. Smooth publication curves are obtained from
the fitted generalized-Maxwell model, not by heavily smoothing the measured
curve.

Filtering is deliberately performed only by `dmaFilterMasterCurve`.
`dmaFitProny` does not apply a hidden second filter.

## 6. WLF relation

The implemented Williams-Landel-Ferry relation is

$$
\log_{10}(a_T) =
-\frac{C_1(T-T_{\mathrm{ref}})}{C_2+(T-T_{\mathrm{ref}})}.
$$

The fit constrains $C_1$ and $C_2$ to the interval $[0,5000]$, following
the bounds used by PyVisco. The output includes residuals, RMSE, $R^2$,
covariance, standard errors, approximate 95% confidence intervals, and solver
status.

WLF should only be used over a temperature range for which the shift relation
and time-temperature superposition are physically defensible. Temperatures can
be excluded from the fit without deleting the original shift factors.

## 7. Generalized-Maxwell/Prony model

For sinusoidal tension/compression,

$$
\omega = 2\pi f.
$$

Angular frequency describes temporal phase rate and does not imply torsional
loading. Relaxation times and their collocation frequencies are related by

$$
\tau_i = \frac{1}{\omega_i}=\frac{1}{2\pi f_i}.
$$

The generalized-Maxwell storage and loss moduli are

$$
E'(\omega)=E_\infty+
\sum_i E_i\frac{(\omega\tau_i)^2}{1+(\omega\tau_i)^2},
$$

$$
E''(\omega)=
\sum_i E_i\frac{\omega\tau_i}{1+(\omega\tau_i)^2}.
$$

The normalized and dimensional branch coefficients satisfy

$$
E_i=E_0g_i,
\qquad
E_\infty=E_0\left(1-\sum_i g_i\right).
$$

The time-domain relaxation modulus is

$$
E(t)=E_\infty+\sum_i E_i\exp(-t/\tau_i).
$$

### Fixed-time fitting

For `exact`, `round`, and `manual` modes, relaxation times are placed
logarithmically in the selected window. Storage and loss modulus are
interpolated at the collocation frequencies, and nonnegative generalized
collocation is solved with `lsqnonneg`. Coefficients are normalized if necessary
so that $\sum_i g_i\le 1$.

The number of terms is always supplied by the user. Term-count comparison is
provided as a diagnostic and does not silently select a model.

### Relaxation-time window modes

- `exact`: use the exact $1/(2\pi f)$ endpoints corresponding to the fitting
  data. This is the package default.
- `round`: use the interior power-of-ten endpoints employed by PyVisco.
- `manual`: logarithmically distribute terms over a user-supplied range in
  seconds.
- `min`: start with interior terms and jointly optimize $g_i$ and $\tau_i$
  with `fmincon`.

The frequency data admitted to the fitting objective can be restricted
independently of the relaxation-time window. The code warns when a prescribed
relaxation window does not cover the selected fitting-frequency range under
$f_i=1/(2\pi\tau_i)$.

## 8. Dense constitutive reconstruction

The fitted equations are evaluated on a user-controlled logarithmic grid to
produce smooth storage modulus, loss modulus, tan delta, and relaxation modulus
curves. This calculation is analytic at each requested point. It is separate
from filtering and can also evaluate an explicit list of frequencies.

For identical master curves, filtering settings, term counts, and exact
relaxation windows, the included MATLAB fixed-time Prony implementation has
been compared with PyVisco and agrees to numerical precision.

## 9. Temperature-dependent prediction

At physical frequency $f$ and temperature $T$, the reduced frequency is

$$
f_r=f\,10^{\log_{10}(a_T(T))}.
$$

The generalized-Maxwell model fitted at the reference temperature is evaluated
at $f_r$. The returned table contains physical frequency, reduced frequency,
shift factor, storage modulus, loss modulus, tan delta, and coverage flags.

By default, predictions outside either the temperatures used in the WLF fit or
the frequencies used in the Prony fit are returned as `NaN`. Extrapolation must
be enabled explicitly and should be treated cautiously.

## 10. Scope and limitations

- The package assumes horizontal time-temperature superposition.
- It does not decide whether questionable endpoint measurements are physically
  valid.
- It does not independently shift storage, loss, and tan delta.
- It does not infer WLF parameters from a master curve alone.
- It does not silently select the number of Prony terms.
- Modulus conversion and tensile/shear interpretation remain the user's
  responsibility.
