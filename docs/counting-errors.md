# Counting errors in pBASEX

This derivation uses the notation of the
[original pBASEX note](pbasex.pdf) and the
[half-plane derivation](forward-backward-asymmetry.md). It gives the
counting uncertainty of yields and anisotropies for fixed analysis
settings. It applies to both detector geometries.

For a first reading, [From detected electrons to pBASEX error bars](uncertainty-derivation.md)
starts with a projected electron distribution and follows one measured pixel
through the fit, energy integration and error calculation.

## Calculate counting errors in Julia

Each call to `pbasex_counts` performs **one Abel inversion**. It returns the
native-grid spectrum and betas with their counting errors. In the same call,
you can request integrated results for any number of energy intervals. All
intervals use the same fitted coefficients.

The examples below use an original count histogram, its measured one-based
center, and a basis file for the chosen geometry. Define the shared settings
once, then choose a call below. With `alpha=1.0`, all energies are in px².
For a calibration in eV, supply that calibration and specify all limits in eV.

```julia
using PBASEX

settings = (; basis=basis_path, geometry=:half_plane, center=(column, row),
              alpha=1.0, make_images=false)
```

Use `:half_plane` for a cylindrically symmetric image with distinct forward
and backward hemispheres. Use `:quadrant` with a quadrant basis when both
reflection symmetries hold. `make_images=false` skips constructing the images.

Supply `rotation=90` in `settings` if the original right-hand side of the
detector should point up along the positive fit axis. Both functions use
one-based `(column, row)` bins before rotation. With no center override,
`pbasex_counts` uses the histogram midpoint (the higher index at an even
dimension). `pbasex_shots` uses the center stored in the checkpoint.
`pbasex_counts` can build a basis at a new path using ordinary `pbasex`
basis options. `pbasex_shots` reads an existing basis. For a comparison,
use the same file, detector coverage, quadrant selection and `reg`.

### Spectrum and betas on the native grid

```julia
analysis = pbasex_counts(counts; settings...)
energy, spectrum, spectrum_se = analysis.result.E, analysis.result.IE, analysis.spectrum_se
orders = analysis.legendre_orders[2:end]
betas, beta_se = analysis.result.betas, analysis.beta_se
```

Here `energy[i]` is a sample on the native grid $E_i=\alpha r_i^2$.
Each `spectrum_se[i]` is a pointwise standard error in the units of `IE`.
The dimensionless `beta_se[i,j]` is the standard error of `betas[i,j]` at
`energy[i]` for Legendre order `orders[j]`. It includes the covariance
between the angular component and the isotropic spectrum at that energy.
It is `NaN` where `IE[i] <= spectrum_se[i]`, including the origin.

### Integrated yield and betas for selected peaks

Pass one pair of limits per interval. Any interval within the sampled energy
range is valid if its lower limit is smaller than its upper limit. The limits
can fall between samples. Choose limits appropriate to your data.

```julia
windows = [(32_000, 34_000), (34_500, 36_500)]
analysis = pbasex_counts(counts; settings..., windows)
orders = analysis.legendre_orders[2:end]
peak = analysis.windows[1]       # First interval, in the order supplied
peak.window, peak.yield, peak.yield_se
peak.betas, peak.beta_se       # One value and standard error per order
```

Both `pbasex_counts` and `pbasex_shots` accept `windows` and return these fields in input order. Existing counts analyses can still use `peak_windows` and `analysis.peaks`.

For raw electron counts, `peak.yield` is the fitted number of detected
electrons in the interval. A beta describes their combined anisotropy.
A zero yield gives undefined betas. For a negative fitted yield the
counts fit returns the algebraic ratio, whose denominator has no physical
interpretation as an electron yield. The shot fit returns `NaN` for that
beta. When the fitted yield in the counts fit is at or below
one estimated standard error, `peak.beta_se` is `NaN`. Even above that
threshold, symmetric beta errors are approximate for weak signals.

### Yield-preserving spectrum and betas on a uniform energy grid

Give each bin its own window in one call. The example extracts the binned
spectrum and $\beta_2$ for a basis that includes order 2.

```julia
edges = collect(30_000:500:40_000)
windows = collect(zip(edges[1:end-1], edges[2:end]))
analysis = pbasex_counts(counts; settings..., windows)

energy_mid = (edges[1:end-1] .+ edges[2:end]) ./ 2
spectrum = [bin.spectrum for bin in analysis.windows]
spectrum_se = [bin.spectrum_se for bin in analysis.windows]

orders = analysis.legendre_orders[2:end]
j = only(findall(==(2), orders))
beta2 = [bin.betas[j] for bin in analysis.windows]
beta2_se = [bin.beta_se[j] for bin in analysis.windows]
```

For a bin $[a,b]$, the code integrates each piecewise-linear angular
component $I_\ell(E)$ on the native grid, including the partial segments
at the bin edges. Writing $B_\ell=\int_a^b I_\ell(E)\,dE$, the fitted count
is $B_0$, the bin-averaged spectrum is $B_0/(b-a)$, and the bin anisotropy
is $B_\ell/B_0$. Only the spectrum and its error are divided by bin width.
Use `energy_mid` as the plotting coordinate for these bin averages.

Adjacent bin yields sum to the integral of the fitted spectrum over the
combined interval. Their errors can be correlated, especially for bins
finer than the native spacing. To obtain the counting error for a combined
interval, include that interval as another window in the same call. The
binned outputs live in `analysis.windows`; `analysis.result.E`,
`analysis.result.IE`, `analysis.spectrum_se`, and `analysis.beta_se` always
refer to the native samples.

### Count statistics

This entry point assumes independent Poisson counts in the supplied pixels,
a fixed center, and fixed analysis settings. An integer histogram or its
integer-valued floating-point representation is suitable. Pixel permutations
such as a 90° rotation preserve the model. Interpolated or tapered images,
background-subtracted data, and fixed-total acquisition require the
corresponding propagation described below. The returned uncertainties
cover counting statistics only.

## Shot-resolved uncertainties from a corrected checkpoint

For a worked three-shot example and a derivation of the bootstrap estimator,
see [How shot-resolved pBASEX works](shot-uncertainties.md).

Use `pbasex_shots` when each corrected electron retains its shot ID. The
estimate pools the angular contributions from all shots. The bootstrap
resamples whole shots and takes each beta ratio after summation. Its spread
includes counting fluctuations and the variation between observed shots.

```julia
using PBASEX, Random
analysis = pbasex_shots(checkpoint_path; basis=basis_path,
    windows=[(Ea, Eb), (Ec, Ed)], alpha=1.0, rotation=90,
    replicates=2000, confidence=0.95, rng=Xoshiro(42))
w = analysis.windows[1]
w.yield, w.yield_se, w.yield_ci
w.spectrum, w.spectrum_se           # yield and standard error per unit energy
analysis.legendre_orders[2:end], w.betas, w.beta_se, w.beta_ci
```

Supply any intervals inside the basis energy range. Adjacent equal-width
intervals give a uniform spectrum. The basis file sets the geometry and
Legendre orders. Energy calibration is per squared detector pixel. The
reader accounts for the histogram resolution and uses the stored center
unless `center=(column,row)` overrides it in one-based **unrotated** bins.
Both estimators accept signed rotations in multiples of 90 degrees (default
zero). Use `rotation=90` for a counterclockwise 90° rotation. A half-plane
basis must match the detector coverage after centering and rotation.

The reader accepts GUI schema 1 checkpoints with these HDF5 entries.

| Location | Contents |
| --- | --- |
| Root attribute `schema_version` | Integer 1 |
| `shots/id` | One-dimensional integer roster in acquisition order, including zero-hit shots. At least two IDs, each unique. |
| `events/shot_id`, `events/x_corrected`, `events/y_corrected` | Three one-dimensional columns of the same length. Each row describes one event and its integer shot ID from the roster. Events may be interleaved. |
| `geometry/x_edges`, `geometry/y_edges` | Increasing bin edges in detector pixels with spacing `1/resolution_factor`. |
| Attributes on `geometry` | `center_xy` in detector coordinates for the default center, positive `resolution_factor`, `x_direction="right"`, `y_direction="down"`, `coordinate_units="detector_pixel"`. |

Bins include their lower edge and the final upper edge. The recorded event
coordinates already include the distortion correction. Fixed `shot_weights`
can be supplied in roster order. Geometry, calibration and any supplied
normalization factors are held fixed through all replicates.

The default bootstrap samples individual exchangeable shots. For serial
dependence within stationary runs, pass `block_length` and optionally
`run_ids` in roster order. Consecutive equal IDs define a run. Each run
retains its shot count. Choose a block length shorter than the run and
long enough to capture the measured correlation. The intervals are
pointwise percentile intervals. A nonpositive pooled yield or a lower
yield confidence limit at or below zero gives undefined beta error bars.
`valid_replicates` counts draws with positive yield.

All windows share draws. With `keep_replicates=true`, each window provides
`replicates.yield[draw]` and `replicates.betas[order_index,draw]`. These
support arbitrary paired contrasts. To stream draws to a file instead,
pass a callback `replicate_sink=(window_index, draws) -> ...`.

The defaults process one window at a time, up to 65,536 events per read,
and approximately 64 MiB of singular-vector columns per read (at least
one column). For a batch of $m$ window/order combinations, the response
and shot tables together use about $8m(P+S)$ bytes for $P$ folded pixels
and $S$ shots, in addition to read buffers, bootstrap samples and runtime
overhead. Increase `window_batch` to trade memory for fewer file passes.

## 1. From detector counts to basis coefficients

Let $n_i$ be the number of centroided electron events in original
detector pixel $i$, with mean $\mu_i$. Suppose the total number of detected
events in an exposure is Poisson and each event independently lands in one
pixel with fixed probabilities. The pixel counts are then independent Poisson
variables, so

$$
\mathbb E[\mathbf n]=\boldsymbol\mu,\qquad
\operatorname{Cov}(\mathbf n)=D_{\mathbf n}
=\operatorname{diag}(\boldsymbol\mu).
\tag{1}
$$

With the center, geometry and preprocessing fixed, write the fitted
projection as $\mathbf b=F\mathbf n$. The linear map $F$ includes any
fixed rotation, interpolation, taper, detector selection and folding.
For a fold applied directly to raw counts, the following entries of
$\mathbf b$ illustrate the two conventions. Here $i$ and $j$ label
distinct pixels and $c$ labels the center pixel.

| Fold and position | Fitted pixel value | Variance under (1) |
| --- | --- | --- |
| Half-plane, two mirror pixels | $(n_i+n_j)/2$ | $(\mu_i+\mu_j)/4$ |
| Half-plane, on the symmetry axis | $n_i$ | $\mu_i$ |
| Four quadrants selected, at the center | $n_c$ | $\mu_c$ |
| Four quadrants selected, on either axis away from the center | $(n_i+n_j)/2$ | $(\mu_i+\mu_j)/4$ |

The quadrant fit divides the selected sum by the number of selected quadrants.
At the center all four summands are the same observation, so averaging leaves
its variance intact. Away from both axes four distinct observations contribute,
giving variance $(\mu_1+\mu_2+\mu_3+\mu_4)/16$ when all four are selected.
If interpolation precedes folding, propagate the composite map from raw
pixels. Its output positions can share the same input noise.
The fitted projection therefore has covariance $F D_{\mathbf n}F^\mathsf T$.
For incomplete half-plane coverage, take $F$ and $G$ on the same
retained detector positions.

Write $G=U\Sigma V^\mathsf T$ for the SVD of the forward operator on
the fitted detector grid. Let $H=V\operatorname{diag}(h_j)U^\mathsf T$
be the inverse operator used by the fit. If $s_j$ is a singular value,
then $h_j=s_j/(s_j^2+\lambda)$ for regularization $\lambda>0$. For
$\lambda=0$, $h_j=1/s_j$ above the numerical cutoff and zero otherwise.
The fitted coefficients and their covariance are

$$
\begin{aligned}
\widehat{\mathbf c}&=HF\mathbf n,\\
\operatorname{Cov}(\widehat{\mathbf c})
&=HF D_{\mathbf n}F^\mathsf T H^\mathsf T.
\end{aligned}
\tag{2}
$$

Equation (2) is exact for the fixed linear estimator under (1). It
describes its variation about its expectation. A small variance alone
does not establish that the radial basis or regularization gives an
unbiased fit.

## 2. Energy spectra and peak yields

With $E=\alpha r^2$ and $r=\sqrt{E/\alpha}$, define the fitted component
of order $\ell$ per unit energy as

$$
I_\ell(E)=\frac{r}{2\alpha}\sum_k \widehat c_{k\ell} f_k(r).
\tag{3}
$$

The angle-integrated spectrum is $I(E)=I_0(E)$. For an energy window
$[E_a,E_b]$, define the integrated angular components

$$
B_\ell=\int_{E_a}^{E_b} I_\ell(E)\,dE
=\sum_k\widehat c_{k\ell}
\int_{\sqrt{E_a/\alpha}}^{\sqrt{E_b/\alpha}} r^2f_k(r)\,dr.
\tag{4}
$$

The isotropic component $B_0$ is the yield in the window on the fit's
intensity scale. If only sampled profiles are available, approximate
(4) by integrating their piecewise-linear curves with the actual
energy spacing and splitting segments at the window edges.
Use the same integration rule for every order.

Write $B_\ell=\mathbf q_\ell^\mathsf T\widehat{\mathbf c}$, where
$\mathbf q_\ell$ holds the corresponding radial integrals in the
order-$\ell$ block. Set

$$
\mathbf a_\ell=F^\mathsf T H^\mathsf T\mathbf q_\ell,\qquad
C_{\ell m}=\operatorname{Cov}(B_\ell,B_m)
=\sum_i a_{\ell i}a_{mi}\mu_i.
\tag{5}
$$

Apply the transposes to the requested observables to compute (5)
without forming a full detector covariance matrix. This also gives
covariance between different energy windows by pairing their respective
vectors $\mathbf a_\ell$. For a point on the spectrum, replace the
radial integrals in $\mathbf q_\ell$ with the factors in (3) evaluated
at the requested energy.

## 3. Anisotropy error bars

Write $m_\ell=\mathbb E[B_\ell]=\mathbf a_\ell^\mathsf T\boldsymbol\mu$
for the mean fitted yield of order $\ell$. For $m_0>0$, the ratio of mean
yields is $\beta_\ell^\star=m_\ell/m_0$. A single image gives the estimate
$\widehat\beta_\ell=B_\ell/B_0$ wherever $B_0\ne0$. Expanding this ratio
about $(m_\ell,m_0)$ gives

$$
\widehat\beta_\ell-\beta_\ell^\star
\simeq\frac{(B_\ell-m_\ell)-\beta_\ell^\star(B_0-m_0)}{m_0}.
\tag{6}
$$

For two anisotropy parameters from the same window, this gives

$$
\operatorname{Cov}(\widehat\beta_\ell,\widehat\beta_m)
\simeq
\frac{
C_{\ell m}-\beta_m^\star C_{\ell0}-\beta_\ell^\star C_{0m}
+\beta_\ell^\star\beta_m^\star C_{00}
}{m_0^2}.
\tag{7}
$$

For one image, substitute $B_0$, $\widehat\beta_\ell$, $\widehat\beta_m$
and the estimated yield covariance into the right-hand side of (7).
Set $m=\ell$ and take the square root for the estimated standard error.
For a pointwise beta, use the same formula with $B_\ell=I_\ell(E_i)$ and
$B_0=I_0(E_i)$ at the same sampled energy. For several windows, apply the
Jacobian from (6) to their joint covariance.

The expansion requires denominator fluctuations to be small relative to
$m_0$. Compare the measured $B_0$ with its standard error
$\sqrt{\widehat C_{00}}$. If the uncertainty in the isotropic yield
includes zero, a symmetric beta error bar is unreliable even when the
measured yield is positive. Basis truncation and regularization can also
shift $\beta_\ell^\star$ away from the physical anisotropy.

## 4. Experimental inputs and validation

Under (1), a single original integer count image gives

$$
\widehat C_{\ell m}=\sum_i a_{\ell i}a_{mi}n_i.
$$

Since $\mathbb E[n_i]=\mu_i$, this estimator of the linear yield covariance
is unbiased at every count level. Sparse counts can make it noisy. A pixel
with zero observed events contributes zero in that exposure even if its
mean is positive. Taking a square root or forming a ratio is nonlinear,
so the unbiasedness result for $\widehat C_{\ell m}$ does not carry over.

If acquisition stops after exactly $N$ detected electrons, independent
events with pixel probabilities $\mathbf p$ give multinomial counts.
Then $\boldsymbol\mu=N\mathbf p$ and the input covariance is
$N[\operatorname{diag}(\mathbf p)-\mathbf p\mathbf p^\mathsf T]$.
Propagating this through the same fixed maps gives

$$
C_{\ell m}^{(N)}
=\sum_i a_{\ell i}a_{mi}\mu_i-\frac{m_\ell m_m}{N}.
$$

The correction affects absolute yields. For a ratio from the same count
image, it cancels in (7) because $m_\ell-\beta_\ell^\star m_0=0$.
This cancellation holds to first order for the fixed linear analysis above.

### Experimental checks

The package tests draw 1,000 independent Poisson histograms per geometry
from a known forward projection and refit each histogram. The quadrant
case tests an even anisotropy. The half-plane case includes both odd and
even anisotropies and a cropped detector edge. For the spectrum samples,
pointwise betas, integrated yield, and integrated betas, the tests compare the root mean
square reported error with the observed standard deviation. They also
measure how often the estimate plus or minus one standard error contains
the result from fitting the noiseless mean image.

- For independent Poisson signal and background exposures with subtraction
  $\mathbf n_{\rm sig}-a\mathbf n_{\rm bg}$ and fixed exposure scale
  $a$, replace $D_{\mathbf n}$ in (2) with
  $\operatorname{diag}(\boldsymbol\mu_{\rm sig}
  +a^2\boldsymbol\mu_{\rm bg})$. If the two images have different
  preprocessing, propagate their contributions separately.
- Validate the calculation with Poisson draws of a detector projection
  with known even and odd coefficients. Repeat the analysis and
  compare the spread of the recovered betas with (7). Fit the noiseless
  mean projection to obtain $\beta_\ell^\star$ and check how often the
  interval covers it. Comparing $\beta_\ell^\star$ to the input physical
  beta separately measures reconstruction bias. With shot-resolved data,
  resampling whole shots or acquisition blocks can capture correlations
  beyond the independent-Poisson model.
- Quote center, energy calibration and detector-response uncertainties
  separately from the fixed-setting counting error. If the image itself
  determines the center or the windows, repeat that estimation step in
  resampled analyses to measure its effect.
