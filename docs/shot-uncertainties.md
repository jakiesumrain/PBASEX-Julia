# How shot-resolved pBASEX works

[From detected electrons to pBASEX error bars](uncertainty-derivation.md)
derives the pixel response used below from the projected basis in the
original pBASEX note. Jump to [Running the calculation](#5-running-the-calculation)
for the Julia inputs, settings and outputs.

Imagine repeating a run of 2,763 FEL shots under the same acquisition
conditions. How much would the fitted electron yield and angular anisotropy
change? `pbasex_shots` estimates this spread using the shots you actually
recorded. It keeps events from the same shot together. A bright shot, an
empty shot, or a shot with an unusual angular distribution enters a simulated
repeat as a whole.

The calculation has two parts. First, determine how much **each shot**
contributed to the fitted angular components in each energy interval.
Second, draw shots with replacement, sum those contributions, and form the
observables again. The first part is linear. The division that produces
each beta happens after the sum.

## 1. A three-shot example

Suppose the inversion assigns the following contributions to one interval.
The isotropic component is the fitted yield. The order-two component is the
signed numerator of $\beta_2$.

| Shot | Isotropic component | Order-two component |
| --- | ---: | ---: |
| A | 10 | 2 |
| B | 30 | 12 |
| C | 0 | 0 |

For the recorded run the yield is $40$ and $\widehat\beta_2=14/40=0.35$.
One possible bootstrap draw is **A, A, C**. Its yield is $20$ and its beta
is $4/20=0.20$. Another is **A, B, B**, giving $70$ and $26/70\approx0.371$.
The draw **C, C, C** gives zero yield and an undefined beta. This illustrates
why the code retains empty shots, and why it takes the ratio *after* pooling:
the bright shot contributes more electrons to the observed distribution.

## 2. Where the table entries come from

Let $\mathbf n_s$ count corrected events from shot $s$ in the detector pixels.
For each event the reader locates the column and row using the checkpoint's
stored bin edges. Edges belong to the bin on their right except that the
last bin includes its upper edge. The complete shot roster determines the
columns of the shot table, including columns containing only zeros.

Use the notation of the [counting-error derivation](counting-errors.md).
$F$ rotates, folds and selects the detector pixels. $H=V\,\mathrm{diag}(h)U^\mathsf T$
is the inverse used by the fit, with the same singular-value cutoff or
regularization as `pbasex`. For fixed shot multipliers $g_s$ (default 1),
the pooled coefficients are

$$
\widehat{\mathbf c}=HF\sum_s g_s\mathbf n_s .
$$

For a window $w=[E_a,E_b]$, let $\mathbf q_{w\ell}$ integrate the angular
energy component $I_\ell(E)$ over that window. The implementation gets
these integrals by linearly interpolating between the native energy
samples, including partial segments at the two edges. If the checkpoint
has $d$ histogram bins per detector pixel and $x$ is radius in histogram
bins, the code uses $E=(\alpha/d^2)x^2$ and the corresponding spectrum
factor $x/[2(\alpha/d^2)]$.

The same fitted integral can now be written in two ways:

$$
B_{w\ell}
=\mathbf q_{w\ell}^{\mathsf T}\widehat{\mathbf c}
=\sum_s g_s\,\underbrace{(F^\mathsf T H^\mathsf T\mathbf q_{w\ell})^\mathsf T}
_{\text{response to each detector pixel}}\,\mathbf n_s .
$$

The vector $\mathbf a_{w\ell}=F^\mathsf T H^\mathsf T\mathbf q_{w\ell}$ answers
a concrete question: **if I add one electron to pixel $p$, by how much does
the fitted component in this window change?** The answer is $a_{w\ell,p}$.
These responses can be fractional or negative because the inverse fit
combines detector measurements to separate radial and angular components.
They are different for different windows and Legendre orders.

For the list of pixels hit in shot $s$, its table entry is simply

$$
b_{sw\ell}=g_s\sum_{e\in s} a_{w\ell,p_e}.
$$

The code computes $U\,\mathrm{diag}(h)V^\mathsf T\mathbf q$ on the folded
grid. A lookup supplies the folded index and the fraction contributed by
each original pixel. For example, half-plane folding gives each of two
mirror pixels half the response of their common folded position. A lone
measured pixel at that position contributes the full response. Quadrant
folding accounts for repeated occurrences of pixels on the symmetry axes.
This calculation uses the deterministic inverse response. The Poisson
variance factors used by `pbasex_counts` enter a different calculation.

**Consequence.** For any selected multiset of shots, summing their table
entries gives the same integrals as rebuilding their detector histogram
and applying the fixed inverse. This identity explains both the speed of
the bootstrap and the direct-refit verification below.

## 3. What is reported for each interval

Sum over the recorded shots to get $B_{w\ell}=\sum_s b_{sw\ell}$. The result
fields follow directly.

| Field | Definition | Units with unit event counts |
| --- | --- | --- |
| `yield` | $B_{w0}$ | Fitted detected electrons in the interval |
| `spectrum` | $B_{w0}/(E_b-E_a)$ | Detected electrons per unit energy |
| `betas[j]` | $B_{w\ell}/B_{w0}$ for a positive yield | Dimensionless |

The beta vector follows `analysis.legendre_orders[2:end]`. The spectrum is
a bin average. With `alpha=1` the energy coordinate is detector radius
squared in px². Equal-width adjacent windows give a uniform-grid spectrum.
With fixed nonunit shot multipliers the yield is a weighted yield.

## 4. How a bootstrap draw becomes an error bar

Let $m_s^{(b)}$ count how often shot $s$ appears in draw $b$. For $S$
independent exchangeable shots the code draws $S$ times with replacement.
For each window it computes

$$
B^{*(b)}_{w\ell}=\sum_s m_s^{(b)}b_{sw\ell},\qquad
\beta^{*(b)}_{w\ell}=B^{*(b)}_{w\ell}/B^{*(b)}_{w0}.
$$

All windows and orders use the **same multiplicities**. If one bright shot
feeds two nearby peaks, drawing it twice raises both contributions in the
same replicate. Pairing the replicate values preserves this covariance
for a downstream difference or ratio.

- **Standard errors.** `yield_se` is the sample standard deviation of the
  bootstrap totals. `spectrum_se` divides it by the window width. The
  replicate already represents a full run of $S$ shots, so there is no
  further division by $\sqrt S$.
- **Intervals.** For confidence $1-\gamma$, the code takes the empirical
  $\gamma/2$ and $1-\gamma/2$ quantiles. The reported central value comes
  from the original pooled shots. Intervals are pointwise across windows.
- **Weak yields.** Each draw with nonpositive isotropic yield has an
  undefined beta. `invalid_fraction` reports its frequency. Beta SEs and
  percentile intervals use positive-yield draws only, and are suppressed
  if the pooled yield or the lower yield interval endpoint is nonpositive,
  or if fewer than two valid draws remain. Near zero, even a finite ratio
  interval warrants care.
- **Serial dependence.** With `block_length=L`, the code samples contiguous
  blocks of $L$ shots within each consecutive run, trimming the final block
  to preserve run length. This preserves short-range correlations inside
  the sampled blocks. Acquisition order and the choice of $L$ therefore
  matter. A run of exactly $L$ shots has only one possible block and
  contributes no resampling variation. If this holds for every run,
  the routine raises an error.

These uncertainties cover variation represented in the observed shots.
Center, calibration, exposure multipliers, regularization, and energy
limits stay fixed. A bias shared by every shot can persist across every
bootstrap draw. Counting noise is already present in the shot variation.

## 5. Running the calculation

```julia
using PBASEX, Random

analysis = pbasex_shots(checkpoint_path; basis=basis_path,
    windows=[(Ea, Eb), (Ec, Ed)], rotation=90,
    replicates=2000, rng=Xoshiro(42), keep_replicates=true)

w = analysis.windows[1]
w.window, w.yield, w.yield_se, w.yield_ci
w.spectrum, w.spectrum_se, w.betas, w.beta_se, w.beta_ci
w.valid_replicates, w.invalid_fraction
```

To compare with a pooled histogram, use the same windows and fit settings.
Here `counts` contains one count per corrected event, binned in the original detector orientation at one bin per detector pixel.

```julia
settings = (; basis=basis_path, geometry=:half_plane, windows=[(Ea, Eb)],
            center=(column, row), rotation=90)
pooled = pbasex_counts(counts; settings...)
shots = pbasex_shots(checkpoint_path; settings...)
for analysis in (pooled, shots)
    window = only(analysis.windows)
    @show window.yield window.yield_se window.spectrum window.spectrum_se window.betas window.beta_se
end
```

With matching settings and positive fitted yield, the central values agree
to numerical precision. For a numerical comparison, supply integer counts
or a `Float64` histogram. The shot fit uses double precision. A `Float32`
histogram uses single precision in `pbasex_counts`, including its
precision-dependent singular-value cutoff at `reg=0`.
Each estimator supplies its own standard errors.
The shot estimator also supplies bootstrap intervals and diagnostics.

The basis selects the geometry for shot analysis. Supply `geometry` to check
that it matches. The center override is one-based
`(column,row)` on the unrotated histogram. The default uses the checkpoint
center. Both estimators rotate counterclockwise for positive angles. Accepted
angles are multiples of 90°, so each electron stays in exactly one bin.
Quadrant filters refer to the rotated image. For a common
spectrum confidence level, divide both endpoints of `yield_ci` by the
positive window width. Retained draws use `w.replicates.betas[j,b]` for
order index $j$ and draw $b$.

The reader streams columns of the large HDF5 matrix $U$ and batches of
events. It retains a folded response map per requested component and a
table of component totals per shot for the current window batch. The same
bootstrap seed is replayed for each batch. Thus changing the batch size
preserves the pairing. The default processes one window per pass. See
the [checkpoint schema and memory controls](counting-errors.md#shot-resolved-uncertainties-from-a-corrected-checkpoint).

### Choose settings for an experimental run

| Keyword | Meaning |
| --- | --- |
| `windows` | Required. Supply one or more `(low, high)` pairs within the native energy range in the units of `alpha`. Windows can overlap and their limits can lie between samples. Results follow input order. |
| `alpha`, `reg` | Default to 1 and 0. For shots, `alpha` is energy per squared **detector pixel**. For counts, it is energy per squared **histogram bin**. If the checkpoint has $d$ bins per detector pixel, use `alpha_counts = alpha_shots / d^2` for its pooled histogram. Use the same regularization in both fits. |
| `geometry`, `quadrant_filter` | Shots infer geometry from the basis. Counts default to `:quadrant`. Set `geometry=:half_plane` in shared settings for half-plane fits. Four entries select top-right, top-left, bottom-left and bottom-right **after** rotation. The default selects all four. |
| `replicates`, `confidence`, `rng` | Defaults are 1000 draws, 0.95 and Julia's default RNG. Supply `rng=Xoshiro(42)` for reproducible draws. Increasing the number of draws reduces Monte Carlo variation in estimated SEs and percentiles. |
| `run_ids`, `block_length` | With no run IDs the roster is one run. A change of ID starts a new run. The default block length of 1 draws individual shots. Use a longer block for serial correlation and retain the acquisition order. |
| `shot_weights` | Optional fixed nonnegative multipliers, one per roster entry; default 1 each. With nonunit multipliers the yield and spectrum are on the chosen normalization scale. |
| `window_batch`, `event_chunk`, `basis_chunk_bytes` | Defaults are 1 window per batch, 65,536 events per read and 64 MiB targeted per read of singular-vector columns. Increase the first setting to reduce repeated reads if memory permits. |
| `keep_replicates`, `replicate_sink` | By default the return value holds summaries. Set `keep_replicates=true` to retain draws, or pass a function `(window_index, draws) -> ...` to consume them window by window. |

The return value records `shots` (number of acquired shots), `geometry`,
`legendre_orders`, `bootstrap_replicates`, `block_length`, `confidence`
and `windows`. Each window includes the shared fields `window`, `yield`,
`yield_se`, `spectrum`, `spectrum_se`, `betas` and `beta_se`. For asymmetric
error bars plot the two endpoints of `yield_ci` or `beta_ci[j]`. To plot a
spectrum interval, divide both yield endpoints by the window width.
`valid_replicates` counts positive-yield draws and `invalid_fraction`
reports the fraction with undefined beta.

### A uniform spectrum from either input

Choose `edges` in energy units inside the fitted range. Both calls below
use the same grid and assume one histogram bin per detector pixel and unit
shot multipliers. A checkpoint call requires windows. A counts call with
windows also returns the native samples in `result`.

```julia
windows = collect(zip(edges[1:end-1], edges[2:end]))
fit = (; basis=basis_path, geometry=:half_plane, center=(column, row),
         rotation=90, windows)
pooled = pbasex_counts(counts; fit..., make_images=false)
shots = pbasex_shots(checkpoint_path; fit..., rng=Xoshiro(42))

energy = (edges[1:end-1] .+ edges[2:end]) ./ 2
function plot_values(a)
    j = only(findall(==(2), a.legendre_orders[2:end]))  # basis includes order 2
    spectrum = [w.spectrum for w in a.windows]
    spectrum_se = [w.spectrum_se for w in a.windows]
    beta2 = [w.betas[j] for w in a.windows]
    beta2_se = [w.beta_se[j] for w in a.windows]
    return (; spectrum, spectrum_se, beta2, beta2_se)
end
count_curve = plot_values(pooled)
shot_curve = plot_values(shots)
# Plot count_curve.spectrum against energy with count_curve.spectrum_se as error bars.
```

### A paired beta difference between any two intervals

Request both windows in one call so their draws use identical shot
multiplicities. For example, the code below compares order two in `w1`
and `w2`, each a pair of energy limits. Check the yield and beta diagnostics
for both windows. The interval for the difference uses draws where both
ratios exist.

```julia
using Statistics

a = pbasex_shots(checkpoint_path; fit..., windows=[w1, w2],
                 keep_replicates=true, rng=Xoshiro(42))
first_window, second_window = a.windows
j = only(findall(==(2), a.legendre_orders[2:end]))
if isfinite(first_window.beta_se[j]) && isfinite(second_window.beta_se[j])
    draws = second_window.replicates.betas[j, :] .- first_window.replicates.betas[j, :]
    valid = filter(isfinite, draws)
    if length(valid) >= 2
        difference = second_window.betas[j] - first_window.betas[j]
        se = std(valid)
        tail = (1 - a.confidence) / 2
        interval = quantile(valid, [tail, 1 - tail])
        @show difference se interval length(valid)
    end
end
```

## 6. What the checks establish

**Algebra and implementation.** The regression tests independently build
pooled and resampled images and run explicit Abel refits. They compare
yields and betas for identical draws in quadrant and half-plane geometry,
including a truncated detector, an empty shot, fixed multipliers, and
moving blocks. Repeating a window across a batch boundary must reproduce
its draws. This checks the linear identity, folding conventions and pairing.

**Statistical interpretation.** In repeated synthetic acquisitions with
known anisotropy, the tests compare the spread between independent runs with
the reported bootstrap SE and measure percentile interval coverage. An
independent-Poisson case checks the counting-noise limit. A correlated case
adds a common angular fluctuation lasting three shots. Agreement with
direct refits tests the numerical implementation; repeated acquisitions
test the uncertainty model under the simulated conditions.

[Implementation](../src/shot_uncertainties.jl) · [Executable validation](../test/shot_uncertainties.jl)
