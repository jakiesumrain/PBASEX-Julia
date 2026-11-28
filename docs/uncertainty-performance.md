# Uncertainty performance on the 512-radius basis

The measurements use baseline commit `2fef038` and Julia 1.11.3 on Windows
with an Intel Core i5-7300HQ and eight OpenBLAS threads. Each run used its own
Julia process and ran serially. The warm repeats followed the first runs
without clearing the Windows file cache. Wall time includes estimator
compilation and computation but excludes process startup and writing results.
Peak RAM is the Julia process's Windows peak working set.

The runs used `gaussian_r512_k256_l10_sigma1p6.h5`,
`localizations_circularized.h5`, and `image_circularized_counts.npy`
as read-only inputs. The basis has six angular orders and a 3.0 GiB U matrix.
The checkpoint has 2,763 shots and 1,258,772 events. The count image has one
bin per detector pixel. Both estimates used the original one-based center
`(477, 474)`, a 90° counterclockwise rotation, and zero regularization.
The full grid has 300 windows of width 400 from 0 to 120,000 detector pixel².
The short grid uses the first eight windows. Shot runs used batches of eight
windows and 1,000 whole-shot bootstrap replicates with 95% confidence,
block length 1, and `Xoshiro(42)`. Replicates were retained for comparison.

## Shot bootstrap

`resident_u=true` reads U once and streams the checkpoint for each batch.
`reuse_shot_work=true` bins and folds the events once and reuses the same
bootstrap multiplicities across batches. Both options are independent and
default to `false`. The cached shot-work benchmark also enables resident U.

| Grid | Run | U | Wall (s) | Peak RAM (GiB) | U read (s) | U multiply (s) | Coordinate build (s) | Checkpoint and accumulation (s) | Bootstrap (s) |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 300 | First | Chunked | 291.69 | 0.81 | 61.69 | 111.70 | 1.37 | 98.73 | 2.91 |
| 300 | First | Resident | 118.69 | 3.86 | 1.12 | 26.70 | 1.39 | 72.98 | 2.78 |
| 300 | Warm | Chunked | 285.47 | 0.79 | 58.62 | 109.63 | 1.34 | 100.98 | 2.88 |
| 300 | Warm | Resident | 118.64 | 3.88 | 1.07 | 25.30 | 1.34 | 73.20 | 2.77 |
| 300 | Cached shot work | Resident | 59.38 | 3.94 | 5.04 | 25.34 | 1.44 | 8.30 | 3.16 |
| 8 | Short | Chunked | 20.40 | 0.73 | 1.84 | 3.35 | 0.71 | 3.06 | 0.07 |
| 8 | Short | Resident | 17.86 | 3.89 | 1.34 | 0.97 | 0.73 | 2.18 | 0.07 |
| 8 | Cached shot work | Resident | 16.84 | 3.91 | 1.05 | 0.59 | 0.70 | 1.39 | 0.12 |

Resident U reduced warm wall time by 58% and replaced 38 HDF5 passes
with one read. Checkpoint work accounted for 62% of the resident run.
The stage times exclude work such as checkpoint metadata setup and summary
construction, so they sum to less than the full wall time.

On the dense grid, the cache spent 1.59 seconds reading and folding events
once and 6.71 seconds accumulating them across batches. Preparing the
draws took 0.03 seconds and applying them took 3.13 seconds. The savings
came from checkpoint work. The measured bootstrap stage was slightly
slower with cached draws. Per-shot projection into basis coordinates is
a separate follow-up experiment.

### Paging

Windows memory counters were sampled about every two seconds. These are
system-wide paging counters and include activity from other processes and
file-backed memory. The run's peak RAM comes from the Julia process.

| Grid | Run | U | Minimum available (MiB) | Hard read operations | Pages in | Page write operations | Pages out |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 300 | First | Chunked | 6,920 | 1,382 | 29,528 | 0 | 0 |
| 300 | First | Resident | 4,366 | 179 | 1,121 | 12 | 2,903 |
| 300 | Warm | Chunked | 7,367 | 1,032 | 6,222 | 22 | 5,530 |
| 300 | Warm | Resident | 4,440 | 206 | 1,802 | 32 | 7,937 |
| 300 | Cached shot work | Resident | 4,003 | 10,880 | 199,603 | 0 | 0 |
| 8 | Short | Chunked | 7,766 | 26 | 143 | 0 | 0 |
| 8 | Short | Resident | 4,460 | 108 | 1,084 | 0 | 0 |
| 8 | Cached shot work | Resident | 4,054 | 98 | 645 | 0 | 0 |

## Poisson count propagation

`native_errors=false` skips the 512 native-grid uncertainty points and
keeps the native fit values in the result. With either setting, window
propagation multiplies up to eight windows at once and reuses its buffers.
The beta derivative subtracts the denominator's contribution at each
folded detector pixel before taking the norm, retaining covariance.

The window-only ablation skips the native-grid uncertainty loop in the
baseline and uses the original per-window multiplication. The batched runs
use the production implementation, with native propagation enabled or
disabled as shown.

| Grid | Native errors | Window propagation | Wall (s) | Peak RAM (GiB) | Basis load (s) | Fit (s) | Native errors (s) | Window errors (s) |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | On | Baseline | 130.97 | 4.57 | 6.77 | 2.18 | 108.03 | 3.13 |
| 8 | Off | Baseline | 21.88 | 3.89 | 6.60 | 2.17 | — | 3.91 |
| 8 | On | Batched | 117.97 | 4.57 | 7.71 | 2.73 | 90.35 | 1.62 |
| 8 | Off | Batched | 23.58 | 3.87 | 6.48 | 2.07 | — | 1.18 |
| 300 | On | Baseline | 225.54 | 4.57 | 6.46 | 2.09 | 97.02 | 108.39 |
| 300 | Off | Baseline | 142.49 | 3.86 | 7.22 | 2.57 | — | 122.37 |
| 300 | On | Batched | 193.71 | 4.57 | 6.41 | 2.01 | 117.65 | 53.91 |
| 300 | Off | Batched, initial | 90.64 | 3.87 | 28.13 | 2.71 | — | 42.79 |
| 300 | Off | Batched, repeat | 56.56 | 3.86 | 6.27 | 2.04 | — | 35.37 |

On the short grid, batching reduced the window propagation stage from
3.91 to 1.18 seconds with native errors off. The complete window-only
batched process was 1.70 seconds slower in this single fresh-process
comparison, which includes compilation.

The dense window-only repeat was 4.0× faster than the baseline with
native errors, with a peak RAM reduction of 0.71 GiB. Other Julia activity
was observed during the dense unbatched window-only run and initial
window-only batched run. The repeat sampled for other Julia processes
about every two seconds and detected none. The dense native-plus-batched
run also had higher system paging than the baseline.

### Gram matrix experiment

The experimental path built the noise-weighted Gram matrix in detector row
tiles. It evaluated each beta variance with the difference of the angular
and beta-scaled isotropic coordinates, preserving the covariance term.
Both prototypes retained direct detector propagation as a fallback for
quadratic forms close to numerical cancellation. The batched run needed
zero fallbacks on the supplied histogram.

| 300-window variant | Wall (s) | Peak RAM (GiB) | Gram build (s) | Total window errors (s) |
| --- | ---: | ---: | ---: | ---: |
| Detector propagation, batched repeat | 56.56 | 3.86 | — | 35.37 |
| Gram, per-window products | 89.58 | 3.87 | 13.67 | 68.20 |
| Gram, batched products | 31.74 | 3.90 | 8.30 | 8.98 |

Gram was evaluated as a Float64 benchmark prototype. The production
implementation uses batched detector propagation.

## Numerical agreement

Finite results passed the per-element tolerance
`abs(a - b) <= 1e-8 + 1e-8 * max(abs(a), abs(b))`. Undefined values matched
exactly. All count estimates for integrated yield, spectrum, and beta
matched exactly, as did the native-grid standard errors when requested.
The counts window-only ablation also matched all window standard errors
exactly.

| 300-window comparison to baseline | Max yield difference | Max spectrum difference | Max beta difference | Max yield SE difference | Max spectrum SE difference | Max beta SE difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Shots, resident U | 9.09e-12 | 2.84e-14 | 1.51e-14 | 7.96e-13 | 2.00e-15 | 1.17e-10 |
| Shots, resident U and cached shot work | 9.09e-12 | 2.84e-14 | 1.51e-14 | 2.16e-12 | 5.33e-15 | 1.09e-10 |
| Counts, batched with or without native errors | 0 | 0 | 0 | 1.14e-13 | 4.44e-16 | 8.88e-16 |
| Counts, both Gram prototypes | 0 | 0 | 0 | 3.07e-12 | 7.77e-15 | 2.59e-14 |

| 300-window shot comparison | Max yield percentile endpoint difference | Max beta percentile endpoint difference | Max paired yield difference | Max paired beta difference | Max paired adjacent yield contrast difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Resident U | 1.46e-11 | 2.98e-13 | 4.37e-11 | 3.69e-9 | 5.82e-11 |
| Resident U and cached shot work | 8.00e-11 | 3.30e-13 | 1.75e-10 | 3.45e-9 | 2.69e-10 |

Every shot window had the same number of valid bootstrap replicates and
the same invalid fraction. The dense runs matched 10 undefined beta
standard errors, 20 undefined beta interval endpoints, and 1,325
undefined paired beta draws. For counts, the native grid matched five
undefined beta estimates and 270 undefined beta standard errors when
native errors were requested. The eight-window comparisons also passed.
All 854 repository tests passed.
