# Counting errors in pBASEX

This derivation uses the notation of the
[original pBASEX note](pbasex.pdf) and the
[half-plane derivation](forward-backward-asymmetry.md). It gives the
counting uncertainty of yields and anisotropies for fixed analysis
settings. It applies to both detector geometries.

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
For several windows, apply the Jacobian from (6) to their joint covariance.

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
