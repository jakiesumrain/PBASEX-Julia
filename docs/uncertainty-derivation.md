# From detected electrons to pBASEX error bars

Suppose we repeat a measurement with the same detector, the same exposure and the same analysis settings. How much would the height of a spectral peak change? How much would its angular anisotropy change? A **standard error** estimates the standard deviation of such repeated estimates. It has the same units as the quantity we report. We can estimate it from a single acquisition if we specify what varies between repetitions.

An accumulated image tells us how many electrons reached each pixel. A shot checkpoint also tells us which electrons arrived together. We first derive the quantity both measurements estimate, then use the information in each kind of record to calculate its uncertainty. The equation numbers cited as “PDF” refer to the [original pBASEX derivation](pbasex.pdf).

## 1. What a pixel measures

Use the PDF's axes. The detector lies in the $yz$ plane, $x$ runs along the line of sight, and $z$ is the physical symmetry axis. Define $R^2=y^2+z^2$, $r^2=x^2+R^2$ and $\cos\theta=z/r$. The three-dimensional density $f(r,\theta)$ gives the expected number of detected particles per unit volume in the reconstructed distribution. Axial symmetry means that turning about $z$ leaves this density unchanged. A detector cannot distinguish two electrons at different $x$ that arrive at the same $(y,z)$. It measures their combined projection (PDF 10–12),

$$
\widetilde f(y,z)=\int_{-\infty}^{\infty} f(r,\theta)\,dx
=2\int_0^\infty f(r,\theta)\,dx .
$$

Changing the sign of $x$ leaves $r$ and $\cos\theta$ unchanged. This explains the factor of two even when emission toward $+z$ differs from emission toward $-z$. Each raw histogram entry $n_i$ counts events in the area of one detector bin. The model samples the projected basis on the detector grid.

The PDF describes the angular dependence with Legendre polynomials $P_\ell$. For example, $P_0(u)=1$, $P_1(u)=u$ and $P_2(u)=(3u^2-1)/2$. Choose radial functions $f_k(r)$ indexed by $k$ and fit their amplitudes $c_{k\ell}$ (PDF 14–19):

$$
f(r,\theta)=\frac{1}{4\pi}\sum_\ell\sum_k c_{k\ell}f_k(r)P_\ell(\cos\theta).
$$

The isotropic term gives $\sigma(r)=\sum_k c_{k0}f_k(r)$. Each higher-order amplitude gives $\sigma(r)\beta_\ell(r)=\sum_k c_{k\ell}f_k(r)$. Thus $\beta_\ell$ measures an angular modulation relative to the isotropic component. Under $z\mapsto-z$, $P_\ell$ acquires a factor $(-1)^\ell$. Quadrant folding averages the two axial directions and requires even orders. Half-plane folding averages left and right while keeping the two axial directions separate, so odd orders survive.

## 2. From the projection to the fitted coefficients

Think of a basis image as a prediction for a single coefficient set to one. Project each term in the expansion along $x$ (PDF 20–27). The entry in row $p$, column $(k,\ell)$ of the forward matrix is the sampled projection

$$
G_{p,(k,\ell)}=\frac{1}{2\pi}\int_0^\infty
 f_k\!\left(\sqrt{R_p^2+x^2}\right)
 P_\ell\!\left(\frac{z_p}{\sqrt{R_p^2+x^2}}\right)\,dx .
$$

Let $\vec n$ list the original pixel counts. Let $F$ perform the fixed rotation, pixel selection and fold so that $\vec b=F\vec n$ lists the fitted detector pixels. Here $F$ denotes a matrix; the PDF's $F(E,\theta,\phi)$ denotes a density. If two mirror pixels contribute equally, one row of $F$ contains two entries of $1/2$. On a symmetry axis a pixel can be its own mirror. Its contributions add in the **same** column of $F$. This distinction matters when we square them to calculate a variance.

Each row of $G\vec c$ predicts one entry of $\vec b$. For $n_p$ fit pixels and $n_b$ basis functions, $G$ has size $n_p\times n_b$ and $\vec c$ has $n_b$ entries. Counting noise generally prevents an exact solution. The fit minimizes the sum of squared residuals, optionally adding the penalty $\lambda\|\vec c\|^2$, where $\lambda$ is Julia's scalar `reg`. Write $G=U\Sigma V^\mathsf T$ as in PDF 28–34. The columns of $U$ and $V$ are orthonormal directions in pixel space and coefficient space. A singular value $\sigma_j$ tells how much the forward projection amplifies direction $j$. This $\sigma_j$ is distinct from the radial function $\sigma(r)$.

Set $d_j=(V^\mathsf T\vec c)_j$ and $t_j=(U^\mathsf T\vec b)_j$. In these coordinates the part of the objective that depends on $d_j$ is $(\sigma_j d_j-t_j)^2+\lambda d_j^2$. Setting its derivative to zero gives $(\sigma_j^2+\lambda)d_j=\sigma_j t_j$. Rotating the fitted coordinates back gives

$$
\widehat{\vec c}=H\vec b=HF\vec n,\qquad
H=V\,\mathrm{diag}(h_j)\,U^\mathsf T,\qquad
h_j=\frac{\sigma_j}{\sigma_j^2+\lambda}\quad(\lambda>0).
$$

For $\lambda=0$ the code uses $1/\sigma_j$ above its numerical cutoff and zero below or at it. A small singular value would otherwise turn a small detector fluctuation into a large coefficient fluctuation. With a fixed basis, geometry and regularization, $HF$ is a fixed linear operation. This is the property both error calculations use.

## 3. What the spectrum and a window measure

Consider a thin spherical shell of thickness $dr$. Its volume element is $r^2\,dr\,d\Omega$, with $d\Omega=\sin\theta\,d\theta\,d\phi$. Integrating over all directions gives $4\pi$ for $P_0$ and zero for each higher order. The shell contains $dN=\sigma(r)r^2dr$ fitted electrons. Energy calibration gives $E=\alpha r^2$ and $dE=2\alpha r\,dr$. Dividing the shell yield by its energy width gives PDF 45,

$$
I(E)=\frac{dN}{dE}=\frac{r}{2\alpha}\sum_k\widehat c_{k0} f_k(r),
\qquad r=\sqrt{E/\alpha}.
$$

Keeping the angular term of order $\ell$ in the same change of variables gives $I_\ell(E)=\frac{r}{2\alpha}\sum_k\widehat c_{k\ell}f_k(r)$. Thus $I_0=I$, and $I_\ell/I_0$ gives PDF 38 wherever the denominator is nonzero. Only the isotropic component counts all directions. Higher components carry signed information about the shape of the distribution. With $\alpha=1$ in detector-pixel units the spectrum has units of detected electrons per px². An interval integral has units of electrons. Integrating the fitted density gives a model yield, which can differ from the number of detector hits in a ring because several three-dimensional shells project onto the same ring.

For a window $w=[E_a,E_b]$, integrate each component before dividing:

$$
B_{w\ell}=\int_{E_a}^{E_b}I_\ell(E)\,dE,\qquad
Y_w=B_{w0},\qquad \overline I_w=\frac{Y_w}{E_b-E_a},\qquad
\widehat\beta_{w\ell}=\frac{B_{w\ell}}{B_{w0}}.
$$

Since $I_\ell=I_0\beta_\ell$ wherever the ratio exists, the integrated beta averages the angular modulation according to the electrons contributed at each energy. An arithmetic average of sampled betas generally gives a different value.

**How the code integrates.** The native samples are $E_j=\alpha r_j^2$. An evenly spaced radius therefore produces unequal energy intervals. Between two samples the code draws a straight line. On the overlap $[u,v]$ of the window with $[E_j,E_{j+1}]$, set $t(E)=(E-E_j)/(E_{j+1}-E_j)$ and $\bar t=(t(u)+t(v))/2$. The integral of this line is exactly

$$
(v-u)\big[(1-\bar t)I_\ell(E_j)+\bar t I_\ell(E_{j+1})\big].
$$

Sum over nonempty overlaps. For example, a window $[2,4]$ on the segment $[1,4]$ uses $t(u)=1/3$, $t(v)=1$ and contributes $(2/3)I_\ell(1)+(4/3)I_\ell(4)$. Its two factors sum to its width, 2. Adjacent windows therefore preserve the integral of the interpolated spectrum. This is exact for the interpolant and approximates the continuous basis integral.

Collect the resulting sample factors into $d_{wj}$. All depend only on energies and limits. Substituting the expression for $I_\ell$ gives a vector $\vec q_{w\ell}$ with entries

$$
(q_{w\ell})_{km}=
\begin{cases}
\displaystyle\sum_j d_{wj}\frac{r_j}{2\alpha}f_k(r_j), &m=\ell,\\
0,&m\ne\ell .
\end{cases}
$$

Consequently $B_{w\ell}=\vec q_{w\ell}^{\,\mathsf T}\widehat{\vec c}$. For a native sample use just that sample's radial factor in place of the sum. This lets us derive one uncertainty formula for both kinds of observation.

## 4. Counting errors from one accumulated image

Consider one window and order and abbreviate its vector as $\vec q$. Combining the fit and the integral gives a weighted sum over the **original** detector pixels,

$$
B=\vec q^{\,\mathsf T}HF\vec n=\vec a^{\,\mathsf T}\vec n
=\sum_i a_i n_i,\qquad \vec a=F^\mathsf T H^\mathsf T\vec q.
$$

Adding one detected electron to pixel $i$ changes the fitted observable by $a_i$. The transpose simply reads the same product from the other end. It asks how each input pixel contributes to one selected output. An inverse fit can give fractional or negative responses.

Assume that independent electrons arrive with fixed mean intensity during the exposure, their total is Poisson, and each is assigned to exactly one bin. Then $n_i$ is Poisson with mean $\mu_i$, variance $\mu_i$, and counts in distinct pixels are independent. Write $\delta n_i=n_i-\mu_i$. For an observable $B$ its fluctuation about its mean is $\sum_i a_i\delta n_i$. Squaring and averaging removes cross terms from different pixels because $\mathbb E[\delta n_i\delta n_j]=0$ for $i\ne j$. For two observables with responses $\vec a$ and $\vec a'$ this gives

$$
\mathrm{Var}(B)=\sum_i a_i^2\mu_i,\qquad
\mathrm{Cov}(B,B')=\sum_i a_i a'_i\mu_i .
$$

Covariance measures how two estimates move together across repetitions. It need not vanish for different energies or angular orders because the same pixels feed the inversion. In one exposure substitute $n_i$ for the unknown $\mu_i$. Since $\mathbb E[n_i]=\mu_i$, the resulting estimate of the **linear covariance** is unbiased. The standard error is the square root of the estimated variance. For a spectrum bin, divide the yield SE by the positive bin width.

For example, two mirror pixels containing 9 and 16 electrons give a folded value $(9+16)/2=12.5$ and estimated variance $(9+16)/4=6.25$. Its SE is 2.5. A single axis pixel with 9 electrons, encountered twice in the fold, gives $(9+9)/2=9$. The two terms represent the same measured fluctuation, so its SE is 3. The code combines duplicate contributions first. Distinct folded positions then use disjoint raw pixels, so the implementation can apply the folded pixel standard deviations to $H^\mathsf T\vec q$ and take the vector norm. This computes the same sum of squares without storing a detector covariance matrix.

**Error of a ratio.** Write $X=B_{w\ell}$, $Y=B_{w0}$ and let $x,y$ be their expectations, with $y>0$. For small fluctuations about these means, expand the quotient to first order:

$$
\frac{x+\delta X}{y+\delta Y}
\simeq \frac{x}{y}+\frac{\delta X}{y}-\frac{x\,\delta Y}{y^2}.
$$

Set $\beta^\star=x/y$. The fluctuation in the ratio is approximately $(\delta X-\beta^\star\delta Y)/y$. Squaring and averaging gives $(\mathrm{Var}(X)-2\beta^\star\mathrm{Cov}(X,Y)+(\beta^\star)^2\mathrm{Var}(Y))/y^2$. Substitute the observed $Y$, $\widehat\beta=X/Y$ and pixel counts to obtain the implemented error:

$$
\widehat{\mathrm{SE}}(\widehat\beta)
=\frac{\sqrt{\sum_i n_i(a_{\ell i}-\widehat\beta\,a_{0i})^2}}{Y},
\qquad Y>\widehat{\mathrm{SE}}(Y).
$$

This subtraction retains the common noise in numerator and denominator. For a native beta use the responses at a single energy. As $Y$ approaches zero, even a small fluctuation can change the sign of the ratio or make its magnitude enormous. The first-order expansion then fails. Julia returns `NaN` for the symmetric beta SE when $Y\le\widehat{\mathrm{SE}}(Y)$. Exceeding that threshold permits an estimate; it does not guarantee a Gaussian ratio distribution.

## 5. What a shot ID adds

Let $s$ label one acquired shot, $\vec n_s$ its histogram and $S$ the number of shots, including shots with zero electrons. With one count per electron, $\vec n=\sum_s\vec n_s$. No independence assumption is needed to distribute a linear operation over this sum. The contribution of shot $s$ to window $w$, order $\ell$, is

$$
b_{sw\ell}=\vec a_{w\ell}^{\,\mathsf T}\vec n_s
=\sum_{e\ \mathrm{in\ shot}\ s}a_{w\ell,i(e)},\qquad
B_{w\ell}=\sum_s b_{sw\ell}.
$$

Here $i(e)$ is the bin of the corrected event. A fixed exposure multiplier $g_s$, if supplied, multiplies the entire shot contribution. Setting $g_s=1$ recovers detected-electron units. The scalar $b_{sw\ell}$ is one fitted component contributed by one shot; $\vec b$ in PDF 27 is the whole detector vector.

The code evaluates $H^\mathsf T\vec q=U\,\mathrm{diag}(h_j)V^\mathsf T\vec q$ in blocks of singular vectors. It looks up each event's folded position and fold fraction, then adds the corresponding response to that shot's entry in a table. Its rows represent window/order pairs and its columns represent shots. An empty shot supplies a column of zeros. This table gives the same pooled linear observables as inverting the pooled histogram for matching settings. Keeping its columns separate also keeps correlations between different pixels, energies and orders **within each shot**.

**A repeated acquisition from the observed shots.** For independent, exchangeable shots, approximate the unknown population of shots by the empirical population assigning probability $1/S$ to each observed shot. Draw $S$ times with replacement. If shot $s$ is drawn $m_s^{(b)}$ times in simulated run $b$, then $\sum_s m_s^{(b)}=S$. Every observable in that run uses the same multiplicities:

$$
B_{w\ell}^{*(b)}=\sum_s m_s^{(b)}b_{sw\ell},\qquad
Y_w^{*(b)}=B_{w0}^{*(b)},\qquad
\beta_{w\ell}^{*(b)}=B_{w\ell}^{*(b)}/Y_w^{*(b)}
\quad (Y_w^{*(b)}>0).
$$

This follows by inserting $\sum_s m_s^{(b)}\vec n_s$ into the same inverse. Hence a bootstrap draw gives the result of a pooled refit using just sums and a division. Averaging individual-shot beta values would give equal influence to bright and dim shots and would require dividing by zero for an empty shot.

For example, suppose three shots contribute $(B_0,B_2)=(10,2),(30,12),(0,0)$. The observed yield is 40 and $\widehat\beta_2=14/40=0.35$. Drawing the first shot twice and the empty shot once gives yield 20 and beta $4/20=0.20$. Drawing the empty shot three times gives a valid zero yield and an undefined beta. Removing empty shots from the roster would change the empirical population.

**From repeated draws to a reported error.** For $M$ bootstrap yields $Y^{*(b)}$, first calculate their mean $\overline Y^*=M^{-1}\sum_b Y^{*(b)}$. Then compute

$$
\widehat{\mathrm{SE}}_{\rm shot}(Y)=
\sqrt{\frac{1}{M-1}\sum_{b=1}^M (Y^{*(b)}-\overline Y^*)^2}.
$$

Each draw already represents an entire acquisition of the original size. Its spread estimates the uncertainty of the measured total. Dividing once more by $\sqrt M$ would estimate the numerical error of the bootstrap mean instead. A pointwise 95% percentile interval takes the 2.5th and 97.5th percentiles of the draws. Julia reports the original pooled central value. For beta it computes the SD and percentiles on positive-yield draws, provided the pooled and lower percentile yields are positive and at least two valid draws remain. The fraction of excluded draws indicates how often the denominator failed.

**A check on the bootstrap normalization.** For one run with independent shots, let $t_s=b_{sw0}$, and let $\bar t=S^{-1}\sum_s t_s$. Conditional on the observed table, one draw of a single shot has mean $\bar t$ and variance $S^{-1}\sum_s(t_s-\bar t)^2$. A simulated acquisition adds $S$ independent selections, so the variance of its total is

$$
\mathrm{Var}_*(Y^*)=S\,\frac1S\sum_s(t_s-\bar t)^2=\sum_s(t_s-\bar t)^2 .
$$

For the three-shot example this gives $1400/3$ electrons squared, or an SE of about 21.6 electrons. With many bootstrap draws, the reported yield SE approaches this value for that table. The finite empirical population matters: estimating the unknown population variance of a single shot with the usual $S-1$ denominator would give a total-variance estimate larger by $S/(S-1)$. The implementation reports the empirical bootstrap spread. Increasing the number of draws samples it more precisely; collecting more independent shots improves how well the table represents the population.

**Shots correlated in time.** If $T_s$ denotes a random yield from shot $s$, expanding the variance of its sum gives $\mathrm{Var}(\sum_s T_s)=\sum_s\mathrm{Var}(T_s)+2\sum_{s<t}\mathrm{Cov}(T_s,T_t)$. Suppose the FEL intensity stays high for several adjacent shots. Sampling individual shots would break these sequences and lose their common fluctuation. With `block_length=L`, draw a starting position uniformly from the full length-$L$ blocks within a run, copy the consecutive shot columns together, and repeat until the original run length is reached. Shorten the last block as needed. Apply the same procedure independently within each specified run. This approximates repetition of a run with short-range dependence when blocks capture the correlation and the process is approximately stationary within each run. Very long blocks leave few distinct sequences. A sustained drift across the run needs separate experimental consideration.

**Comparing observables.** Label two chosen windows $u$ and $v$ and fix one order $\ell$. Set $X=\widehat\beta_{u\ell}$ and $Y=\widehat\beta_{v\ell}$. Compute $D^{*(b)}=Y^{*(b)}-X^{*(b)}$ for each paired draw with both ratios defined. The spread of $D^*$ includes the covariance term in $\mathrm{Var}(Y-X)=\mathrm{Var}(Y)+\mathrm{Var}(X)-2\,\mathrm{Cov}(X,Y)$. This works for any two windows. The individual uncertainties alone cannot recover their covariance.

## 6. Which uncertainty belongs on the plot?

For a pooled raw histogram, `pbasex_counts` answers what would fluctuate if independent electrons sampled the same mean detector distribution again. For a checkpoint, `pbasex_shots` answers what would fluctuate if another acquisition of the same number of shots sampled the observed shot population or its within-run blocks. A change in pulse intensity can move many pixel counts together. That changes the shot variance even when both calculations produce the same fitted spectrum. Counting fluctuations already appear in the observed shot spread, so adding the two reported SEs in quadrature would count that contribution twice.

Both calculations condition on the fitted model, calibration, center, distortion correction and window limits. For example, if every event is shifted by the same erroneous center, every bootstrap run inherits that shift. A narrow error bar measures repeatability under the stated assumptions. It cannot establish the absolute accuracy of the reconstructed energy or angular distribution.

[Count-image examples](counting-errors.md) · [Shot-checkpoint examples](shot-uncertainties.md) · [Julia usage guide](usage.md)
