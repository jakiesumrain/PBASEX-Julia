# Line-of-sight integration

A detector pixel at projected radius $R$ receives contributions from 3D
radii $r = \sqrt{R^2 + u^2}$, where $u$ is the distance along the line of sight.
For a radial basis $f_k(r)$ centered at $k$, the isotropic forward-matrix
entry, before division by $2\pi$, is

$$
J(R,k) = \int_0^\infty f_k\!\left(\sqrt{R^2 + u^2}\right)\,du.
$$

For angular order $\ell$, the integrand also contains
$P_\ell(z/\sqrt{R^2+u^2})$, where $z$ is the signed coordinate along the
vertical physical axis. Half-plane geometry samples both signs of $z$
using the same radial integration bounds. Reversing $z$ preserves even
columns and reverses odd columns, while the integrand remains even in the
line-of-sight coordinate $u$. [Half-plane forward operator](forward-backward-asymmetry.md#3-the-pbasex-approach)

At $R<k$, the line of sight crosses the shell at $u_* = \sqrt{k^2-R^2}$.
A Gaussian shell of width $\sigma$ contributes mostly near this distance.
For a central pixel, $R=0$ and $u_*=k$. A central pixel can therefore receive
a substantial contribution from a distant shell.

## Why an infinite interval can miss a shell

QuadGK maps $u\in[0,\infty)$ onto $t\in[0,1)$ using $u=t/(1-t)$. For the
center pixel, the Gaussian peak maps to $t_*=k/(1+k)$ and its approximate
width maps to $\sigma_t=\sigma/(1+k)^2$. As $k$ increases, the peak moves
toward $t=1$ and gets narrower.

For the center pixel with $k=509.5$ px and $\sigma=2.8$ px, the first 15
samples reached at most $u=233.1$ px. The Gaussian rounded to zero at every
sample. The embedded quadrature rules both returned zero, so their
difference gave no reason to subdivide. Tightening the relative tolerance
cannot recover a peak that both rules miss.
[QuadGK infinite-interval examples](https://juliamath.github.io/QuadGK.jl/stable/quadgk-examples/#Improper-integrals:-Infinite-limits)

| Integration of this entry | Result before division by $2\pi$ |
| --- | ---: |
| Finite Gaussian cutoff at $k+\sqrt{20}\,\sigma$ | 7.018531992 |
| Unsplit interval from zero to infinity | 0, with estimated error 0 |
| Infinite interval split around the Gaussian shell | 7.018559162 |

At $R=0$ the infinite Gaussian integral is
$\sigma\sqrt{2\pi}\,\Phi(k/\sigma)$, with $\Phi$ the standard normal CDF.
The finite cutoff in this example loses about 3.9 parts per million of
that integral. The unsplit call loses the entire entry.

## How the error reaches the energy spectrum

Consider an inner pixel with signal $b=g c_L+a c_H$, where $c_L$ and $c_H$
are low- and high-radius coefficients. Suppose the outer pixels determine
$c_H$ and a missed integral sets the high-radius contribution $a$ to
zero in the inner pixel's model. The fit can match the measured inner
pixel by assigning it a low-radius coefficient

$$
\widehat{c}_L = b/g = c_L + (a/g)c_H.
$$

This gives a spurious low-radius component even if the original $c_L$ is
zero. Radius maps to energy through $E=\alpha r^2$. In the full SVD solve,
small singular values can further amplify errors in the forward matrix.

## Choose integration bounds for the radial function

A radial feature at $r_i>R$ occurs at $u_i=\sqrt{(r_i-R)(r_i+R)}$.
Convert relevant support boundaries or knots to these line-of-sight
coordinates and pass them as ordered endpoints to a split quadrature.
For a Gaussian or exponential shell, split around its peak and extend
the finite limit until the spectrum is stable, or integrate the final
interval to infinity. For a Lorentzian, check the slower tail separately.
At the center pixel and for a shell well away from zero, stopping ten
widths above the shell drops about
$(\pi/2-\arctan(10))/\pi \approx 3.2\%$ of its integral.

In the package, `rBF` defines the radial profile and `zIP` sets the upper
limit of the forward integral. Supply both when you change the radial
family through custom callbacks. The built-in
[cubic B-spline](bases.md#use-a-compact-cubic-b-spline) selects its own
support bounds and splits the forward integral at its knots.
