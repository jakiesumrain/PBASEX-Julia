# pBASEX with forward–backward asymmetry

This derivation extends the [pBASEX theory note](pbasex.pdf) to a cylindrically
symmetric distribution whose forward and backward hemispheres can differ.
It uses the same coordinate and normalization conventions. PBASEX implements
this derivation as `geometry=:half_plane` in `pbasex` and `BASISData`.
[Run a half-plane fit](usage.md#fit-forwardbackward-asymmetry)

## 1. The Abel transform

Let $x$ be the line-of-sight coordinate, and let $(y,z)$ be the detector
coordinates, with the physical symmetry axis along $z$ in the detector
plane. Write $\rho=\sqrt{x^2+y^2}$ for the distance from that axis. Axial
symmetry means the three-dimensional density can be written as
$f_{\mathrm{cyl}}(\rho,z)$. The axial coordinate $z$ retains its sign.

At each fixed $z$, the detector records

$$
\begin{aligned}
\widetilde f(y,z)
&=\int_{-\infty}^{\infty}
f_{\mathrm{cyl}}\!\left(\sqrt{x^2+y^2},z\right)\,dx \\
&=2\int_{|y|}^{\infty}
\frac{\rho\,f_{\mathrm{cyl}}(\rho,z)}{\sqrt{\rho^2-y^2}}\,d\rho.
\end{aligned}
\tag{1}
$$

This is the Abel transform of the transverse radial profile at that $z$.
For a sufficiently smooth, decaying projection, the inverse is

$$
f_{\mathrm{cyl}}(\rho,z)
=-\frac{1}{\pi}\int_\rho^\infty
\frac{\partial\widetilde f(y,z)/\partial y}{\sqrt{y^2-\rho^2}}\,dy.
\tag{2}
$$

Equations (1) and (2) act independently at each $z$. The projection has
$\widetilde f(-y,z)=\widetilde f(y,z)$, so the left and right detector halves
can be combined. The measurement at $+z$ and the one at $-z$ carry distinct
information whenever the source has forward–backward asymmetry.

## 2. Axial symmetry and Legendre orders

Use the spherical radius $r=\sqrt{\rho^2+z^2}$ and the polar angle
$\cos\theta=z/r$. As in the original pBASEX derivation, define $\sigma(r)$
as the integral of $f(r,\theta)$ over solid angle at fixed $r$. Then

$$
f(r,\theta)=\frac{\sigma(r)}{4\pi}
\left[1+\sum_{\ell=1}^{\ell_{\max}}\beta_\ell(r)P_\ell(\cos\theta)\right].
\tag{3}
$$

The truncation $\ell_{\max}$ specifies the angular model. Reflection through
the plane $z=0$ maps $\theta$ to $\pi-\theta$, and
$P_\ell(-\mu)=(-1)^\ell P_\ell(\mu)$. The even terms therefore describe the
part common to the two hemispheres. The odd terms describe their signed
difference. Choosing the opposite positive direction for the $z$ axis
reverses the signs of the odd coefficients.

## 3. The pBASEX approach

Use radial basis functions $f_k(r)$ to expand the isotropic term and the
unnormalized angular terms.

$$
\begin{aligned}
\sigma(r)&=\sum_k c_{k0} f_k(r),\\
\sigma(r)\beta_\ell(r)&=\sum_k c_{k\ell} f_k(r)
\qquad(\ell>0).
\end{aligned}
\tag{4}
$$

With $R=\sqrt{y^2+z^2}$, the projected basis function is

$$
\begin{aligned}
\widetilde f(y,z)
&=\sum_{k,\ell}c_{k\ell}\widetilde f_{k\ell}(y,z),\\
\widetilde f_{k\ell}(y,z)
&=\frac{1}{2\pi}\int_0^\infty
f_k\!\left(\sqrt{R^2+x^2}\right)
P_\ell\!\left(\frac{z}{\sqrt{R^2+x^2}}\right)\,dx .
\end{aligned}
\tag{5}
$$

The factor $1/(2\pi)$ comes from the angular normalization in (3) and the
two identical halves of the line-of-sight integral. The radial
quadrature is the same for a given shell on either side of the equatorial
plane. Its angular factor gives

$$
\widetilde f_{k\ell}(-y,z)=\widetilde f_{k\ell}(y,z),\qquad
\widetilde f_{k\ell}(y,-z)=(-1)^\ell\widetilde f_{k\ell}(y,z).
\tag{6}
$$

To see what four-quadrant folding does, consider a pair of projected
intensities at $z$ and $-z$. Their half-sum contains only even orders, and
their half-difference contains only odd orders. A fit to their sum alone
cannot recover the odd coefficients. At $z=0$, the odd columns contribute
zero.

### Discretization and fitting

On a centered detector, take $n_r$ transverse samples from $y=0$ outward
and $2n_r-1$ signed axial samples, including $z=0$ once. For a paired
left–right measurement, average the two detector values at fixed $z$.
Use the single measured value on the $y=0$ axis. This defines a half-plane
in detector-pixel intensity units. PBASEX stores rows from positive $z$ at
the top to negative $z$ at the bottom. This fixes the sign of every odd
coefficient relative to the input image orientation.

The data vector $\mathbf b$ then has $n_r(2n_r-1)$ entries. Sample
(5) at the same coordinates and with the same intensity convention to
make the rows of $G$. Each column corresponds to a pair $(k,\ell)$.
With an optional diagonal matrix of pixel multipliers $W$ and a scalar
regularization parameter $\lambda\ge0$, fit

$$
\widehat{\mathbf c}
=\arg\min_{\mathbf c}\left(
\lVert W(G\mathbf c-\mathbf b)\rVert_2^2
+\lambda\lVert\mathbf c\rVert_2^2\right).
\tag{7}
$$

The usual SVD solution applies to $WG$ and $W\mathbf b$. With paired
$+z$ and $-z$ pixels and equal multipliers within each pair, the even
and odd columns are orthogonal in the weighted detector inner product.
A mask or unequal multipliers can couple them. A joint fit on the signed
half-plane covers both cases.

The fitted projection is $G\widehat{\mathbf c}$. For the central velocity
slice, evaluate the unprojected basis
$f_k(\sqrt{y^2+z^2})P_\ell(z/\sqrt{y^2+z^2})/(4\pi)$ at $x=0$ and use
the same coefficients. Reflect both outputs across $y=0$ to reconstruct
the full detector plane while retaining both signs of $z$.

## 4. Energy spectrum and anisotropy

With the VMI calibration $E=\alpha r^2$, the volume element gives
$r^2\,dr=r\,dE/(2\alpha)$. Integrating over all solid angles leaves
only the isotropic term. The energy spectrum and anisotropy parameters are

$$
\begin{aligned}
I(E)=\frac{dN}{dE}
&=\frac{r}{2\alpha}\sum_k\widehat c_{k0}f_k(r),\\
\beta_\ell(E)
&=\frac{\sum_k\widehat c_{k\ell} f_k(r)}
        {\sum_k\widehat c_{k0} f_k(r)},\qquad
r=\sqrt{E/\alpha},\quad \ell>0 .
\end{aligned}
\tag{8}
$$

The ratio is undefined wherever the fitted isotropic term vanishes.
PBASEX also reports undefined anisotropies at the origin, where the
direction of emission is undefined.
For an energy with nonzero yield, the forward–backward asymmetry formed
from the yields $I_+(E)$ and $I_-(E)$ over their respective hemispheres is

$$
\frac{I_+(E)-I_-(E)}{I_+(E)+I_-(E)}
=\sum_{\ell\ \mathrm{odd}}\beta_\ell(E)\int_0^1 P_\ell(\mu)\,d\mu
=\frac{\beta_1(E)}{2}-\frac{\beta_3(E)}{8}
+\frac{\beta_5(E)}{16}+\cdots .
\tag{9}
$$

Thus the angle-integrated spectrum and the hemispheric asymmetry probe
different components of the same fitted distribution.

## 5. An analytic check

A nonnegative density with a known projection and a finite Legendre
expansion is

$$
f(x,y,z)=A\,e^{-(x^2+y^2+z^2)/(2s^2)}
\left(1+\eta z/s\right)^2,\qquad A>0,\quad s>0.
\tag{10}
$$

Its detector projection is

$$
\widetilde f(y,z)=A\sqrt{2\pi}\,s\,e^{-(y^2+z^2)/(2s^2)}
\left(1+\eta z/s\right)^2.
\tag{11}
$$

Let $q=\eta r/s$ and $D=1+q^2/3$. Expanding $(1+q\cos\theta)^2$ gives
$\sigma(r)=4\pi A e^{-r^2/(2s^2)}D$, $\beta_1(r)=2q/D$, and
$\beta_2(r)=2q^2/(3D)$, with all higher orders zero. Reversing $\eta$
reflects the image in $z$ and reverses $\beta_1$. It preserves $\sigma$
and $\beta_2$. Setting $\eta=0$ gives an isotropic check. Recovery with a
finite radial basis is subject to that basis's approximation of the
analytic radial profiles.
