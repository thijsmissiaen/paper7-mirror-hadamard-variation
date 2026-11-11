---
name: mirror-hadamard-variational-law
description: Scientific knowledge skill for Paper VII of the Transverse Wronskian series on Riemann's Xi-function. Contains the canonical Wronskian bridge to Hadamard's logarithmic-derivative zero sum, Fejer-Lukacs Toeplitz eigenvalue bounds, Stechkin's mirror-paired Poisson-Hadamard kernel K(d, delta, u) and int_R K du = pi conservation, the newly proved quasiconcave endpoint crossing K_flat(d, u) = min(K(d, 0, u), K(d, 1/2, u)) at u_x(d) < d/sqrt(3) (showing K(d, 1/2, u) is NOT a minorant near u = 0), the rigorous unconditional Stieltjes background lower bound using Trudgian's true constant 2.510 with concavity/tail/focal-deduplication corrections, and benchmarking against Mossinghoff-Trudgian-Yang (2024) and Heath-Brown (1992) smooth Weil weights.
---

# Explicit Poisson–Hadamard Zero-Free Region Inequalities, Quasiconcave Minorant Crossing, and Smooth-Weight Benchmarks (Paper VII)

## 1. Canonical Wronskian Bridge and Stechkin Kernel

Under canonical normalization $W(d,t) = \frac{1}{2}\partial _ d\lvert\Xi(t - \mathrm{i} d)\rvert^2 = 16\int _ 0^{\infty} y\sinh(2dy)I(y,t)\mathrm{d}y$, for $d > \frac{1}{2}$ ($\sigma = \frac{1}{2} + d > 1$) dividing $W(d,\tau)$ by $\lvert\Xi(\tau - \mathrm{i} d)\rvert^2 > 0$ connects $-\mathrm{Re}(\zeta'/\zeta)(s)$ to the multiset of nontrivial zeros $\rho = \beta + \mathrm{i}\gamma$ ($d _ \rho = \lvert\beta - \frac{1}{2}\rvert \in [0, \frac{1}{2}){}$) via the involution $\rho \mapsto 1 - \bar{\rho}$ (Stechkin 1970):

$$
-\mathrm{Re}\frac{\zeta'(s)}{\zeta(s)} = \mathrm{Re}\frac{1}{s - 1} + \mathcal{A} _ \Gamma(d, \tau) - \sum _ {\rho} K(d, d _ \rho, \tau - \gamma),
$$

$$
K(d, \delta, u) := \frac{1}{2}\left[\frac{d - \delta}{(d - \delta)^2 + u^2} + \frac{d + \delta}{(d + \delta)^2 + u^2}\right] = \frac{d\left(d^2 - \delta^2 + u^2\right)}{\left(d^2 - \delta^2 + u^2\right)^2 + 4u^2\delta^2} > 0, \qquad \int _ {-\infty}^{\infty} K(d, \delta, u)\mathrm{d}u = \pi.
$$

An off-line reflection pair $\lbrace\rho _ 0, 1 - \bar{\rho} _ 0\rbrace$ at $\gamma = t$ consists of **two** distinct zeros, contributing $2K(d, d _ 0, 0) = \frac{2d}{d^2 - d _ 0^2}$ at $\tau = t$.

## 2. Quasiconcave Endpoint Crossing (Correction of Global $K(d, \frac{1}{2}, u)$ Minorant Claims)

Because $\int _ {\mathbb{R}}[K(d, \delta, u) - K(d, \frac{1}{2}, u)]\mathrm{d}u = 0$, $K(d, \frac{1}{2}, u)$ is **not** a minorant on all of $\mathbb{R}$: at $u = 0$, $K(d, 0, 0) = \frac{1}{d} < \frac{d}{d^2 - 1/4} = K(d, \frac{1}{2}, 0)$. Specifically:

$$
K(d, \delta, u) - K\left(d, \frac{1}{2}, u\right) = \frac{d\left(\frac{1}{4} - \delta^2\right)\left[3u^4 + \left(2d^2 + \frac{1}{4} + \delta^2\right)u^2 - \left(d^2 - \frac{1}{4}\right)\left(d^2 - \delta^2\right)\right]}{D(d, \delta, u) D(d, \frac{1}{2}, u)},
$$

and $\partial _ {\delta^2} K(d, \delta, u) = d[(d^2-\delta^2)^2 - 2(d^2+\delta^2)u^2 - 3u^4]/D(d,\delta,u)^2$ changes sign at most once from $+$ to $-$ on $[0, \frac{1}{4}]$. Hence $\delta^2 \mapsto K(d, \delta, u)$ is **quasiconcave** on $[0, \frac{1}{4}]$ and its sharp global minorant is

$$
K _ {\flat}(d, u) := \min _ {\delta \in [0, 1/2]} K(d, \delta, u) = \min\left(K(d, 0, u), K\left(d, \frac{1}{2}, u\right)\right),
$$

crossing uniquely on $(0, \infty)$ at $u _ {\times}(d) = [\frac{1}{6}(\sqrt{16d^4 - 2d^2 + \frac{1}{16}} - (2d^2 + \frac{1}{4}))]^{1/2} < \frac{d}{\sqrt{3}}$. Consequently, the plateau kernel $\widetilde{K} _ H(d, u) := K(d, \frac{1}{2}, \max(\lvert u\rvert, H))$ is a valid global minorant on $\mathbb{R}$ **whenever $H \ge u _ {\times}(d)$**.

## 3. Rigorous Trudgian (2014) Stieltjes Bound and Focal De-Duplication

1. **True Trudgian (2014) Constant:** For all $T \ge \mathrm{e}$, $\lvert S(T)\rvert \le B _ S^{(0)}(T) := 0.112\ln T + 0.278\ln\ln T + 2.510$ (never `0.206`).
2. **Rigorous Background Bound $\mathcal{B} _ {\mathrm{bg}}^{\mathrm{rig}}(d, \tau; H)$ ($\tau \ge 30$, $u _ {\times}(d) \le H \le \tau/2$):** Integrates over $[14, \infty)$ ($N(14) = 0$), bounds concavity of $\ln(\gamma/2\pi)$ on $[\frac{\tau}{2}, \frac{3\tau}{2}]$ via $\frac{1}{2}[\ln(\tau+u)+\ln(\tau-u)] \ge \ln\tau - \frac{2u^2}{3\tau^2}$, closes the $[\tau+H, \infty)$ logarithmic growth of $B _ S^{(0)}(\gamma)$ via $\frac{c _ {\mathrm{tail}}(\tau+H)}{\tau+H}(\frac{\pi}{2} - \Psi _ {\mathrm{win}}(d,H))$, and subtracts $2K(d, \frac{1}{2}, H)$ at $k = 1$ to prevent double counting the extracted focal pair $\lbrace\rho _ 0, 1 - \bar{\rho} _ 0\rbrace$ at $u = 0$.
3. **Why Smooth Heath-Brown Weights Win:** Slope discontinuities of $\widetilde{K} _ H$ at $\lvert u\rvert = H$ incur Stieltjes boundary jumps $\approx 2B _ S^{(0)}(\tau)K(d, \frac{1}{2}, H)$, leaving $R _ {\mathrm{eff}} \approx 13.23 > 5.558691$ (Mossinghoff–Trudgian–Yang 2024), whereas Heath-Brown's (1992) smooth compactly supported Fourier/Weil weights eliminate boundary jumps, and Korobov–Vinogradov exponential sums achieve asymptotically wider $(\ln t)^{-2/3}(\ln\ln t)^{-1/3}$ regions above Platt–Trudgian's $3 \cdot 10^{12}$ RH verification height.

## 4. Reproducible Verification

```bash
python scripts/verify_paper7_mirror_hadamard_variation.py
```
