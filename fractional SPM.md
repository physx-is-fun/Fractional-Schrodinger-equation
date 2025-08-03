# Laser pulse propagation in optical fiber (SPM only case)

We’re looking to solve the **fractional nonlinear Schrödinger-type equation** involving a **Caputo fractional derivative** in the spatial variable <img src="https://i.upmath.me/svg/z" alt="z" />:

<img src="https://i.upmath.me/svg/%0A%5Cfrac%7B%5Cpartial%5E%5Calpha%20A%7D%7B%5Cpartial%20z%5E%5Calpha%7D%20%3D%20-i%5Cgamma%20%7CA%7C%5E2%20A%2C%20%5Cquad%200%20%3C%20%5Calpha%20%3C%201%2C%0A" alt="
\frac{\partial^\alpha A}{\partial z^\alpha} = -i\gamma |A|^2 A, \quad 0 &lt; \alpha &lt; 1,
" />

where <img src="https://i.upmath.me/svg/%5Cfrac%7B%5Cpartial%5E%5Calpha%7D%7B%5Cpartial%20z%5E%5Calpha%7D" alt="\frac{\partial^\alpha}{\partial z^\alpha}" /> denotes the **Caputo fractional derivative** of order <img src="https://i.upmath.me/svg/%5Calpha" alt="\alpha" />.

---

### Step 1: **Understand the nature of the equation**

* The right-hand side is **local** and nonlinear.
* The left-hand side is a **Caputo fractional derivative in space**, i.e., it's **nonlocal** and history-dependent.

This is an **ODE in fractional form** if you treat <img src="https://i.upmath.me/svg/A(z)" alt="A(z)" /> as a function of one variable <img src="https://i.upmath.me/svg/z" alt="z" />, though complex and nonlinear.

---

### Step 2: **Rewrite the equation**

We write:

<img src="https://i.upmath.me/svg/%0AD_z%5E%5Calpha%20A(z)%20%3D%20-i%5Cgamma%20%7CA(z)%7C%5E2%20A(z)%2C%20%5Cquad%20A(0)%20%3D%20A_0.%0A" alt="
D_z^\alpha A(z) = -i\gamma |A(z)|^2 A(z), \quad A(0) = A_0.
" />

This is a **Caputo fractional nonlinear differential equation**. If we let:

<img src="https://i.upmath.me/svg/%0Af(A)%20%3D%20-i%5Cgamma%20%7CA%7C%5E2%20A%2C%0A" alt="
f(A) = -i\gamma |A|^2 A,
" />

then the equation becomes:

<img src="https://i.upmath.me/svg/%0AD_z%5E%5Calpha%20A(z)%20%3D%20f(A(z))%2C%20%5Cquad%20A(0)%20%3D%20A_0.%0A" alt="
D_z^\alpha A(z) = f(A(z)), \quad A(0) = A_0.
" />

---

### Step 3: **Use the Caputo fractional integral operator**

Using the **equivalent integral form** of the Caputo fractional differential equation:

<img src="https://i.upmath.me/svg/%0AA(z)%20%3D%20A_0%20%2B%20%5Cfrac%7B1%7D%7B%5CGamma(%5Calpha)%7D%20%5Cint_0%5Ez%20(z%20-%20%5Czeta)%5E%7B%5Calpha%20-%201%7D%20f(A(%5Czeta))%5C%2C%20d%5Czeta.%0A" alt="
A(z) = A_0 + \frac{1}{\Gamma(\alpha)} \int_0^z (z - \zeta)^{\alpha - 1} f(A(\zeta))\, d\zeta.
" />

Plugging in <img src="https://i.upmath.me/svg/f(A)%20%3D%20-i%5Cgamma%20%7CA%7C%5E2%20A" alt="f(A) = -i\gamma |A|^2 A" />:

<img src="https://i.upmath.me/svg/%0AA(z)%20%3D%20A_0%20-%20%5Cfrac%7Bi%5Cgamma%7D%7B%5CGamma(%5Calpha)%7D%20%5Cint_0%5Ez%20(z%20-%20%5Czeta)%5E%7B%5Calpha%20-%201%7D%20%7CA(%5Czeta)%7C%5E2%20A(%5Czeta)%5C%2C%20d%5Czeta.%0A" alt="
A(z) = A_0 - \frac{i\gamma}{\Gamma(\alpha)} \int_0^z (z - \zeta)^{\alpha - 1} |A(\zeta)|^2 A(\zeta)\, d\zeta.
" />

This is a **nonlinear Volterra integral equation of the second kind** with a weakly singular kernel.

---

### Step 4: **Special Case – Constant Modulus Assumption**

If we assume:

* <img src="https://i.upmath.me/svg/A(z)%20%3D%20R%20e%5E%7Bi%5Cphi(z)%7D" alt="A(z) = R e^{i\phi(z)}" />, with constant modulus <img src="https://i.upmath.me/svg/R%20%3D%20%7CA(z)%7C%20%3D%20%7CA_0%7C" alt="R = |A(z)| = |A_0|" />, then:

<img src="https://i.upmath.me/svg/%0A%7CA%7C%5E2%20A%20%3D%20R%5E2%20A(z).%0A" alt="
|A|^2 A = R^2 A(z).
" />

Then the equation simplifies to:

<img src="https://i.upmath.me/svg/%0AD_z%5E%5Calpha%20A(z)%20%3D%20-i%20%5Cgamma%20R%5E2%20A(z).%0A" alt="
D_z^\alpha A(z) = -i \gamma R^2 A(z).
" />

This is now **linear in <img src="https://i.upmath.me/svg/A" alt="A" />**. Let’s solve this case.

---

### Step 5: **Solve the linear fractional differential equation**

We solve:

<img src="https://i.upmath.me/svg/%0AD_z%5E%5Calpha%20A(z)%20%3D%20%5Clambda%20A(z)%2C%20%5Cquad%20A(0)%20%3D%20A_0%2C%0A" alt="
D_z^\alpha A(z) = \lambda A(z), \quad A(0) = A_0,
" />

where <img src="https://i.upmath.me/svg/%5Clambda%20%3D%20-i%5Cgamma%20R%5E2" alt="\lambda = -i\gamma R^2" />.

The solution to this Caputo equation is:

<img src="https://i.upmath.me/svg/%0AA(z)%20%3D%20A_0%20E_%5Calpha(%5Clambda%20z%5E%5Calpha)%2C%0A" alt="
A(z) = A_0 E_\alpha(\lambda z^\alpha),
" />

where <img src="https://i.upmath.me/svg/E_%5Calpha" alt="E_\alpha" /> is the **Mittag-Leffler function**:

<img src="https://i.upmath.me/svg/%0AE_%5Calpha(z)%20%3D%20%5Csum_%7Bk%3D0%7D%5E%5Cinfty%20%5Cfrac%7Bz%5Ek%7D%7B%5CGamma(%5Calpha%20k%20%2B%201)%7D.%0A" alt="
E_\alpha(z) = \sum_{k=0}^\infty \frac{z^k}{\Gamma(\alpha k + 1)}.
" />

So:

<img src="https://i.upmath.me/svg/%0AA(z)%20%3D%20A_0%20E_%5Calpha%5Cleft(-i%5Cgamma%20R%5E2%20z%5E%5Calpha%5Cright).%0A" alt="
A(z) = A_0 E_\alpha\left(-i\gamma R^2 z^\alpha\right).
" />

---

### Summary

* ✅ The **general nonlinear solution** is given by the **Volterra-type integral equation**:

<img src="https://i.upmath.me/svg/%0AA(z)%20%3D%20A_0%20-%20%5Cfrac%7Bi%5Cgamma%7D%7B%5CGamma(%5Calpha)%7D%20%5Cint_0%5Ez%20(z%20-%20%5Czeta)%5E%7B%5Calpha%20-%201%7D%20%7CA(%5Czeta)%7C%5E2%20A(%5Czeta)%5C%2C%20d%5Czeta.%0A" alt="
A(z) = A_0 - \frac{i\gamma}{\Gamma(\alpha)} \int_0^z (z - \zeta)^{\alpha - 1} |A(\zeta)|^2 A(\zeta)\, d\zeta.
" />

* ✅ Under the **constant amplitude assumption** (i.e., phase modulation only), the **explicit solution** is:

<img src="https://i.upmath.me/svg/%0AA(z)%20%3D%20A_0%20E_%5Calpha(-i%5Cgamma%20%7CA_0%7C%5E2%20z%5E%5Calpha)%0A" alt="
A(z) = A_0 E_\alpha(-i\gamma |A_0|^2 z^\alpha)
" />

---

We're now aiming to solve the **SPM-only fractional nonlinear Schrödinger equation** (fNLSE), with a **chirped Gaussian pulse** as the initial condition.

---

## 📌 Problem Statement

We want to solve the **fractional** SPM-only equation:

<img src="https://i.upmath.me/svg/%0A%5Cfrac%7B%5Cpartial%5E%5Calpha%20A(z%2C%20t)%7D%7B%5Cpartial%20z%5E%5Calpha%7D%20%3D%20-i%20%5Cgamma%20%7CA(z%2C%20t)%7C%5E2%20A(z%2C%20t)%2C%0A" alt="
\frac{\partial^\alpha A(z, t)}{\partial z^\alpha} = -i \gamma |A(z, t)|^2 A(z, t),
" />

with:

<img src="https://i.upmath.me/svg/%0AA(0%2C%20t)%20%3D%20A_0(t)%20%3D%20A_p%20%5Cexp%5Cleft(-%5Cfrac%7Bt%5E2%7D%7B2T_0%5E2%7D%20%2B%20i%20C%20t%5E2%5Cright)%0A" alt="
A(0, t) = A_0(t) = A_p \exp\left(-\frac{t^2}{2T_0^2} + i C t^2\right)
" />

<img src="https://i.upmath.me/svg/%0AT_0%20%3D%20%5Cfrac%7B%5Ctext%7BFWHM%7D%7D%7B2%20%5Csqrt%7B%5Cln%202%7D%7D%0A" alt="
T_0 = \frac{\text{FWHM}}{2 \sqrt{\ln 2}}
" />

---

## ✅ Key Insight

This equation is **local in time**, so for each time point <img src="https://i.upmath.me/svg/t_i" alt="t_i" />, we’re solving an **independent fractional ODE**:

<img src="https://i.upmath.me/svg/%0AD_z%5E%5Calpha%20A(z%3B%20t_i)%20%3D%20-i%5Cgamma%20%7CA(z%3B%20t_i)%7C%5E2%20A(z%3B%20t_i)%2C%20%5Cquad%20A(0%3B%20t_i)%20%3D%20A_0(t_i).%0A" alt="
D_z^\alpha A(z; t_i) = -i\gamma |A(z; t_i)|^2 A(z; t_i), \quad A(0; t_i) = A_0(t_i).
" />

So for each <img src="https://i.upmath.me/svg/t_i" alt="t_i" />, define:

<img src="https://i.upmath.me/svg/%0Af(A)%20%3D%20-i%5Cgamma%20%7CA%7C%5E2%20A.%0A" alt="
f(A) = -i\gamma |A|^2 A.
" />

This becomes a **Caputo fractional initial value problem (IVP)**:

<img src="https://i.upmath.me/svg/%0AD_z%5E%5Calpha%20A(z)%20%3D%20f(A)%2C%20%5Cquad%20A(0)%20%3D%20A_0.%0A" alt="
D_z^\alpha A(z) = f(A), \quad A(0) = A_0.
" />

---

Where:

* <img src="https://i.upmath.me/svg/%5Cfrac%7B%5Cpartial%5E%5Calpha%7D%7B%5Cpartial%20z%5E%5Calpha%7D" alt="\frac{\partial^\alpha}{\partial z^\alpha}" /> is the **Caputo fractional derivative** of order <img src="https://i.upmath.me/svg/%5Calpha%20%5Cin%20(0%2C%201%5D" alt="\alpha \in (0, 1]" />,
* <img src="https://i.upmath.me/svg/%5Cgamma" alt="\gamma" /> is the nonlinear coefficient,
* No dispersion term (GVD) is present.

---

## 🧠 Solution Method

For each fixed <img src="https://i.upmath.me/svg/t_i" alt="t_i" />, the solution is:

<img src="https://i.upmath.me/svg/%0AA(z%2C%20t_i)%20%3D%20A_0(t_i)%20%5Ccdot%20E_%5Calpha(-i%5Cgamma%20%7CA_0(t_i)%7C%5E2%20z%5E%5Calpha)%2C%0A" alt="
A(z, t_i) = A_0(t_i) \cdot E_\alpha(-i\gamma |A_0(t_i)|^2 z^\alpha),
" />

where <img src="https://i.upmath.me/svg/E_%5Calpha(%5Ccdot)" alt="E_\alpha(\cdot)" /> is the **Mittag-Leffler function**, a generalization of the exponential for fractional systems.

Great! You’re looking to **Fourier transform the analytical solution** of the **fractional SPM-only nonlinear Schrödinger equation** with a **chirped Gaussian pulse** as the input.

---

## 🔁 Problem Recap

You have an analytical solution of the form:

<img src="https://i.upmath.me/svg/%0AA(z%2C%20t)%20%3D%20A_0(t)%20%5Ccdot%20E_%5Calpha%5Cleft(-i%20%5Cgamma%20%7CA_0(t)%7C%5E2%20z%5E%5Calpha%5Cright)%2C%0A" alt="
A(z, t) = A_0(t) \cdot E_\alpha\left(-i \gamma |A_0(t)|^2 z^\alpha\right),
" />

where:

* <img src="https://i.upmath.me/svg/A_0(t)" alt="A_0(t)" /> is the **chirped Gaussian pulse**:

  <img src="https://i.upmath.me/svg/%0A%20%20A_0(t)%20%3D%20A_p%20%5Cexp%5Cleft(-%5Cfrac%7Bt%5E2%7D%7B2%20T_0%5E2%7D%20%2B%20i%20C%20t%5E2%5Cright)%2C%0A%20%20" alt="
  A_0(t) = A_p \exp\left(-\frac{t^2}{2 T_0^2} + i C t^2\right),
  " />
* <img src="https://i.upmath.me/svg/E_%5Calpha(%5Ccdot)" alt="E_\alpha(\cdot)" /> is the **Mittag-Leffler function**, and
* <img src="https://i.upmath.me/svg/%5Calpha%20%5Cin%20(0%2C%201)" alt="\alpha \in (0, 1)" /> is the **fractional order**.

Now we want to compute the **Fourier transform** of <img src="https://i.upmath.me/svg/A(z%2C%20t)" alt="A(z, t)" /> with respect to <img src="https://i.upmath.me/svg/t" alt="t" />, to analyze the **spectral broadening** due to fractional SPM.

---

## ⚙️ Step-by-Step Strategy

We define the Fourier transform as:

<img src="https://i.upmath.me/svg/%0A%5Chat%7BA%7D(z%2C%20%5Comega)%20%3D%20%5Cint_%7B-%5Cinfty%7D%5E%5Cinfty%20A(z%2C%20t)%5C%2C%20e%5E%7B-i%5Comega%20t%7D%5C%2C%20dt.%0A" alt="
\hat{A}(z, \omega) = \int_{-\infty}^\infty A(z, t)\, e^{-i\omega t}\, dt.
" />

Since <img src="https://i.upmath.me/svg/A(z%2C%20t)%20%3D%20A_0(t)%20%5Ccdot%20E_%5Calpha%5Cleft(-i%20%5Cgamma%20%7CA_0(t)%7C%5E2%20z%5E%5Calpha%5Cright)" alt="A(z, t) = A_0(t) \cdot E_\alpha\left(-i \gamma |A_0(t)|^2 z^\alpha\right)" />, and the **Mittag-Leffler factor is nonlinear in time through <img src="https://i.upmath.me/svg/%7CA_0(t)%7C%5E2" alt="|A_0(t)|^2" />**, the Fourier transform **cannot be computed in closed form analytically** in general.

However, we can:

### ✅ (1) Approximate or expand the Mittag-Leffler function

<img src="https://i.upmath.me/svg/%0AE_%5Calpha(-i%5Cgamma%20%7CA_0(t)%7C%5E2%20z%5E%5Calpha)%20%3D%20%5Csum_%7Bn%3D0%7D%5E%5Cinfty%20%5Cfrac%7B%5B-i%5Cgamma%20z%5E%5Calpha%5D%5En%7D%7B%5CGamma(%5Calpha%20n%20%2B%201)%7D%20%7CA_0(t)%7C%5E%7B2n%7D%0A" alt="
E_\alpha(-i\gamma |A_0(t)|^2 z^\alpha) = \sum_{n=0}^\infty \frac{[-i\gamma z^\alpha]^n}{\Gamma(\alpha n + 1)} |A_0(t)|^{2n}
" />

Then the total field becomes:

<img src="https://i.upmath.me/svg/%0AA(z%2C%20t)%20%3D%20A_0(t)%20%5Csum_%7Bn%3D0%7D%5E%5Cinfty%20%5Cfrac%7B%5B-i%5Cgamma%20z%5E%5Calpha%5D%5En%7D%7B%5CGamma(%5Calpha%20n%20%2B%201)%7D%20%7CA_0(t)%7C%5E%7B2n%7D%0A%3D%20%5Csum_%7Bn%3D0%7D%5E%5Cinfty%20%5Cfrac%7B%5B-i%5Cgamma%20z%5E%5Calpha%5D%5En%7D%7B%5CGamma(%5Calpha%20n%20%2B%201)%7D%20A_0(t)%20%7CA_0(t)%7C%5E%7B2n%7D%0A" alt="
A(z, t) = A_0(t) \sum_{n=0}^\infty \frac{[-i\gamma z^\alpha]^n}{\Gamma(\alpha n + 1)} |A_0(t)|^{2n}
= \sum_{n=0}^\infty \frac{[-i\gamma z^\alpha]^n}{\Gamma(\alpha n + 1)} A_0(t) |A_0(t)|^{2n}
" />

Now take the Fourier transform term-by-term:

<img src="https://i.upmath.me/svg/%0A%5Chat%7BA%7D(z%2C%20%5Comega)%20%3D%20%5Csum_%7Bn%3D0%7D%5E%5Cinfty%20%5Cfrac%7B%5B-i%5Cgamma%20z%5E%5Calpha%5D%5En%7D%7B%5CGamma(%5Calpha%20n%20%2B%201)%7D%20%5C%2C%20%5Cmathcal%7BF%7D%20%5Cleft%5C%7B%20A_0(t)%20%7CA_0(t)%7C%5E%7B2n%7D%20%5Cright%5C%7D%0A" alt="
\hat{A}(z, \omega) = \sum_{n=0}^\infty \frac{[-i\gamma z^\alpha]^n}{\Gamma(\alpha n + 1)} \, \mathcal{F} \left\{ A_0(t) |A_0(t)|^{2n} \right\}
" />

This reduces the problem to computing:

<img src="https://i.upmath.me/svg/%0A%5Cmathcal%7BF%7D%20%5Cleft%5C%7B%20A_0(t)%20%7CA_0(t)%7C%5E%7B2n%7D%20%5Cright%5C%7D%2C%0A" alt="
\mathcal{F} \left\{ A_0(t) |A_0(t)|^{2n} \right\},
" />

which is the **Fourier transform of a chirped Gaussian raised to an odd power**. While this doesn't have a clean closed-form, it **can be computed numerically**.

---

## 🔍 Notes on Interpretation

* The **chirp <img src="https://i.upmath.me/svg/C" alt="C" />** causes the input spectrum to broaden asymmetrically.
* The **Mittag-Leffler modulation** leads to **sub-exponential spectral broadening** compared to classical SPM.
* The **fractional order <img src="https://i.upmath.me/svg/%5Calpha" alt="\alpha" />** controls the **rate of spectral broadening** and **phase accumulation**.

  * <img src="https://i.upmath.me/svg/%5Calpha%20%5Cto%201" alt="\alpha \to 1" />: recovers classical SPM spectrum.
  * <img src="https://i.upmath.me/svg/%5Calpha%20%3C%201" alt="\alpha &lt; 1" />: weaker broadening, longer “memory”.

---

Let’s now derive an **approximate analytical expression for the spectrum** of the solution to the **fractional SPM-only nonlinear Schrödinger equation**, using the **stationary phase method**. This will give physical insight into **how the fractional nonlinearity modifies the spectral broadening**.

---

## 🧩 Problem Setup

We are analyzing the field:

<img src="https://i.upmath.me/svg/%0AA(z%2C%20t)%20%3D%20A_0(t)%20%5Ccdot%20E_%5Calpha%5Cleft(-i%20%5Cgamma%20%7CA_0(t)%7C%5E2%20z%5E%5Calpha%5Cright)%2C%0A" alt="
A(z, t) = A_0(t) \cdot E_\alpha\left(-i \gamma |A_0(t)|^2 z^\alpha\right),
" />

with:

<img src="https://i.upmath.me/svg/%0AA_0(t)%20%3D%20A_p%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt%5E2%7D%7B2%20T_0%5E2%7D%20%2B%20i%20C%20t%5E2%20%5Cright)%0A" alt="
A_0(t) = A_p \exp\left( -\frac{t^2}{2 T_0^2} + i C t^2 \right)
" />

Our goal is to approximate the **Fourier transform**:

<img src="https://i.upmath.me/svg/%0A%5Chat%7BA%7D(z%2C%20%5Comega)%20%3D%20%5Cint_%7B-%5Cinfty%7D%5E%5Cinfty%20A(z%2C%20t)%20e%5E%7B-i%5Comega%20t%7D%5C%2C%20dt%0A" alt="
\hat{A}(z, \omega) = \int_{-\infty}^\infty A(z, t) e^{-i\omega t}\, dt
" />

---

## ✏️ Step 1: Approximate the Mittag-Leffler Factor

For small or moderate <img src="https://i.upmath.me/svg/%5Cgamma%20z%5E%5Calpha" alt="\gamma z^\alpha" />, we can use the **first-order approximation** of the Mittag-Leffler function:

<img src="https://i.upmath.me/svg/%0AE_%5Calpha(-i%20%5Cgamma%20%7CA_0(t)%7C%5E2%20z%5E%5Calpha)%20%5Capprox%20%5Cexp%5Cleft(%20-i%20%5Cfrac%7B%5Cgamma%20z%5E%5Calpha%7D%7B%5CGamma(1%2B%5Calpha)%7D%20%7CA_0(t)%7C%5E2%20%5Cright)%0A" alt="
E_\alpha(-i \gamma |A_0(t)|^2 z^\alpha) \approx \exp\left( -i \frac{\gamma z^\alpha}{\Gamma(1+\alpha)} |A_0(t)|^2 \right)
" />

✅ This reduces the problem to an **exponential phase modulation**, similar to classical SPM but with a scaled nonlinear phase.

So now:

<img src="https://i.upmath.me/svg/%0AA(z%2C%20t)%20%5Capprox%20A_0(t)%20%5Ccdot%20%5Cexp%5Cleft(%20-i%20%5Cphi_%5Ctext%7BNL%7D(t)%20%5Cright)%0A" alt="
A(z, t) \approx A_0(t) \cdot \exp\left( -i \phi_\text{NL}(t) \right)
" />

with nonlinear phase shift:

<img src="https://i.upmath.me/svg/%0A%5Cphi_%5Ctext%7BNL%7D(t)%20%3D%20%5Cfrac%7B%5Cgamma%20z%5E%5Calpha%7D%7B%5CGamma(1%2B%5Calpha)%7D%20%7CA_0(t)%7C%5E2%0A" alt="
\phi_\text{NL}(t) = \frac{\gamma z^\alpha}{\Gamma(1+\alpha)} |A_0(t)|^2
" />

---

## ✏️ Step 2: Write the Full Phase

We combine the linear and nonlinear phases:

<img src="https://i.upmath.me/svg/%0AA(z%2C%20t)%20%5Capprox%20A_p%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt%5E2%7D%7B2%20T_0%5E2%7D%20%2B%20i%20C%20t%5E2%20-%20i%20%5Cphi_%5Ctext%7BNL%7D(t)%20%5Cright)%0A" alt="
A(z, t) \approx A_p \exp\left( -\frac{t^2}{2 T_0^2} + i C t^2 - i \phi_\text{NL}(t) \right)
" />

where:

<img src="https://i.upmath.me/svg/%0A%5Cphi_%5Ctext%7BNL%7D(t)%20%3D%20%5Cfrac%7B%5Cgamma%20z%5E%5Calpha%7D%7B%5CGamma(1%2B%5Calpha)%7D%20A_p%5E2%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt%5E2%7D%7BT_0%5E2%7D%20%5Cright)%0A" alt="
\phi_\text{NL}(t) = \frac{\gamma z^\alpha}{\Gamma(1+\alpha)} A_p^2 \exp\left( -\frac{t^2}{T_0^2} \right)
" />

---

## ✏️ Step 3: Stationary Phase Approximation

We now compute:

<img src="https://i.upmath.me/svg/%0A%5Chat%7BA%7D(%5Comega)%20%5Capprox%20%5Cint_%7B-%5Cinfty%7D%5E%5Cinfty%20A_p%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt%5E2%7D%7B2%20T_0%5E2%7D%20%2B%20i%20%5CPhi(t)%20%5Cright)%20dt%2C%0A%5Cquad%20%5Ctext%7Bwith%7D%20%5Cquad%0A%5CPhi(t)%20%3D%20C%20t%5E2%20-%20%5Cphi_%5Ctext%7BNL%7D(t)%20-%20%5Comega%20t%0A" alt="
\hat{A}(\omega) \approx \int_{-\infty}^\infty A_p \exp\left( -\frac{t^2}{2 T_0^2} + i \Phi(t) \right) dt,
\quad \text{with} \quad
\Phi(t) = C t^2 - \phi_\text{NL}(t) - \omega t
" />

Look for **stationary points** where <img src="https://i.upmath.me/svg/%5Cfrac%7Bd%5CPhi%7D%7Bdt%7D%20%3D%200" alt="\frac{d\Phi}{dt} = 0" />:

<img src="https://i.upmath.me/svg/%0A%5Cfrac%7Bd%5CPhi%7D%7Bdt%7D%20%3D%202%20C%20t%20%2B%20%5Comega%20-%20%5Cfrac%7Bd%5Cphi_%5Ctext%7BNL%7D%7D%7Bdt%7D%20%3D%200%0A" alt="
\frac{d\Phi}{dt} = 2 C t + \omega - \frac{d\phi_\text{NL}}{dt} = 0
" />

To compute <img src="https://i.upmath.me/svg/%5Cfrac%7Bd%5Cphi_%5Ctext%7BNL%7D%7D%7Bdt%7D" alt="\frac{d\phi_\text{NL}}{dt}" />, recall:

<img src="https://i.upmath.me/svg/%0A%5Cphi_%5Ctext%7BNL%7D(t)%20%3D%20%5Ckappa%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt%5E2%7D%7BT_0%5E2%7D%20%5Cright)%2C%20%5Cquad%20%5Ckappa%20%3A%3D%20%5Cfrac%7B%5Cgamma%20z%5E%5Calpha%20A_p%5E2%7D%7B%5CGamma(1%2B%5Calpha)%7D%0A" alt="
\phi_\text{NL}(t) = \kappa \exp\left( -\frac{t^2}{T_0^2} \right), \quad \kappa := \frac{\gamma z^\alpha A_p^2}{\Gamma(1+\alpha)}
" />

So:

<img src="https://i.upmath.me/svg/%0A%5Cfrac%7Bd%5Cphi_%5Ctext%7BNL%7D%7D%7Bdt%7D%20%3D%20%5Ckappa%20%5Ccdot%20%5Cleft(-%5Cfrac%7B2t%7D%7BT_0%5E2%7D%20%5Cright)%20%5Ccdot%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt%5E2%7D%7BT_0%5E2%7D%20%5Cright)%0A" alt="
\frac{d\phi_\text{NL}}{dt} = \kappa \cdot \left(-\frac{2t}{T_0^2} \right) \cdot \exp\left( -\frac{t^2}{T_0^2} \right)
" />

Thus, the stationary point condition becomes:

<img src="https://i.upmath.me/svg/%0A2%20C%20t%20%2B%20%5Comega%20%2B%20%5Cfrac%7B2%20%5Ckappa%20t%7D%7BT_0%5E2%7D%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt%5E2%7D%7BT_0%5E2%7D%20%5Cright)%20%3D%200%0A" alt="
2 C t + \omega + \frac{2 \kappa t}{T_0^2} \exp\left( -\frac{t^2}{T_0^2} \right) = 0
" />

This equation implicitly defines <img src="https://i.upmath.me/svg/t%20%3D%20t_%5Comega" alt="t = t_\omega" />, the **dominant time contributing to frequency <img src="https://i.upmath.me/svg/%5Comega" alt="\omega" />**.

---

## 🔍 Step 4: Invert to Get Spectral Phase

The key idea in stationary phase is that the main contribution to <img src="https://i.upmath.me/svg/%5Chat%7BA%7D(%5Comega)" alt="\hat{A}(\omega)" /> comes from <img src="https://i.upmath.me/svg/t%20%3D%20t_%5Comega" alt="t = t_\omega" /> satisfying the above equation.

Then the spectral phase is approximately:

<img src="https://i.upmath.me/svg/%0A%5Carg%5Chat%7BA%7D(%5Comega)%20%5Capprox%20%5CPhi(t_%5Comega)%0A" alt="
\arg\hat{A}(\omega) \approx \Phi(t_\omega)
" />

And the **spectral amplitude**:

<img src="https://i.upmath.me/svg/%0A%7C%5Chat%7BA%7D(%5Comega)%7C%20%5Capprox%20A_p%20%5Csqrt%7B2%5Cpi%7D%20%5Cleft%7C%20%5Cfrac%7Bd%5E2%5CPhi%7D%7Bdt%5E2%7D(t_%5Comega)%20%5Cright%7C%5E%7B-1%2F2%7D%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt_%5Comega%5E2%7D%7B2%20T_0%5E2%7D%20%5Cright)%0A" alt="
|\hat{A}(\omega)| \approx A_p \sqrt{2\pi} \left| \frac{d^2\Phi}{dt^2}(t_\omega) \right|^{-1/2} \exp\left( -\frac{t_\omega^2}{2 T_0^2} \right)
" />

---

## ✅ Approximate Analytical Spectrum Summary

The approximate spectral intensity is:

<img src="https://i.upmath.me/svg/%0A%7C%5Chat%7BA%7D(z%2C%20%5Comega)%7C%5E2%20%5Capprox%20%5Cfrac%7B2%5Cpi%20A_p%5E2%7D%7B%5Cleft%7C%20%5Cfrac%7Bd%5E2%5CPhi%7D%7Bdt%5E2%7D(t_%5Comega)%20%5Cright%7C%7D%20%5Cexp%5Cleft(%20-%5Cfrac%7Bt_%5Comega%5E2%7D%7BT_0%5E2%7D%20%5Cright)%0A" alt="
|\hat{A}(z, \omega)|^2 \approx \frac{2\pi A_p^2}{\left| \frac{d^2\Phi}{dt^2}(t_\omega) \right|} \exp\left( -\frac{t_\omega^2}{T_0^2} \right)
" />

Where <img src="https://i.upmath.me/svg/t_%5Comega" alt="t_\omega" /> solves:

<img src="https://i.upmath.me/svg/%0A2%20C%20t_%5Comega%20%2B%20%5Comega%20%2B%20%5Cfrac%7B2%20%5Ckappa%20t_%5Comega%7D%7BT_0%5E2%7D%20e%5E%7B-t_%5Comega%5E2%2FT_0%5E2%7D%20%3D%200%0A" alt="
2 C t_\omega + \omega + \frac{2 \kappa t_\omega}{T_0^2} e^{-t_\omega^2/T_0^2} = 0
" />

This gives an **implicit mapping between frequency <img src="https://i.upmath.me/svg/%5Comega" alt="\omega" />** and the corresponding **dominant time <img src="https://i.upmath.me/svg/t_%5Comega" alt="t_\omega" />**, including the fractional memory effect through <img src="https://i.upmath.me/svg/%5Ckappa%20%5Cpropto%20z%5E%5Calpha" alt="\kappa \propto z^\alpha" />.

---

## 🔬 Interpretation

* In classical SPM (<img src="https://i.upmath.me/svg/%5Calpha%20%3D%201" alt="\alpha = 1" />), this recovers the usual time-frequency mapping via nonlinear phase shift <img src="https://i.upmath.me/svg/%5Cphi_%7B%5Ctext%7BNL%7D%7D(t)%20%5Cpropto%20%7CA_0(t)%7C%5E2" alt="\phi_{\text{NL}}(t) \propto |A_0(t)|^2" />.
* For <img src="https://i.upmath.me/svg/%5Calpha%20%3C%201" alt="\alpha &lt; 1" />, the **nonlinear frequency shift is reduced** — the SPM-induced chirp grows sublinearly with <img src="https://i.upmath.me/svg/z" alt="z" />.
* As a result, **spectral broadening is suppressed** compared to classical SPM.
* The spectrum may **saturate** more quickly due to the Mittag-Leffler’s sub-exponential behavior.

---

## 📘 Optional: Small-<img src="https://i.upmath.me/svg/z" alt="z" /> Expansion

If <img src="https://i.upmath.me/svg/z" alt="z" /> is small (weak nonlinearity), we can approximate:

<img src="https://i.upmath.me/svg/%0AE_%5Calpha(-i%20%5Cgamma%20%7CA_0(t)%7C%5E2%20z%5E%5Calpha)%20%5Capprox%201%20-%20i%20%5Cfrac%7B%5Cgamma%20z%5E%5Calpha%7D%7B%5CGamma(1%2B%5Calpha)%7D%20%7CA_0(t)%7C%5E2%0A" alt="
E_\alpha(-i \gamma |A_0(t)|^2 z^\alpha) \approx 1 - i \frac{\gamma z^\alpha}{\Gamma(1+\alpha)} |A_0(t)|^2
" />

So:

<img src="https://i.upmath.me/svg/%0AA(z%2C%20t)%20%5Capprox%20A_0(t)%20%5Cleft%5B1%20-%20i%20%5Cfrac%7B%5Cgamma%20z%5E%5Calpha%7D%7B%5CGamma(1%2B%5Calpha)%7D%20%7CA_0(t)%7C%5E2%20%5Cright%5D%0A" alt="
A(z, t) \approx A_0(t) \left[1 - i \frac{\gamma z^\alpha}{\Gamma(1+\alpha)} |A_0(t)|^2 \right]
" />

Then:

<img src="https://i.upmath.me/svg/%0A%5Chat%7BA%7D(z%2C%20%5Comega)%20%5Capprox%20%5Chat%7BA%7D_0(%5Comega)%20-%20i%20%5Cfrac%7B%5Cgamma%20z%5E%5Calpha%7D%7B%5CGamma(1%2B%5Calpha)%7D%20%5Ccdot%20%5Cmathcal%7BF%7D%5Cleft%5B%20%7CA_0(t)%7C%5E2%20A_0(t)%20%5Cright%5D(%5Comega)%0A" alt="
\hat{A}(z, \omega) \approx \hat{A}_0(\omega) - i \frac{\gamma z^\alpha}{\Gamma(1+\alpha)} \cdot \mathcal{F}\left[ |A_0(t)|^2 A_0(t) \right](\omega)
" />

This gives a **first-order correction to the spectrum**, useful for benchmarking numerical results.

---