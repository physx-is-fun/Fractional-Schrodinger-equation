# Laser pulse propagation in optical fiber

This notebook shows a simple, scalar implementation of the finite difference <a href="#1">[1]</a> for solving the Nonlinear Schrödinger Equation <a href="#2">[2]</a>.

$\frac{\partial A}{\partial z}=-\frac{\alpha}{2}A+i \frac{\beta_2}{2} \frac{\partial^2 A}{\partial t^2}-i \gamma(|A|^2A)$

This nonlinear partial differential equation models how the envelope and phase of light pulse changes when propagating through a single mode optical fiber, when taking power attenuation ($\alpha$), group velocity dispersion ($\beta_2$)and waveguide nonlinearity ($\gamma$) causing self-phase modulation (SPM) into account. A is the slowly varying amplitude of the pulse envelope and t is measured in a frame of reference moving with the pulse at the group velocity $v_g$. The nonlinear Schrödinger equation (for optics) can be derived from the wave equation. However we can choose between two (Fourier transform) notations. The derivation with the negative one can be found in Ursula Keller's book <a href="#3">[3]</a>. I used this, because the scipy library uses the negative convention for the DFT. Depending on the initial width $T_0$ and the peak power $P_0$ of the incident pulse, either dispersive or non linear
effects may dominate along the fiber. It is useful to introduce two length scales, known as the dispersion length $L_D$ and the nonlinear length $L_{NL}$. Let us consider a time scale normalized to the input width $T_0$ as:

$\tau=\frac{t}{T_0}$

In addition, we introduce a normalized amplitude U as:

$A(z,\tau)=\sqrt{P_0}e^{\frac{- \alpha Z}{2}}U(z,\tau)$

We now take into consideration a space scale normalized to the fiber length as:

$\zeta=\frac{z}{L}$

Where L is the fiber length. Thus, it turns out that $U(\zeta,\tau)$ satisfies

$\frac{\partial U}{\partial \zeta}=+i \frac{L}{2 L_D} sgn(\beta_2) \frac{\partial^2 U}{\partial \tau^2}-i e^{- \alpha Z \zeta} \frac{L}{L_{NL}} (|U|^2 U)$

Where $sgn(\beta_2)=\pm 1$ depending on the sign of the coefficient $\beta_2$ and

$L_D=\frac{T_0 ^ 2}{|\beta_2|}$

$L_{NL}=\frac{1}{\gamma_0 P_0}$

The concept behind this normalization process is to exclude any kind of overflow error that may occur during solving the PDE with finite difference method. The derivation of the dimensionless transformation for the nonlinear Schrödinger equation can be found here <a href="#4">[4]</a>. The plotting functions originally came from here <a href="#5">[5]</a>.

## References

<div id="1">[1] Wikimedia Foundation. (2024, November 28). Finite difference method. Wikipedia. https://en.wikipedia.org/wiki Finite_difference_method</div>
<div id="2">[2] Wikimedia Foundation. (2024, November 21). Nonlinear Schrödinger equation. Wikipedia. https://en.wikipedia.org/wiki/Nonlinear_Schr%C3%B6dinger_equation</div>
<div id="3">[3] Keller, U. (2023). Ultrafast lasers: A comprehensive introduction to fundamental principles with practical applications. Springer International Publishing.</div>
<div id="4">[4] Felice, D. (2016, December 1). A Study of a Nonlinear Schrödinger Equation for Optical Fibers. arxiv. https://arxiv.org/pdf/1612.00358</div>
<div id="5">[5] Krarup, O. (n.d.). OLEKRARUP123/NLSE-vector-solver: Code for modelling the nonlinear Schrödinger equation for optical fiber pulse propagation. GitHub. https://github.com/OleKrarup123/NLSE-vector-solver</div>