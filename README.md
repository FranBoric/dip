# Diplomski rad - Matematičko modeliranje ateroskleroze

Implementacija u Fenicsu 1D modela

```math
\begin{align*}
 \dfrac{\partial M}{\partial t} &= d_1 \dfrac{\partial^2 M}{\partial x^2} + f_1(A) - \lambda_1M,\\
 \dfrac{\partial A}{\partial t} &= d_2 \dfrac{\partial^2 A}{\partial x^2} + f_2(A)M - \lambda_2A,
\end{align*}
```

i 2D modela

```math
\begin{align*}
 \dfrac{\partial M}{\partial t} &= d_1 \Delta M - \lambda_1M,\\
 \dfrac{\partial A}{\partial t} &= d_2 \Delta A + f_2(A)M - \lambda_2A,
\end{align*}
```

Gdje su 

```math
\begin{align*}
f_1(A) &= \frac{\alpha_1 + \beta_1A}{1 + A/\tau_1}, \\
f_2(A) &= \frac{\alpha_2 A}{1 + A/\tau_2}
\end{align*}
```