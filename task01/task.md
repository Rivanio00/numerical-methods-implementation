# Assignment 01 — Numerical Differentiation

## Description

Develop a formula for the **second derivative** using the central difference approach,
ensuring that the truncation error is of order $O(\Delta x^4)$.

The formula must be derived using two different methods:

- **(a)** Using Taylor series expansions
- **(b)** Using Newton's interpolation polynomial

---

## Application

Use the derived formula to compute the second derivative of:

$$f(x) = \sqrt{e^{3x} + 4x^2}$$

at the point $x = 2$.

---

## Convergence Study

Perform a convergence study by halving $\Delta x$ iteratively until the relative error
between successive approximations satisfies:

$$\left|\frac{f''_k - f''_{k-1}}{f''_k}\right| < 10^{-5}$$

---

## Objectives

- Derive the 5-point central difference stencil from first principles
- Confirm that both derivation methods yield the same formula
- Demonstrate $O(\Delta x^4)$ convergence through numerical experiment