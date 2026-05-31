# SSRN Submission Package

**Paper:** Swap Curve Construction and Optimization of the Tension Parameter: An Energy-Partition Formulation
**Author:** Ashutosh Sinha
**Contact:** ajsinha@gmail.com
**Classification (SSRN networks):** Financial Engineering; Econometrics: Mathematical Methods; Capital Markets: Asset Pricing & Valuation

---

## 1. Title (paste into SSRN "Title")

Swap Curve Construction and Optimization of the Tension Parameter: An Energy-Partition Formulation

## 2. Abstract (paste into SSRN "Abstract" — plain text, no LaTeX)

Polynomial splines have long been used for the construction of swap curves. They suffer, however, from excessive convexity, a lack of locality under perturbation of the input rates, and the possibility of negative implied forwards. Tension splines remedy these defects through a single shape parameter, but the choice of that parameter has remained largely heuristic. In this paper we place the tension spline on its proper physical footing as the equilibrium shape of a thin elastic beam held under longitudinal tension, and we derive in closed form the two energies stored by such a beam: a flexural (bending) energy and an axial (tension) energy. We show that the curvature "smoothness" norm used in the literature is the flexural energy, and that the first-derivative norm is the axial energy, so that the two are the complementary halves of one elastic energy rather than independent quantities. We then propose selecting the tension parameter by the bending-energy fraction phi(sigma) = E_bend / (E_bend + E_tens), a dimensionless quantity that runs monotonically from one (the cubic limit) to zero (the linear limit). The neutral choice phi(sigma*) = 1/2 corresponds to equipartition of elastic energy; the target phi* provides a single, transparent dial through which a desk may express its trade-off between smoothness and locality. The criterion is invariant under rescaling of the rate and time axes, removing the unit dependence of earlier energy-based schemes. On the U.S. Treasury par-yield curve of 22 May 2026 the criterion yields a unique optimum (sigma* = 4.86 at equipartition) and, relative to the natural cubic spline, reduces the hedge-locality radius of the benchmark 5-year node from 6.7 years to 0.8 year — a factor of eight — while suppressing the spurious long-end forward hump that the cubic manufactures from a flat 20-30 year segment. The framework is extended to the instantaneous forward curve under area constraints that reproduce discount factors exactly, to per-segment non-uniform tension via energy equipartition, to a multi-curve SOFR/OIS construction, and to analytic risk Jacobians whose existence makes the selected tension a smooth (non-jumping) function of the input rates. A fully reproducible reference implementation accompanies the paper, including a module that fits the criterion to live FRED / U.S. Treasury data over arbitrary historical windows.

## 3. Keywords (paste into SSRN "Keywords")

yield curve construction; term-structure interpolation; tension spline; energy-partition criterion; hedge locality; instantaneous forward curve; multi-curve SOFR/OIS; variational method; Euler-Bernoulli beam; reproducible research

## 4. JEL Classification (paste into SSRN "JEL")

C63 (Computational Techniques; Simulation Modeling); G12 (Asset Pricing; Bond Interest Rates); E43 (Interest Rates: Determination, Term Structure); C61 (Optimization Techniques); G17 (Financial Forecasting)

## 5. Highlights (optional — useful for the abstract footer or a covering note)

- Recasts the tension spline as the exact equilibrium of a tensioned elastic beam, and proves the "smoothness" and "length" penalties are the two halves of one stored elastic energy (resolving a circularity in the prior beam-energy literature).
- Introduces a dimensionless, scale-invariant selection criterion — the bending-energy fraction phi(sigma) — with a transparent default at equipartition (phi* = 1/2); proves existence/uniqueness of the optimum and characterises the dial across phi* in [0.3, 0.7].
- Demonstrates an eight-fold compression of the 5-year hedge-locality radius (6.7y to 0.8y) on the real 22-May-2026 U.S. Treasury curve, with a Green's-function (inverse-screening-length) interpretation.
- Extends the method to the instantaneous forward curve (exact discount-factor reproduction), to per-segment non-uniform tension, and to a multi-curve SOFR/OIS construction with the discounting curve held exact.
- Derives analytic risk Jacobians (O(N) curve delta) and, via the implicit function theorem, shows the selected tension is a C^1 function of the inputs — it cannot jump under market moves.
- Benchmarks even-handedly against Hagan-West monotone-convex, Smith-Wilson, log-discount and piecewise-flat forwards; ships a fully reproducible code base and a FRED/Treasury retrieval module for historical studies at scale.

## 6. Brief cover note (SSRN does not require one; include if useful, e.g. for a working-paper series)

This working paper introduces an energy-partition criterion for selecting the tension parameter
in tension-spline term-structure construction. Its contribution is not the tension spline itself
— which is well established — but a dimensionless, physically grounded, scale-invariant rule for
choosing its single parameter, together with a transparent link between that choice and hedge
locality and risk management. All empirical results use the real U.S. Treasury par-yield curve of
22 May 2026 (Federal Reserve H.15); no data are simulated or fabricated. A complete, reproducible
reference implementation accompanies the paper, and a companion module fits the criterion to live
FRED / U.S. Treasury data so the diagnostics can be reproduced over arbitrary historical windows.

---

### Submission checklist
- [ ] Upload `swap_curve_tension_ssrn.pdf` as the full text.
- [ ] Paste Title, Abstract, Keywords, JEL from sections 1-4 above.
- [ ] Author: Ashutosh Sinha; affiliation "Independent Researcher" (or as preferred); email ajsinha@gmail.com.
- [ ] (Optional) Add a link to the code repository / Zenodo DOI in the abstract footer once minted.
- [ ] Confirm you hold the rights and select the desired license on SSRN.
