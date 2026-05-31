# Energy-Partition Tension Spline for Swap-Curve Construction

Reference implementation accompanying the paper

> **Swap Curve Construction and Optimization of the Tension Parameter: An Energy-Partition Formulation**, A. Sinha.

The paper places the tension spline on a variational footing (the equilibrium of a thin elastic
beam under longitudinal tension) and selects the tension parameter by the dimensionless
**bending-energy fraction** `phi(sigma) = E_bend / (E_bend + E_tens)`, with the neutral default at
equipartition `phi* = 1/2`. This repository reproduces every quantitative result in the paper and
provides tools to apply the method to your own (or live historical) curves.

All empirical results use the **real** U.S. Treasury par-yield curve of 22 May 2026 (Federal
Reserve H.15). Nothing in this repository is simulated or fabricated; the multi-curve SOFR levels
used for illustration are clearly labelled as representative, and the bootstrap/discount-factor
reproduction is exact regardless of the levels.

## Contents

| File | Purpose |
|------|---------|
| `tension_curve.py` | **Core module.** Cubic, uniform- and non-uniform-tension splines; closed-form + quadrature energies; `phi(sigma)` and `sigma*`; hedge locality; forward-space KKT; multi-curve SOFR/OIS; risk Jacobians; benchmark methods; full `validation_report()`. |
| `benchmark.py` | Head-to-head comparison vs. Hagan-West monotone-convex, Smith-Wilson, log-discount, piecewise-flat forward; `phi*` sweep; stress cases. Draws the benchmark figure. |
| `multicurve.py` | SOFR/OIS dual-curve example: crisp OIS discount bootstrap + tensioned projection forward. |
| `jacobian_stability.py` | Analytic curve-delta and `dsigma*/dy` (implicit function theorem); `sigma*` stress-stability. |
| `fred_data.py` | **Run locally.** Retrieves real Treasury curves from FRED (or the U.S. Treasury CSV) and fits the criterion over arbitrary historical windows. Needs internet + a free FRED API key. |
| `analyze_history.py` | Summarises a `tension_history.csv` and regenerates the historical-study figure + per-year table. |
| `benchmark_history.py` | **Run locally.** Fits all benchmark methods to every historical curve; writes per-method metrics over 2006-2026 (resumable, checkpointed). |
| `run_all.py` | **Master script.** Runs the whole pipeline: validation, figures, and (with a FRED key) both historical studies. |

The core module is self-contained; the others import from it.

## Installation

```bash
python -m venv venv && source venv/bin/activate    # optional
pip install -r requirements.txt
```

Requires Python 3.9+ and numpy / matplotlib (plus pandas + requests for `fred_data.py`).

## Reproducing the paper

```bash
python tension_curve.py            # prints the full validation report + writes figures
python tension_curve.py --no-figures
```

### One command for everything (`run_all.py`)

```bash
python run_all.py --offline                 # steps 1-3: validation + all figures (no data needed)

export FRED_API_KEY=xxxxxxxx                 # free key: https://fred.stlouisfed.org/docs/api/api_key.html
python run_all.py --full --start 2006-01-01 --end 2026-05-22 --freq W
```

The `--full` run additionally produces `tension_history.csv` and `benchmark_history.csv` (the two
files to send back for the paper's historical sections). `benchmark_history.py` is resumable: if a
long run is interrupted, re-run the same command and it skips weeks already done.

The validation report recomputes, from scratch, every headline number in the paper, e.g.:

```
[V/B] closed-form vs quadrature: E_bend rel.err=4e-15, E_tens rel.err=1e-12
[V  ] phi(0+)=1.0000, phi(inf)=0.0116; strictly decreasing on grid: True
[V  ] sigma*: phi*=0.75 -> 2.76; phi*=0.50 -> 4.86; phi*=0.25 -> 8.69
[V  ] phi invariance: decimals=0.500000, bps=0.500000, time-rescaled=0.500000
[VI ] hedge-locality radius at 5y: cubic=6.66 yr, tension=0.79 yr (8.4x tighter)
[X  ] SOFR/OIS: OIS bootstrap reprices swaps to 6e-17; ...
[XI ] analytic curve-delta df/dy vs finite-diff: max err=2e-11; sigma* smooth under stress
[XII] benchmark forward roughness + phi* sweep
```

The companion scripts reproduce the corresponding figures:

```bash
python benchmark.py          # Table/figure: benchmarks, phi* sweep, stress cases
python multicurve.py         # Figure: SOFR/OIS dual curve
python jacobian_stability.py # Figure: analytic Jacobians + sigma* stability
```

## Using it on your own curve

```python
import numpy as np
from tension_curve import solve_sigma_star, tension_solve, tension_eval

tenors = np.array([0.25, 1, 2, 5, 10, 30.0])
yields = np.array([0.040, 0.041, 0.042, 0.044, 0.046, 0.047])

sigma = solve_sigma_star(tenors, yields, phi_target=0.5)   # equipartition tension
z     = tension_solve(tenors, yields, sigma)
grid  = np.linspace(tenors[0], tenors[-1], 500)
curve = tension_eval(tenors, yields, z, sigma, grid)        # interpolated curve
```

## Historical study on live data (`fred_data.py`)

This runs on **your** machine (it needs outbound internet). Get a free FRED API key at
<https://fred.stlouisfed.org/docs/api/api_key.html>, then:

```bash
export FRED_API_KEY=xxxxxxxx
python fred_data.py --start 2006-01-01 --end 2026-05-22 --freq W --plot
```

It downloads complete Treasury curves over the window, fits the equipartition tension spline to
each, and writes `tension_history.csv` (per-date sigma*, phi distribution, forward roughness, 5y
locality, min forward) plus year-by-year sigma* stability statistics and a history plot.

## Data provenance

- Base curve: U.S. Treasury CMT par yields, 22 May 2026, Federal Reserve H.15 release.
- Historical curves (`fred_data.py`): downloaded live from FRED (series `DGS1MO`...`DGS30`) or the
  U.S. Treasury daily par-yield CSV. No curve is bundled or synthesised.
- SOFR/OIS multi-curve levels: representative of the period and clearly labelled as illustrative;
  the bootstrap and discount-factor reproduction are exact for any input levels.

## Citing

See `CITATION.cff`. If you mint a Zenodo release, cite the archived DOI for the exact code version.

## License

Code is released under the MIT License (`LICENSE`). The accompanying paper text/figures are
the author's and are not covered by the MIT License.
