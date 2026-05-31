# Energy-Partition Tension Spline — Complete Project Archive (v2.1.1)

Everything for the paper **"Swap Curve Construction and Optimization of the Tension Parameter:
An Energy-Partition Formulation"** (A. Sinha): manuscripts, reproducible code, figures, and
publication materials.

New in v2.x: a 20-year historical validation (2006--2026, 1,064 real weekly Treasury curves)
and a cross-regime benchmark comparison added to the paper; a master pipeline script
(`run_all.py`); the benchmark-history and analysis tooling; and the author's real result CSVs
(`tension_history.csv`, `benchmark_history.csv`) bundled in `code/`.

All empirical results use the **real** U.S. Treasury par-yield curve of 22 May 2026 (Federal
Reserve H.15). No data are simulated or fabricated; the illustrative SOFR/OIS levels are labelled
as such, and the bootstrap/discount-factor reproduction is exact for any inputs.

## Layout

```
energy-partition-tension-spline/
├── MANIFEST.md                  ← this file
├── paper/
│   ├── swap_curve_tension.pdf        Main paper (22 pp) — the master version
│   ├── swap_curve_tension.tex        LaTeX source of the above
│   ├── swap_curve_tension_ssrn.pdf   SSRN copy (adds Keywords + JEL on title page)
│   ├── swap_curve_tension_ssrn.tex   LaTeX source of the SSRN copy
│   └── journal_versions/
│       ├── swap_curve_tension_AMF.tex  Applied Math. Finance (T&F "interact" class)
│       └── jcf_references.tex          J. Computational Finance house-style references
├── code/
│   ├── tension_curve.py        Core module: splines, energies, phi/sigma*, locality,
│   │                           forward KKT, multi-curve, risk Jacobians, benchmarks,
│   │                           and the full validation_report().
│   ├── benchmark.py            Benchmarks (Hagan-West, Smith-Wilson, log-disc, p/w-flat),
│   │                           phi* sweep, stress cases; draws fig_benchmark.
│   ├── multicurve.py           SOFR/OIS dual-curve example; draws fig_multicurve.
│   ├── jacobian_stability.py   Analytic Jacobians + sigma* stability; draws fig_jacobian.
│   ├── fred_data.py            Run LOCALLY: fits the criterion to live FRED/Treasury data
│   │                           over arbitrary historical windows (needs internet + API key).
│   ├── benchmark_history.py    Run LOCALLY: fits all benchmark methods to every historical
│   │                           curve (resumable); writes benchmark_history.csv.
│   ├── analyze_history.py      Summarises tension_history.csv; redraws the historical figure.
│   ├── run_all.py              MASTER script: runs validation, figures, and (with a FRED key)
│   │                           both historical studies in one command.
│   ├── tension_history.csv     The author's real 1,064-week FRED run (2006-2026), bundled.
│   ├── validation_report.txt   Saved output of `python tension_curve.py --no-figures`.
│   ├── README.md               How to install, reproduce, and use on your own curve.
│   ├── requirements.txt        numpy, matplotlib (+ pandas, requests for the data scripts).
│   ├── LICENSE                 MIT (code only).
│   ├── CITATION.cff            Citation metadata (v2.0.0).
│   └── .zenodo.json            Zenodo metadata for minting a code DOI.
├── figures/                    All 14 figures (PNG) used in the paper.
└── submission/
    ├── SSRN_submission.md       Paste-ready title/abstract/keywords/JEL + highlights.
    ├── HOW_TO_PUBLISH.md        Step-by-step: SSRN preprint + Zenodo code DOI.
    ├── AMF_submission.md        Applied Mathematical Finance (Taylor & Francis) guide.
    └── JCF_submission.md        Journal of Computational Finance (Risk) guide.
```

## Quick starts

**Run everything with one command:**
```bash
cd code && pip install -r requirements.txt
python run_all.py --offline                  # validation + all figures (no data needed)
export FRED_API_KEY=xxxxxxxx                  # then, for the historical studies:
python run_all.py --full --start 2006-01-01 --end 2026-05-22 --freq W
```

**Reproduce just the paper's numbers:**
```bash
cd code && python tension_curve.py
```

**Historical study on real data (your machine):** produces `tension_history.csv` and
`benchmark_history.csv` — the two files to send back for a paper update.

**Publish:** start with `submission/HOW_TO_PUBLISH.md` (SSRN + Zenodo), then
`submission/AMF_submission.md` or `submission/JCF_submission.md` for a journal.

## Notes
- `paper/swap_curve_tension.pdf` is the definitive read; the SSRN and AMF versions are the same
  content reformatted for those venues. The JCF route requires condensing to ~5,500 words
  (see `submission/JCF_submission.md`).
- The Taylor & Francis `interact.cls` is not bundled (it is the publisher's copyrighted file);
  `AMF_submission.md` says where to get it.
- License: code under MIT (`code/LICENSE`); the paper text and figures are the author's and are
  not covered by the MIT License.

Author: Ashutosh Sinha · ajsinha@gmail.com
