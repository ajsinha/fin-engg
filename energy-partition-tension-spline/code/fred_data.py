#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
fred_data.py  --  retrieve real yield-curve data and run the energy-partition
                  tension-spline analysis over large historical samples.

This module is designed to run on YOUR machine (it needs outbound internet to
FRED / the U.S. Treasury). It does NOT ship any data; every curve it analyses is
downloaded live from a public source, so there is nothing fabricated here.

Two data sources are supported:

  1. FRED (https://fred.stlouisfed.org) Treasury par-yield series (DGS*).
     Requires a free API key: https://fred.stlouisfed.org/docs/api/api_key.html
     Set it via the environment variable FRED_API_KEY or pass api_key=...

  2. U.S. Treasury par-yield CSV (no key needed):
     https://home.treasury.gov/.../TextView?...  (daily Treasury par yield curve)

Typical use
-----------
    export FRED_API_KEY=xxxxxxxx
    python fred_data.py --start 2006-01-01 --end 2026-05-22 --freq W

That downloads ~20y of weekly curves, fits the equipartition tension spline to
each, and writes a CSV of per-date diagnostics (sigma*, phi distribution,
forward roughness, 5y locality, min forward) plus summary statistics by year.

Dependencies: numpy, pandas, requests, matplotlib, and tension_curve.py
(the reference implementation) on the PYTHONPATH.
"""
from __future__ import annotations
import os, sys, time, argparse, io
import numpy as np

# FRED Treasury constant-maturity series and their tenors (in years)
FRED_SERIES = {
    "DGS1MO": 1/12, "DGS3MO": 0.25, "DGS6MO": 0.5, "DGS1": 1.0, "DGS2": 2.0,
    "DGS3": 3.0, "DGS5": 5.0, "DGS7": 7.0, "DGS10": 10.0, "DGS20": 20.0, "DGS30": 30.0,
}

# ----------------------------------------------------------------------
# Data retrieval
# ----------------------------------------------------------------------
def fetch_fred_series(series_id, start, end, api_key):
    """Download one FRED series as (dates, values) using the public API."""
    import requests
    url = "https://api.stlouisfed.org/fred/series/observations"
    params = dict(series_id=series_id, api_key=api_key, file_type="json",
                  observation_start=start, observation_end=end)
    r = requests.get(url, params=params, timeout=30); r.raise_for_status()
    obs = r.json()["observations"]
    import pandas as pd
    s = pd.Series({o["date"]: (float(o["value"]) if o["value"] != "." else np.nan)
                   for o in obs})
    s.index = pd.to_datetime(s.index)
    return s

def fetch_fred_curves(start, end, api_key=None, freq="W"):
    """Assemble a DataFrame of par yields (decimals); index=date, cols=tenor (yr)."""
    import pandas as pd
    api_key = api_key or os.environ.get("FRED_API_KEY")
    if not api_key:
        raise RuntimeError("Set FRED_API_KEY env var or pass api_key=... "
                           "(free: https://fred.stlouisfed.org/docs/api/api_key.html)")
    cols = {}
    for sid, ten in FRED_SERIES.items():
        cols[ten] = fetch_fred_series(sid, start, end, api_key) / 100.0  # % -> decimal
        time.sleep(0.2)                                                  # be polite
    df = pd.DataFrame(cols).sort_index(axis=1)
    df = df.resample(freq).last().dropna(how="any")    # complete curves only
    return df

def fetch_treasury_csv(year):
    """Key-free fallback: daily Treasury par yield curve for one calendar year."""
    import requests, pandas as pd
    url = ("https://home.treasury.gov/resource-center/data-chart-center/"
           "interest-rates/daily-treasury-rates.csv/" + str(year) +
           "/all?type=daily_treasury_yield_curve&field_tdr_date_value=" + str(year))
    r = requests.get(url, timeout=30); r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text))
    return df

# ----------------------------------------------------------------------
# Per-curve analysis (delegates to the reference implementation)
# ----------------------------------------------------------------------
def analyse_curve(tenors, yields, phi_star=0.5):
    """Fit the equipartition tension spline to one curve; return diagnostics."""
    import tension_curve as tc
    x = np.asarray(tenors, float); y = np.asarray(yields, float)
    order = np.argsort(x); x, y = x[order], y[order]
    s = tc.solve_sigma_star(x, y, phi_star)
    z = tc.tension_solve(x, y, s)
    h = np.diff(x)
    # per-segment phi distribution (uniform sigma)
    phiseg = np.array([
        tc.bending_energy_closed(z[i], z[i+1], h[i], s) /
        (tc.bending_energy_closed(z[i], z[i+1], h[i], s)
         + tc.axial_energy_closed(z[i], z[i+1], y[i], y[i+1], h[i], s))
        for i in range(len(x)-1)])
    tg = np.linspace(x[0], x[-1], 1500)
    zg = tc.tension_eval(x, y, z, s, tg)
    fwd = np.gradient(tg*zg, tg)
    k5 = int(np.argmin(np.abs(x-5.0)))
    loc, _, _ = tc.locality_radius(x, y, k5, "tension", s=s)
    return dict(sigma=s, phi_seg_min=phiseg.min(), phi_seg_max=phiseg.max(),
                fwd_roughness=float(np.trapezoid(np.gradient(fwd, tg)**2, tg)),
                min_fwd=float(fwd.min()), loc5y=float(loc),
                inverted=bool(y[0] > y[-1]))

def run_historical(start, end, api_key=None, freq="W", phi_star=0.5, out_csv="tension_history.csv"):
    """Download a historical window and fit every curve; write a diagnostics CSV."""
    import pandas as pd
    df = fetch_fred_curves(start, end, api_key, freq)
    print(f"retrieved {len(df)} complete curves from {df.index.min().date()} "
          f"to {df.index.max().date()}")
    tenors = df.columns.to_numpy(float)
    recs = []
    for dt, row in df.iterrows():
        d = analyse_curve(tenors, row.to_numpy(float), phi_star)
        d["date"] = dt.date(); recs.append(d)
    res = pd.DataFrame(recs).set_index("date")
    res.to_csv(out_csv)
    print(f"wrote {out_csv}")
    print("\nsummary (sigma* and 5y locality):")
    print(res[["sigma", "loc5y", "fwd_roughness", "min_fwd"]].describe().round(4))
    print(f"\nfraction of dates with inverted curve: {res['inverted'].mean():.1%}")
    print(f"fraction with any negative implied forward: {(res['min_fwd']<0).mean():.1%}")
    # sigma* stability: year-over-year variation
    res2 = res.copy(); res2.index = pd.to_datetime(res2.index)
    by_year = res2.groupby(res2.index.year)["sigma"].agg(["mean", "std", "min", "max"])
    print("\nsigma* by year (stability check):"); print(by_year.round(3))
    return res

def plot_history(res, fname="fig_history.png"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import pandas as pd
    res = res.copy(); res.index = pd.to_datetime(res.index)
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    ax[0].plot(res.index, res["sigma"], color="#C8102E", lw=1)
    ax[0].set_ylabel(r"$\sigma^\star$ (equipartition)")
    ax[0].set_title("Path stability of the energy-partition tension over history")
    ax[1].plot(res.index, res["loc5y"], color="#1f77b4", lw=1)
    ax[1].set_ylabel("5y hedge-locality (yr)"); ax[1].set_xlabel("date")
    plt.tight_layout(); plt.savefig(fname, dpi=140)
    print(f"wrote {fname}")

# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="2006-01-01")
    ap.add_argument("--end", default="2026-05-22")
    ap.add_argument("--freq", default="W", help="pandas resample freq: D, W, M, Q")
    ap.add_argument("--phi", type=float, default=0.5, help="target energy fraction phi*")
    ap.add_argument("--api-key", default=None, help="FRED API key (or set FRED_API_KEY)")
    ap.add_argument("--out", default="tension_history.csv")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    res = run_historical(args.start, args.end, args.api_key, args.freq, args.phi, args.out)
    if args.plot:
        plot_history(res)

if __name__ == "__main__":
    main()
