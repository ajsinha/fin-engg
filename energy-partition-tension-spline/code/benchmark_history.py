#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
benchmark_history.py -- for every historical Treasury curve, fit the tension spline
AND the benchmark interpolators, and record comparison metrics per week.

This produces the per-method DISTRIBUTION of metrics across 2006-2026 and every rate
regime -- the statistical comparison reviewers ask for, as opposed to a single-curve
snapshot.

Run LOCALLY (needs internet + a free FRED API key); it is heavier than fred_data.py
(six methods x ~1000 weeks, Hagan-West integrates per grid point), so it CHECKPOINTS
to CSV after each week and is fully RESUMABLE: re-run the same command and it skips
weeks already done.

Usage
-----
    export FRED_API_KEY=xxxxxxxx
    python benchmark_history.py --start 2006-01-01 --end 2026-05-22 --freq W
    # resume after interruption: just run the same line again.

Outputs
-------
    benchmark_history.csv   one row per (date, method): fwd_roughness, min_fwd, loc5y
    (upload this back for the paper's cross-regime comparison.)

Dependencies: numpy, pandas, requests, and tension_curve.py + fred_data.py on the path.
"""
from __future__ import annotations
import os, sys, argparse, time
import numpy as np
import pandas as pd

import tension_curve as tc
from fred_data import fetch_fred_curves

METHODS = ["cubic", "tension", "hagan_west", "smith_wilson", "log_discount", "pw_flat"]


def _fwd_roughness(tg, zg):
    f = np.gradient(tg * zg, tg)
    return float(np.trapezoid(np.gradient(f, tg) ** 2, tg))

def _min_fwd(tg, zg):
    return float(np.gradient(tg * zg, tg).min())


def metrics_for_curve(x, y, ngrid=1200):
    """Return {method: (fwd_roughness, min_fwd, loc5y or nan)} for one curve."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    order = np.argsort(x); x, y = x[order], y[order]
    tg = np.linspace(x[0], x[-1], ngrid)
    k5 = int(np.argmin(np.abs(x - 5.0)))
    s = tc.solve_sigma_star(x, y, 0.5)
    out = {}

    zc = tc.cubic_eval(x, y, tc.cubic_solve(x, y), tg)
    rc, _, _ = tc.locality_radius(x, y, k5, "cubic")
    out["cubic"] = (_fwd_roughness(tg, zc), _min_fwd(tg, zc), rc)

    zt = tc.tension_eval(x, y, tc.tension_solve(x, y, s), s, tg)
    rt, _, _ = tc.locality_radius(x, y, k5, "tension", s=s)
    out["tension"] = (_fwd_roughness(tg, zt), _min_fwd(tg, zt), rt)

    zhw = tc.bench_hagan_west(x, y, tg)
    out["hagan_west"] = (_fwd_roughness(tg, zhw), _min_fwd(tg, zhw), np.nan)

    try:
        zsw = tc.bench_smith_wilson(x, y, tg)
        out["smith_wilson"] = (_fwd_roughness(tg, zsw), _min_fwd(tg, zsw), np.nan)
    except Exception:
        out["smith_wilson"] = (np.nan, np.nan, np.nan)

    zld = tc.bench_log_discount(x, y, tg)
    out["log_discount"] = (_fwd_roughness(tg, zld), _min_fwd(tg, zld), np.nan)

    zpf = tc.bench_piecewise_flat_forward(x, y, tg)
    out["pw_flat"] = (_fwd_roughness(tg, zpf), _min_fwd(tg, zpf), np.nan)

    return out, s


def load_done(path):
    if os.path.exists(path):
        prev = pd.read_csv(path)
        return prev, set(prev["date"].unique())
    return None, set()


def run(start, end, api_key, freq, out_csv, ngrid):
    df = fetch_fred_curves(start, end, api_key, freq)
    tenors = df.columns.to_numpy(float)
    print(f"retrieved {len(df)} curves {df.index.min().date()} -> {df.index.max().date()}")

    prev, done = load_done(out_csv)
    if done:
        print(f"resuming: {len(done)} dates already in {out_csv}")
    rows = []
    t0 = time.time()
    for i, (dt, row) in enumerate(df.iterrows()):
        ds = str(dt.date())
        if ds in done:
            continue
        try:
            m, s = metrics_for_curve(tenors, row.to_numpy(float), ngrid)
        except Exception as e:
            print(f"  {ds}: skipped ({e})"); continue
        for meth in METHODS:
            r, mf, loc = m[meth]
            rows.append(dict(date=ds, method=meth, sigma=s,
                             fwd_roughness=r, min_fwd=mf, loc5y=loc))
        # checkpoint every 20 curves
        if len(rows) >= 20 * len(METHODS):
            _flush(rows, out_csv); rows = []
            el = time.time() - t0
            print(f"  ...{i+1}/{len(df)} curves  ({el:.0f}s)")
    if rows:
        _flush(rows, out_csv)
    print(f"done -> {out_csv}")
    _summary(out_csv)


def _flush(rows, path):
    new = pd.DataFrame(rows)
    if os.path.exists(path):
        new.to_csv(path, mode="a", header=False, index=False)
    else:
        new.to_csv(path, index=False)


def _summary(path):
    df = pd.read_csv(path)
    print("\nper-method medians across all weeks:")
    g = df.groupby("method").agg(
        fwd_roughness_med=("fwd_roughness", "median"),
        loc5y_med=("loc5y", "median"),
        neg_fwd_frac=("min_fwd", lambda s: float((s < 0).mean())),
    )
    print(g.round(4))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="2006-01-01")
    ap.add_argument("--end", default="2026-05-22")
    ap.add_argument("--freq", default="W")
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--out", default="benchmark_history.csv")
    ap.add_argument("--ngrid", type=int, default=1200, help="forward grid resolution")
    a = ap.parse_args()
    key = a.api_key or os.environ.get("FRED_API_KEY")
    run(a.start, a.end, key, a.freq, a.out, a.ngrid)


if __name__ == "__main__":
    main()
