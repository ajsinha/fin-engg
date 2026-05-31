#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_history.py -- summarise the output of fred_data.py (tension_history.csv)
and regenerate the historical-study figure used in the paper.

Usage:
    python analyze_history.py tension_history.csv

Produces:
    fig_history.png          3-panel figure (sigma* series; locality-vs-tension; by-year boxes)
    history_summary.csv      per-year medians (sigma*, locality, inverted share, neg-fwd share)

The CSV is whatever fred_data.py wrote: columns
    date, sigma, phi_seg_min, phi_seg_max, fwd_roughness, min_fwd, loc5y, inverted
Every row is a real Treasury curve fitted by the equipartition criterion; nothing here
is simulated.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED = "#C8102E"; BLUE = "#1f77b4"


def load(path):
    df = pd.read_csv(path, parse_dates=["date"]).set_index("date").sort_index()
    return df


def summarise(df):
    out = {}
    out["n"] = len(df)
    out["start"] = df.index.min().date()
    out["end"] = df.index.max().date()
    out["sigma_median"] = df["sigma"].median()
    out["sigma_p05"] = df["sigma"].quantile(0.05)
    out["sigma_p95"] = df["sigma"].quantile(0.95)
    out["loc5y_median"] = df["loc5y"].median()
    out["neg_fwd_frac"] = (df["min_fwd"] < 0).mean()
    out["worst_fwd"] = df["min_fwd"].min()
    out["inverted_frac"] = df["inverted"].mean()
    out["wow_abs_median"] = df["sigma"].diff().abs().median()
    out["wow_abs_p95"] = df["sigma"].diff().abs().quantile(0.95)
    out["corr_log_sigma_loc"] = np.corrcoef(np.log(df["sigma"]), np.log(df["loc5y"]))[0, 1]
    out["sigma_loc_product_median"] = (df["sigma"] * df["loc5y"]).median()
    out["phi_spread_median"] = (df["phi_seg_max"] - df["phi_seg_min"]).median()
    return out


def by_year(df):
    g = df.groupby(df.index.year).agg(
        sigma_median=("sigma", "median"),
        loc5y_median=("loc5y", "median"),
        inverted=("inverted", "mean"),
        neg_fwd=("min_fwd", lambda s: (s < 0).mean()),
        n=("sigma", "size"),
    )
    return g.round(4)


def make_figure(df, fname="fig_history.png"):
    fig = plt.figure(figsize=(11, 7.2))

    ax1 = fig.add_subplot(2, 2, (1, 2))
    ax1.plot(df.index, df["sigma"], color=RED, lw=0.8)
    ax1.set_ylabel(r"$\sigma^\star$ (equipartition)")
    ax1.set_title(r"(a) Equipartition tension $\sigma^\star$, weekly U.S. Treasury curves")
    for s, e, lab in [("2008-09-01", "2009-03-01", "crisis"),
                      ("2009-03-01", "2015-12-01", "ZIRP"),
                      ("2022-03-01", "2023-08-01", "hiking"),
                      ("2023-08-01", "2025-01-01", "inversion")]:
        ax1.axvspan(pd.Timestamp(s), pd.Timestamp(e), color="0.85", alpha=0.5)
        ax1.text(pd.Timestamp(s), df["sigma"].max() * 0.92, lab, fontsize=7, color="0.4")

    ax2 = fig.add_subplot(2, 2, 3)
    sc = ax2.scatter(df["sigma"], df["loc5y"], s=6, c=df.index.year, cmap="viridis", alpha=0.6)
    c = np.corrcoef(np.log(df["sigma"]), np.log(df["loc5y"]))[0, 1]
    k = (df["sigma"] * df["loc5y"]).median()
    xx = np.linspace(df["sigma"].min(), df["sigma"].max(), 100)
    ax2.plot(xx, k / xx, color=RED, lw=1.2, ls="--", label=r"$\rho_5=c/\sigma^\star$")
    ax2.set_xlabel(r"$\sigma^\star$"); ax2.set_ylabel("5y hedge-locality radius (yr)")
    ax2.set_title(rf"(b) Locality vs. tension (corr$_{{\log}}={c:.2f}$)")
    ax2.legend(fontsize=8); cb = plt.colorbar(sc, ax=ax2); cb.set_label("year", fontsize=7)

    ax3 = fig.add_subplot(2, 2, 4)
    yrs = sorted(df.index.year.unique())
    data = [df[df.index.year == y]["sigma"].values for y in yrs]
    ax3.boxplot(data, positions=yrs, widths=0.6, showfliers=False,
                medianprops=dict(color=RED))
    ax3.set_xticks(yrs[::2]); ax3.set_xticklabels([str(y) for y in yrs[::2]], rotation=45, fontsize=7)
    ax3.set_ylabel(r"$\sigma^\star$"); ax3.set_title(r"(c) $\sigma^\star$ distribution by year")

    plt.tight_layout(); plt.savefig(fname, dpi=150)
    return fname


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "tension_history.csv"
    df = load(path)
    s = summarise(df)
    print(f"curves: {s['n']}  ({s['start']} -> {s['end']})")
    print(f"sigma*: median={s['sigma_median']:.2f}, 5-95 pct=[{s['sigma_p05']:.2f},{s['sigma_p95']:.2f}]")
    print(f"week-over-week |d sigma*|: median={s['wow_abs_median']:.3f}, 95th={s['wow_abs_p95']:.2f}")
    print(f"corr(log sigma*, log loc5y) = {s['corr_log_sigma_loc']:.3f}; median sigma*xloc = {s['sigma_loc_product_median']:.2f}")
    print(f"5y locality median = {s['loc5y_median']:.2f} yr")
    print(f"negative-forward weeks = {s['neg_fwd_frac']:.2%} (worst {s['worst_fwd']:.4%})")
    print(f"inverted (1m>30y proxy) = {s['inverted_frac']:.1%}")
    print(f"per-segment phi spread under global sigma* (median max-min) = {s['phi_spread_median']:.2f}")
    by_year(df).to_csv("history_summary.csv")
    print("wrote history_summary.csv")
    make_figure(df)
    print("wrote fig_history.png")


if __name__ == "__main__":
    main()
