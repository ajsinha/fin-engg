#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py -- one command to run the entire pipeline for the
energy-partition tension-spline project.

It runs, in order:
  1. validation_report           reproduce every number in the paper (no data needed)
  2. paper figures               regenerate the in-paper figures (no data needed)
  3. companion figures           benchmark / multicurve / jacobian figures (no data needed)
  4. historical retrieval        fred_data.py  -> tension_history.csv      (NEEDS FRED key)
  5. historical analysis         analyze_history.py -> fig_history.png + summary
  6. benchmark history           benchmark_history.py -> benchmark_history.csv (NEEDS key; SLOW)

Steps 1-3 are offline and always run. Steps 4-6 need internet and a FRED API key; they
are skipped automatically if no key is found, unless you pass --offline to force skip or
--full to require them.

Usage
-----
    # offline only (steps 1-3):
    python run_all.py --offline

    # everything (set your free key first):
    export FRED_API_KEY=xxxxxxxx
    python run_all.py --full --start 2006-01-01 --end 2026-05-22 --freq W

What to upload back
-------------------
    After a --full run, send back:
        tension_history.csv
        benchmark_history.csv
    (these two CSVs are all that is needed to update the paper's historical sections.)

Requires: numpy, matplotlib, pandas, requests; all project modules in the same folder.
"""
import os, sys, subprocess, argparse, time

HERE = os.path.dirname(os.path.abspath(__file__))


def banner(msg):
    print("\n" + "=" * 70 + f"\n  {msg}\n" + "=" * 70)


def run_py(args):
    """Run a python subprocess in this folder, streaming output; return exit code."""
    cmd = [sys.executable] + args
    print("$ " + " ".join(cmd))
    return subprocess.call(cmd, cwd=HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="2006-01-01")
    ap.add_argument("--end", default="2026-05-22")
    ap.add_argument("--freq", default="W")
    ap.add_argument("--offline", action="store_true", help="run only offline steps 1-3")
    ap.add_argument("--full", action="store_true", help="require the online steps 4-6")
    ap.add_argument("--skip-figures", action="store_true")
    a = ap.parse_args()

    key = os.environ.get("FRED_API_KEY")
    do_online = (a.full or (key and not a.offline))
    if a.full and not key:
        print("ERROR: --full requested but FRED_API_KEY is not set."); sys.exit(2)

    t0 = time.time()

    banner("STEP 1/6  Validation report (offline)")
    run_py(["tension_curve.py", "--no-figures"])

    if not a.skip_figures:
        banner("STEP 2/6  Paper figures (offline)")
        run_py(["tension_curve.py"])           # writes the core figures
        banner("STEP 3/6  Companion figures (offline)")
        for script in ("benchmark.py", "multicurve.py", "jacobian_stability.py"):
            if os.path.exists(os.path.join(HERE, script)):
                run_py([script])
    else:
        banner("STEPS 2-3  Figures skipped (--skip-figures)")

    if not do_online:
        banner("STEPS 4-6  Skipped (no FRED_API_KEY; offline mode)")
        print("Set FRED_API_KEY and re-run with --full to produce the historical CSVs.")
        print(f"\nTotal time: {time.time()-t0:.0f}s")
        return

    common = ["--start", a.start, "--end", a.end, "--freq", a.freq]

    banner("STEP 4/6  Historical retrieval -> tension_history.csv  (online)")
    rc = run_py(["fred_data.py"] + common + ["--plot"])
    if rc != 0:
        print("retrieval failed; stopping before analysis."); sys.exit(rc)

    banner("STEP 5/6  Historical analysis -> fig_history.png + history_summary.csv")
    run_py(["analyze_history.py", "tension_history.csv"])

    banner("STEP 6/6  Benchmark history -> benchmark_history.csv  (online, SLOW, resumable)")
    run_py(["benchmark_history.py"] + common)

    banner("DONE")
    print("Upload these back for the paper update:")
    print("   - tension_history.csv")
    print("   - benchmark_history.csv")
    print(f"\nTotal time: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
