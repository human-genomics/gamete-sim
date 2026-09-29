"""Command line: python -m gamete_sim {build-map,simulate,report}."""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import genetic_map, pedsim, plots, qc
from . import simulate as sim
from .constants import SEXES


def cmd_build_map(args) -> None:
    summary = genetic_map.build(args.bherer, args.chain, args.out_dir)
    total = summary.drop(columns="chrom").sum()
    kept = total["markers_hg38"] / total["markers_hg19"]
    print(f"  markers: {int(total['markers_hg19']):,} hg19 -> {int(total['markers_hg38']):,} hg38 ({kept:.3%} kept)")
    print(f"  dropped: unmapped {int(total['unmapped'])}, other chromosome {int(total['other_chrom'])}, "
          f"duplicate {int(total['duplicate'])}, out of order {int(total['out_of_order'])}")
    for sex in SEXES:
        print(f"  {sex} map length: {total[f'{sex}_cM_hg19']:.2f} cM hg19 -> {total[f'{sex}_cM_hg38']:.2f} cM hg38")


def cmd_simulate(args) -> None:
    gmap = genetic_map.read_map(args.map)
    intf = sim.read_interference(args.interference)
    t0 = time.time()
    table = sim.simulate(gmap, intf, args.sex, args.n, args.seed)
    sim.write(table, args.out)
    xo = len(table) / args.n - 22
    print(f"  {args.sex}: {args.n:,} gametes, {len(table):,} segments, "
          f"{xo:.2f} crossovers per gamete, {time.time() - t0:.1f} s -> {args.out}")


def cmd_report(args) -> int:
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    gmap = genetic_map.read_map(args.map)

    summaries, totals, breaks = {}, {}, {}
    for sex in SEXES:
        path = Path(args.gametes_dir) / f"gametes_{sex}_hg38.csv"
        table = pd.read_csv(path)
        qc.validate_gametes(table, args.n)
        print(f"  [ok] {path.name}: {args.n:,} gametes, all segments tile chr1-22 (1-based, inclusive)")
        summaries[sex] = qc.summarize(table, gmap, sex)
        totals[sex] = qc.crossover_counts(table).sum(axis=1)
        breaks[sex] = pedsim.breakpoints_from_segments(table)
    summary = pd.concat(summaries.values(), ignore_index=True)
    summary.to_csv(out / "crossover_summary.tsv", sep="\t", index=False, float_format="%.4f")
    for sex in SEXES:
        g = summaries[sex].iloc[-1]
        print(f"  {sex}: {g['mean_crossovers']:.2f} crossovers per gamete (map length {g['map_length_M']:.2f} M), "
              f"variance {g['var_crossovers']:.2f}")

    comparison, pedsim_totals, pedsim_breaks = None, None, None
    pedsim_bin = args.pedsim_bin or shutil.which("ped-sim")
    if pedsim_bin:
        comparison, pedsim_totals, pedsim_breaks = {}, {}, {}
        work = out / "pedsim"
        work.mkdir(exist_ok=True)
        simmap = work / "genetic_map_hg38.simmap"
        genetic_map.write_simmap(gmap, simmap)
        for sex in SEXES:
            bp_file = pedsim.run(pedsim_bin, simmap, Path(args.interference), sex, args.n, args.seed, work)
            theirs = pedsim.read_breakpoints(bp_file, sex, args.n)
            pedsim_breaks[sex] = theirs
            comparison[sex] = pedsim.compare(breaks[sex], theirs, gmap, sex, args.n)
            pedsim_totals[sex] = np.bincount(theirs["gamete"] - 1, minlength=args.n)
            print(f"  ped-sim {sex}: {pedsim_totals[sex].mean():.2f} crossovers per gamete, "
                  f"variance {pedsim_totals[sex].var(ddof=1):.2f}")
        table = pd.concat(comparison.values(), ignore_index=True)
        se = np.sqrt((table["gamete_sim_var"] + table["pedsim_var"]) / args.n)
        table["mean_diff_z"] = (table["gamete_sim_mean"] - table["pedsim_mean"]) / se
        table.to_csv(out / "pedsim_concordance.tsv", sep="\t", index=False, float_format="%.4f")
        print(f"  ped-sim concordance ({pedsim_bin}): max |z| of per-chromosome mean difference "
              f"{table['mean_diff_z'].abs().max():.2f} over {len(table)} tests; "
              f"min KS p (positions) {table['position_ks_p'].min():.3f}, "
              f"(spacing) {table['spacing_ks_p'].min():.3f}")
    else:
        print("  [skip] ped-sim not found: concordance check not run (Docker image includes it)")

    plots.crossover_figure(summaries, totals, comparison, pedsim_totals, out / "crossovers.png")
    plots.interference_figure(breaks, pedsim_breaks, gmap, args.n, out / "interference.png")
    plots.landscape_figure(breaks, gmap, args.n, out / "landscape_chr1.png")
    print(f"  figures -> {out}/")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="gamete-sim", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("build-map", help="lift the Bherer et al. map to hg38")
    p.add_argument("--bherer", required=True)
    p.add_argument("--chain", required=True)
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=cmd_build_map)

    p = sub.add_parser("simulate", help="simulate gametes")
    p.add_argument("--map", required=True)
    p.add_argument("--interference", required=True)
    p.add_argument("--sex", choices=SEXES, required=True)
    p.add_argument("--n", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_simulate)

    p = sub.add_parser("report", help="validate outputs, compare with ped-sim, draw figures")
    p.add_argument("--map", required=True)
    p.add_argument("--interference", required=True)
    p.add_argument("--gametes-dir", required=True)
    p.add_argument("--n", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--pedsim-bin", default=None)
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
