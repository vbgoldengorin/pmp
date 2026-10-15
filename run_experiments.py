import argparse
import json
import os
import sys
import time

import numpy as np

import pmed_tolerance as pt


def run_instance(path, name, number, eps=1e-2):
    t0 = time.time()
    c, p = pt.read_orlib(path)
    mip = pt.PMedianMIP(c, p)
    Y_star, z_star = mip.solve()
    rec = dict(name=name, n=int(c.shape[0]), p=int(p), z=z_star,
               known=pt.KNOWN_OPT.get(number),
               open=[int(i) + 1 for i in Y_star],
               duplicates=pt.duplicate_edges(path),
               M=float(c.max() * c.shape[0] + 1), sites=[])
    for i in Y_star:
        y_ub = np.ones(c.shape[0])
        y_ub[i] = 0
        _, z_not = mip.solve(y_ub=y_ub)
        lam, trace = pt.lambda_minus_cutting(mip, i, z_star)
        phi_eps, cg_it = pt.phi_column_generation(
            mip, i, 1 - eps, [(1.0, z_star), (0.0, z_not)])
        rec["sites"].append(dict(
            i=int(i) + 1, z_not=z_not, u=z_not - z_star, lam=lam,
            lam_fd=(phi_eps - z_star) / eps, cut_iterations=len(trace),
            cg_iterations=cg_it, trace=trace))
    rec["time"] = time.time() - t0
    return rec


def main():
    ap = argparse.ArgumentParser(description="Reproduce the results on pmed1..pmed5.")
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--instances", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--out", default=os.path.join("results", "results.json"))
    args = ap.parse_args()

    results, n_bad = [], 0
    for k in args.instances:
        path = os.path.join(args.data_dir, f"pmed{k}.txt")
        if not os.path.exists(path):
            sys.exit(f"missing {path} -- see README.md for where to download it")
        rec = run_instance(path, f"pmed{k}", k)
        bad = sum(abs(s["u"] - s["lam"]) > 1e-6 or abs(s["u"] - s["lam_fd"]) > 1e-4
                  for s in rec["sites"])
        n_bad += bad
        check = "" if rec["known"] is None else \
            (" (matches OR-Library)" if round(rec["z"]) == rec["known"] else " (DOES NOT match OR-Library!)")
        print(f"{rec['name']}: n={rec['n']}, p={rec['p']}, z*={rec['z']:.0f}{check}, "
              f"open sites {len(rec['sites'])}, mismatches {bad}, {rec['time']:.1f} s", flush=True)
        results.append(rec)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=1)
    total = sum(len(r["sites"]) for r in results)
    print(f"\nopen sites checked: {total}, u(i) = lambda_i^-(1) in {total - n_bad} of them")
    print(f"results written to {args.out}")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
