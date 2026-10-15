# Upper tolerances and shadow prices in the p-median problem

Code accompanying the paper *"Upper tolerances and shadow prices in the p-median problem"*.

For a site `i` that is open in a fixed optimal solution of the p-median problem, the paper proves that the
**upper tolerance** `u(i) = z_{¬i} − z*` (the increase of the optimum when site `i` is forbidden) equals the
**left shadow price** `λ_i^-(1)` of the constraint `y_i ≤ 1` over the convex hull `conv(S)` of feasible
solutions, i.e. the largest optimal value of the corresponding dual variable.

This repository verifies the result numerically on the OR-Library instances `pmed1`–`pmed5`.

## What is computed

For every open site `i` of the optimal solution:

1. **Upper tolerance** `u(i) = z_{¬i} − z*`, where `z_{¬i}` is obtained by solving the p-median MIP with `y_i = 0`.
2. **Left shadow price** `λ_i^-(1)` by a cutting-plane method on the dual of
   `min { d : (x,y) ∈ conv(S), y_i ≤ 1 }` (Algorithm 1 in the paper). The separation problem is a
   p-median MIP with an opening fee on site `i`.
3. **Cross-check:** the finite difference `(φ_i(1−ε) − z*)/ε` with `ε = 0.01`, where `φ_i(1−ε)` is computed
   by column generation on the primal problem (Algorithm 2 in the paper).

All MIPs and LPs are solved with [HiGHS](https://highs.dev) through SciPy.

## Files

| File | Content |
|---|---|
| `pmed_tolerance.py` | Instance reader, p-median MIP oracle, Algorithm 1 (cutting planes), Algorithm 2 (column generation) |
| `run_experiments.py` | Runs the experiments and writes `results/results.json` |
| `make_tables.py` | Produces the LaTeX tables of the paper (`results/tables.tex`) from `results/results.json` |
| `results/` | Results and tables as reported in the paper |
| `data/` | Place the instance files here (see below) |

## Data

The instances `pmed1.txt`–`pmed5.txt` are from J. E. Beasley's OR-Library:

* description: <http://people.brunel.ac.uk/~mastjjb/jeb/orlib/pmedinfo.html>
* files: `http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/pmed1.txt`, …, `pmed5.txt`

Download them into `data/`. They are not redistributed here.

Note: some files contain repeated edges. As specified in the OR-Library description, the **last** occurrence of
an edge is used. With this rule the optimal values coincide with the published ones
(5819, 4093, 4250, 3034, 1355).

## Installation and usage

Requires Python ≥ 3.9.

```bash
pip install -r requirements.txt
python run_experiments.py          # all five instances, about 2–3 minutes
python make_tables.py              # writes results/tables.tex
```

Options:

```bash
python run_experiments.py --instances 1 2      # only pmed1 and pmed2
python run_experiments.py --data-dir path/to/files --out results/my_results.json
```

Expected output:

```
pmed1: n=100, p=5, z*=5819 (matches OR-Library), open sites 5, mismatches 0
pmed2: n=100, p=10, z*=4093 (matches OR-Library), open sites 10, mismatches 0
pmed3: n=100, p=10, z*=4250 (matches OR-Library), open sites 10, mismatches 0
pmed4: n=100, p=20, z*=3034 (matches OR-Library), open sites 20, mismatches 0
pmed5: n=100, p=33, z*=1355 (matches OR-Library), open sites 33, mismatches 0

open sites checked: 78, u(i) = lambda_i^-(1) in 78 of them
```

The results in the paper were obtained with Python 3.13, NumPy 2.5 and SciPy 1.18.1.

## Citation

If you use this code, please cite the paper (reference to be added after publication).

## License

MIT, see `LICENSE`.
