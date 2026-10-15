import argparse
import json


def fmt(x):
    x = 0.0 if abs(x) < 1e-9 else x
    return f"{x:.0f}" if abs(x - round(x)) < 1e-9 else f"{x:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/results.json")
    ap.add_argument("--out", default="results/tables.tex")
    args = ap.parse_args()
    R = json.load(open(args.results))

    L = []
    tot = sum(len(r["sites"]) for r in R)
    ok = sum(abs(s["u"] - s["lam"]) < 1e-6 and abs(s["u"] - s["lam_fd"]) < 1e-4
             for r in R for s in r["sites"])
    L.append(f"\\newcommand{{\\TotSites}}{{{tot}}}\\newcommand{{\\TotOK}}{{{ok}}}")

    L.append(r"\newcommand{\DetailTable}{%")
    L.append(r"\begin{longtable}{@{}lrrrrrc@{}}")
    L.append(r"\caption{Per-site results. $z_{\neg i}$: optimum with site $i$ closed; "
             r"$u(i)=z_{\neg i}-z^*$; $\lambda_i^-(1)$: maximal optimal $\sigma_i$ obtained by "
             r"\cref{alg:cut}; $\hat\lambda_i$: finite-difference estimate~\eqref{eq:fd}; "
             r"it.: iterations of \cref{alg:cut}.}\label{tab:detail}\\")
    hdr = (r"\toprule Instance & site $i$ & $z_{\neg i}$ & $u(i)$ & $\lambda_i^-(1)$ & "
           r"$\hat\lambda_i$ & it.\\\midrule")
    L.append(hdr + r"\endfirsthead")
    L.append(r"\multicolumn{7}{@{}l}{\emph{(continued)}}\\" + hdr + r"\endhead\bottomrule\endfoot")
    for r in R:
        for k, s in enumerate(r["sites"]):
            name = "\\texttt{" + r["name"] + "}" if k == 0 else ""
            L.append(f"{name} & {s['i']} & {fmt(s['z_not'])} & {fmt(s['u'])} & "
                     f"{fmt(s['lam'])} & {fmt(s['lam_fd'])} & {s['cut_iterations']}\\\\")
        if r is not R[-1]:
            L.append(r"\midrule")
    L.append(r"\end{longtable}}")

    L.append(r"\newcommand{\DupTable}{%")
    L.append(r"\begin{tabular}{@{}lll@{}}\toprule Instance & edge & "
             r"lengths in order of appearance (last one used)\\\midrule")
    for r in R:
        for a, b, w1, w2 in r["duplicates"]:
            L.append(f"\\texttt{{{r['name']}}} & $\\{{{a},{b}\\}}$ & {w1} $\\to$ \\textbf{{{w2}}}\\\\")
    L.append(r"\bottomrule\end{tabular}}")

    r = R[0]
    L.append(r"\newcommand{\TraceTable}{%")
    L.append(r"\begin{tabular}{@{}rrrrrrl@{}}\toprule site $i$ & it. & $\sigma$ & "
             r"$v=z^*+\sigma$ & $d(Y)$ & $y_i(Y)$ & separation outcome\\\midrule")
    for s in r["sites"]:
        for k, t in enumerate(s["trace"], 1):
            out = "violated; cut added" if t["violated"] else "no violated cut; stop"
            L.append(f"{s['i'] if k == 1 else ''} & {k} & {fmt(t['sigma'])} & {fmt(t['v'])} & "
                     f"{fmt(t['dY'])} & {int(t['a'])} & {out}\\\\")
        if s is not r["sites"][-1]:
            L.append(r"\midrule")
    L.append(r"\bottomrule\end{tabular}}")
    L.append(f"\\newcommand{{\\BigM}}{{{fmt(r['M'])}}}")

    with open(args.out, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"written {args.out}: {tot} open sites, {ok} with u(i) = lambda_i^-(1)")


if __name__ == "__main__":
    main()
