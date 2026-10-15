import numpy as np
from scipy.optimize import milp, linprog, LinearConstraint, Bounds
from scipy.sparse import coo_matrix, csr_matrix


KNOWN_OPT = {1: 5819, 2: 4093, 3: 4250, 4: 3034, 5: 1355, 6: 7824, 7: 5631,
             8: 4445, 9: 2734, 10: 1255, 11: 7696, 12: 6634, 13: 4374,
             14: 2968, 15: 1729, 16: 8162, 17: 6999, 18: 4809, 19: 2845,
             20: 1789}

TOL = 1e-6


def read_orlib(path):
    with open(path) as f:
        tok = f.read().split()
    n, ne, p = int(tok[0]), int(tok[1]), int(tok[2])
    c = np.full((n, n), np.inf)
    np.fill_diagonal(c, 0)
    edges = np.array(tok[3:3 + 3 * ne], dtype=float).reshape(ne, 3)
    for a, b, w in edges:
        a, b = int(a) - 1, int(b) - 1
        c[a, b] = c[b, a] = w
    for k in range(n):
        c = np.minimum(c, c[:, [k]] + c[[k], :])
    return c, p


def duplicate_edges(path):
    with open(path) as f:
        tok = f.read().split()
    ne = int(tok[1])
    edges = np.array(tok[3:3 + 3 * ne], dtype=int).reshape(ne, 3)
    seen, dup = {}, []
    for a, b, w in edges:
        key = (min(a, b), max(a, b))
        if key in seen:
            dup.append((int(key[0]), int(key[1]), int(seen[key]), int(w)))
        seen[key] = w
    return dup


class PMedianMIP:

    def __init__(self, c, p):
        m, n = c.shape
        self.c, self.m, self.n, self.p = c, m, n, p
        nx = m * n
        idx = np.arange(nx)
        I, J = idx // n, idx % n
        A1 = coo_matrix((np.ones(nx), (J, idx)), shape=(n, nx + m))
        A2 = coo_matrix((np.r_[np.ones(nx), -np.ones(nx)],
                         (np.r_[idx, idx], np.r_[idx, nx + I])),
                        shape=(nx, nx + m))
        A3 = coo_matrix((np.ones(m), (np.zeros(m), nx + np.arange(m))),
                        shape=(1, nx + m))
        self.constraints = [LinearConstraint(csr_matrix(A1), 1, 1),
                            LinearConstraint(csr_matrix(A2), -np.inf, 0),
                            LinearConstraint(csr_matrix(A3), p, p)]
        self.integrality = np.r_[np.zeros(nx), np.ones(m)]

    def solve(self, f=None, y_ub=None):
        m, n = self.m, self.n
        f = np.zeros(m) if f is None else f
        ub = np.r_[np.ones(m * n), np.ones(m) if y_ub is None else y_ub]
        res = milp(np.r_[self.c.ravel(), f], constraints=self.constraints,
                   integrality=self.integrality, bounds=Bounds(0, ub),
                   options={"mip_rel_gap": 0, "disp": False})
        if res.status != 0:
            raise RuntimeError(res.message)
        Y = np.flatnonzero(np.rint(res.x[m * n:]).astype(int))
        return Y, float(self.c[Y].min(axis=0).sum())


def lambda_minus_cutting(mip, i, z_star, max_iter=1000):
    M = float(mip.c.max() * mip.n + 1)
    cuts = [(1.0, z_star)]
    trace = []
    for _ in range(max_iter):
        r = linprog([0, -1],
                    A_ub=[[1.0, -a] for a, _ in cuts],
                    b_ub=[d for _, d in cuts],
                    A_eq=[[1, -1]], b_eq=[z_star],
                    bounds=[(None, None), (0, M)], method="highs")
        v, sigma = r.x
        f = np.zeros(mip.m)
        f[i] = sigma
        Y, dY = mip.solve(f=f)
        a = float(i in Y)
        violated = dY + sigma * a < v - TOL
        trace.append(dict(sigma=float(sigma), v=float(v), dY=dY, a=a,
                          violated=bool(violated)))
        if not violated:
            return float(sigma), trace
        cuts.append((a, dY))
    raise RuntimeError("cutting-plane method did not converge")


def phi_column_generation(mip, i, t, init_columns, max_iter=1000):
    cols = list(init_columns)
    for it in range(1, max_iter + 1):
        a = np.array([col[0] for col in cols])
        d = np.array([col[1] for col in cols])
        r = linprog(d, A_ub=a[None, :], b_ub=[t],
                    A_eq=np.ones((1, len(d))), b_eq=[1],
                    bounds=(0, None), method="highs")
        v, rho = r.eqlin.marginals[0], r.ineqlin.marginals[0]
        f = np.zeros(mip.m)
        f[i] = -rho
        Y, dY = mip.solve(f=f)
        aY = float(i in Y)
        if dY - rho * aY - v >= -TOL:
            return float(r.fun), it
        cols.append((aY, dY))
    raise RuntimeError("column generation did not converge")
