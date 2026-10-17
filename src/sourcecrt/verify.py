"""Independent certificate verifier: direct CRT sum, no decoder import.

The support family and intervals are part of the trusted problem declaration,
not assertions that a certificate may change. JSON integer arithmetic is exact.
"""
from math import prod, gcd


def require(condition):
    if not condition:
        raise ValueError("certificate check failed")


def verify(moduli, residues, bounds, supports, cert):
    try:
        n, d = len(moduli), len(bounds)
        require(n == len(residues) and d > 0 and supports)
        require(all(type(m) is int and m >= 2 for m in moduli))
        require(all(gcd(moduli[i], moduli[j]) == 1 for i in range(n) for j in range(i)))
        require(all(len(row) == d and all(type(v) is int for v in row) for row in residues))
        require(all(type(l) is int and type(u) is int and l <= u for l, u in bounds))
        canonical = sorted(set(tuple(sorted(set(s))) for s in supports))
        require(all(all(type(i) is int and 0 <= i < n for i in F) for F in canonical))
        require(cert['status'] == 'unique')
        x = cert['candidate']
        require(len(x) == d and all(type(v) is int and l <= v <= u for v, (l,u) in zip(x, bounds)))
        records = cert['records']
        require(len(records) == len(canonical))
        require(sorted(tuple(r['support']) for r in records) == canonical)
        exists = False
        for rec in records:
            F = set(rec['support'])
            good = [i for i in range(n) if i not in F]
            M = prod(moduli[i] for i in good)
            j = rec['exclude_coordinate']
            if j is not None:
                require(type(j) is int and 0 <= j < d)
                a = sum(residues[i][j]*(M//moduli[i])*pow(M//moduli[i], -1, moduli[i]) for i in good) % M
                l, u = bounds[j]
                require((l-a+M-1)//M > (u-a)//M)
            else:
                exists = True
                for j, (l,u) in enumerate(bounds):
                    require(all(x[j] % moduli[i] == residues[i][j] % moduli[i] for i in good))
                    # Since x is in the interval, its adjacent progression points
                    # being outside proves this entire box is the singleton x.
                    require(x[j]-M < l and x[j]+M > u)
        return exists
    except (AssertionError, KeyError, TypeError, ValueError, IndexError, ZeroDivisionError):
        return False
