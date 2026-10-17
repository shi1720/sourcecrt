"""Exact bounded-vector CRT with a finite family of admissible fault supports.

No numerical tolerances, statistical error assumptions, or unverified proposals.
Moduli must be pairwise coprime; primality is not needed.
"""
from itertools import combinations
from math import gcd, prod


def validate(moduli, residues, bounds, supports):
    if not bounds or any(type(l) is not int or type(u) is not int or l > u for l, u in bounds):
        raise ValueError('nonempty integer intervals required')
    n = len(moduli)
    if len(residues) != n or any(len(r) != len(bounds) for r in residues):
        raise ValueError('residue dimensions')
    if any(type(v) is not int for row in residues for v in row):
        raise ValueError('integer residues required')
    if any(type(m) is not int or m < 2 for m in moduli):
        raise ValueError('invalid modulus')
    if any(gcd(moduli[i], moduli[j]) != 1 for i in range(n) for j in range(i)):
        raise ValueError('moduli must be pairwise coprime')
    if not supports or any(any(type(i) is not int or i < 0 or i >= n for i in s) for s in supports):
        raise ValueError('invalid support family')


def crt(rs, ms):
    """Incremental CRT; empty system has residue zero modulo one."""
    a, M = 0, 1
    for r, m in zip(rs, ms):
        a += M * (((r-a)*pow(M, -1, m)) % m)
        M *= m
    return a, M


def progression(a, M, lower, upper):
    lo = -((a-lower)//M)
    hi = (upper-a)//M
    return lo, hi


def grouped_supports(groups, f):
    """Maximal supports suffice for at-most-f group failures (also overlap)."""
    if not 0 <= f <= len(groups):
        raise ValueError('fault budget')
    return sorted(set(tuple(sorted(set().union(*(groups[i] for i in c))))
                      for c in combinations(range(len(groups)), f)))


def analyze(moduli, residues, bounds, supports):
    """Return exact Cartesian-progression boxes, without enumerating points."""
    validate(moduli, residues, bounds, supports)
    boxes = []
    for support in supports:
        bad = set(support)
        good = [i for i in range(len(moduli)) if i not in bad]
        ms = [moduli[i] for i in good]
        rows = []
        for j, (l, u) in enumerate(bounds):
            a, M = crt([residues[i][j] for i in good], ms)
            lo, hi = progression(a, M, l, u)
            rows.append((a, M, lo, hi))
        boxes.append({'support': list(support), 'rows': rows})
    return boxes


def classify(boxes):
    candidate = None
    for box in boxes:
        rows = box['rows']
        if any(lo > hi for a, M, lo, hi in rows):
            continue
        if any(lo < hi for a, M, lo, hi in rows):
            return 'ambiguous', None
        point = tuple(a+M*lo for a, M, lo, hi in rows)
        if candidate is not None and candidate != point:
            return 'ambiguous', None
        candidate = point
    return ('inconsistent', None) if candidate is None else ('unique', candidate)


def certificate(moduli, residues, bounds, supports):
    boxes = analyze(moduli, residues, bounds, supports)
    status, x = classify(boxes)
    if status != 'unique':
        return {'status': status}
    records = []
    for box in boxes:
        excluded = next((j for j, (_, _, lo, hi) in enumerate(box['rows']) if lo > hi), None)
        records.append({'support': box['support'], 'exclude_coordinate': excluded})
    return {'status': status, 'candidate': list(x), 'records': records}


def separation_modulus(moduli, supports):
    """Exact universal injectivity threshold, for an interval box."""
    return min(prod(m for i, m in enumerate(moduli) if i not in set(F) | set(G))
               for F in supports for G in supports)


def scalar_status(moduli, residues, bounds, supports):
    return [classify(analyze(moduli, [[r[j]] for r in residues], [b], supports))[0]
            for j, b in enumerate(bounds)]


def greedy_layout(moduli, h):
    """Classical least-load placement, with load log(product of moduli)."""
    if h < 1:
        raise ValueError('positive source count required')
    loads = [1]*h
    groups = [set() for _ in range(h)]
    for i,m in enumerate(moduli):
        g=min(range(h),key=lambda j:loads[j])
        groups[g].add(i)
        loads[g]*=m
    return groups
