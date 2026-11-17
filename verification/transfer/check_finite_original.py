#!/usr/bin/env python3
"""Exact-rational checks for the finite transfer-and-generation propositions.

These checks corroborate the displayed identities; they are not a proof for
arbitrary sample sizes, a statistical simulation study, or a novelty audit.
Run with Python 3 and SymPy installed.
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path

import sympy as sp

CHECKS: dict[str, int] = {}


def check(group: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(f"Failed check in {group}")
    CHECKS[group] = CHECKS.get(group, 0) + 1


def weights(a: list[sp.Expr], n: int) -> list[sp.Expr]:
    """Weight of each specified m-subset, m=0,...,n."""
    if len(a) <= n:
        raise ValueError("Need coefficients a_0,...,a_n")
    return [sp.expand(sum((-1)**j * comb(n-m, j) * a[m+j]
                          for j in range(n-m+1)))
            for m in range(n+1)]


def binomial_channel(nmax: int, s: sp.Rational) -> sp.Matrix:
    return sp.Matrix(nmax+1, nmax+1,
        lambda n, k: sp.Rational(comb(n,k))*s**k*(1-s)**(n-k)
        if k <= n else sp.S.Zero)


def verify_subsets() -> None:
    for n in range(0, 13):
        for s in [sp.Rational(0), sp.Rational(1,4), sp.Rational(2,3),
                  sp.Rational(1), sp.Rational(3,2), sp.Rational(2)]:
            a = [s**k for k in range(n+1)]
            w = weights(a, n)
            for m in range(n+1):
                check('thinning_weights', w[m] == s**m*(1-s)**(n-m))
            check('normalization', sum(comb(n,m)*w[m] for m in range(n+1)) == 1)
            gamma = sum(comb(n,m)*abs(w[m]) for m in range(n+1))
            check('exact_amplification', gamma == (abs(s)+abs(1-s))**n)
        for retained in range(n+1):
            a = [sp.Rational(comb(retained,k), comb(n,k)) if k <= retained
                 else sp.S.Zero for k in range(n+1)]
            w = weights(a, n)
            for m in range(n+1):
                expected = sp.Rational(1,comb(n,retained)) if m == retained else 0
                check('uniform_subsampling_weights', w[m] == expected)
        # A genuine mixture of two retention rates.
        a = [sp.Rational(2,5)*sp.Rational(1,3)**k +
             sp.Rational(3,5)*sp.Rational(4,5)**k for k in range(n+1)]
        w = weights(a, n)
        for m in range(n+1):
            expected = (sp.Rational(2,5)*sp.Rational(1,3)**m*sp.Rational(2,3)**(n-m)
                       +sp.Rational(3,5)*sp.Rational(4,5)**m*sp.Rational(1,5)**(n-m))
            check('mixture_representation', w[m] == expected and w[m] >= 0)
    # Arbitrary signed multipliers: expansion agrees on all subsets.
    for n in range(7):
        a = [sp.Rational((-1)**k*(k+2), k+1) for k in range(n+1)]
        a[0] = sp.S.One
        f = {mask: sp.Rational((mask*17+5)%23-11,3) for mask in range(1<<n)}
        h = {}
        for mask in range(1<<n):
            h[mask] = sum((-1)**(mask.bit_count()-sub.bit_count()) * f[sub]
                          for sub in range(1<<n) if sub & mask == sub)
        w = weights(a,n)
        lhs = sum(a[mask.bit_count()]*h[mask] for mask in range(1<<n))
        rhs = sum(w[mask.bit_count()]*f[mask] for mask in range(1<<n))
        check('arbitrary_signed_mobius_expansion', lhs == rhs)


def verify_inversions() -> None:
    for nmax in range(1,11):
        for s in [sp.Rational(1,4), sp.Rational(1,2), sp.Rational(3,4), sp.Rational(1)]:
            P = binomial_channel(nmax,s)
            S = binomial_channel(nmax,1/s)
            check('binomial_inverse', P*S == sp.eye(nmax+1))
            check('binomial_inverse_normalization', S*sp.ones(nmax+1,1) == sp.ones(nmax+1,1))
            gamma = max(sum(abs(S[i,j]) for j in range(nmax+1)) for i in range(nmax+1))
            check('binomial_inverse_amplification', gamma == (2/s-1)**nmax)
            parity = sp.Matrix([(-1)**j for j in range(nmax+1)])
            h = S*parity
            check('parity_inverse', h == sp.Matrix([(1-2/s)**j for j in range(nmax+1)]))
            second = P*sp.Matrix([value**2 for value in h])
            check('parity_variance', second-sp.ones(nmax+1,1) ==
                  sp.Matrix([(4/s-3)**j-1 for j in range(nmax+1)]))
    for upper in range(1,16):
        for k in range(1,upper+1):
            size = upper-k+1
            P = sp.Matrix(size,size,
                lambda i,j: sp.Rational(comb(k+j-1,k-1),comb(k+i,k)) if j<=i else 0)
            S = sp.zeros(size)
            for i in range(size):
                m = k+i
                S[i,i] = sp.Rational(m,k)
                if i:
                    S[i,i-1] = -sp.Rational(m-k,k)
            check('tank_inverse', P*S == sp.eye(size))
            check('tank_inverse_normalization', S*sp.ones(size,1) == sp.ones(size,1))
            gamma = max(sum(abs(S[i,j]) for j in range(size)) for i in range(size))
            check('tank_inverse_amplification', gamma == sp.Rational(2*upper,k)-1)
    eps = sp.Rational(1,5)
    P = sp.Matrix([[1-eps,eps],[eps,1-eps]])
    S = P.inv()
    check('binary_signed_inverse', P*S == sp.eye(2))
    check('binary_negative_entry', S[0,1] < 0)
    check('binary_amplification', max(sum(abs(S[i,j]) for j in range(2)) for i in range(2)) == 1/(1-2*eps))


def main() -> None:
    verify_subsets()
    verify_inversions()
    report = {
        'arithmetic': 'Exact rational and integer arithmetic using SymPy',
        'groups': CHECKS,
        'total_checks': sum(CHECKS.values()),
        'all_passed': True,
        'scope': 'Finite corroboration of identities, not universal proof or novelty verification',
    }
    path = Path(__file__).with_name('verification.json')
    path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
