#!/usr/bin/env python3
"""
Regenerate every number printed in the paper and check it against the page.

Each block recomputes a figure or table from scratch and compares with the
value as typeset.  Nothing here trusts the paper; the paper is the thing being
tested.  Table 11 (wall-clock cost) is machine dependent and is reported but
not asserted -- see Remark 13.

    python reproduce/reproduce_all.py
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F
from mpmath import (mp, mpf, mpc, sqrt, log, pi, qfrom, kfrom, findroot, nstr)

sys.path.insert(0, ".")
from vimana import mu, V, digits, run_schedule, order_for, efficiency_index, coefficients

getcontext().prec = 200
FAILS = []


def check(label, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {label}{('   ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(label)


def head(t):
    print(f"\n{t}\n{'-' * len(t)}")


# ---------------------------------------------------------------- Figure 1
head("Figure 1  the coordinate raised to p, over and over")
def plain(fr, N=30):
    d = Decimal(fr.numerator) / Decimal(fr.denominator)
    return format(d.quantize(Decimal(1).scaleb(-N)), "f")

fig1 = {(2, 2): "0.012345679012345679012345679012",
        (2, 3): "0.000152415790275872580399329371",
        (2, 4): "0.000000023230573125418774637910",
        (3, 3): "0.000000000000131137265239709251",
        (5, 2): "0.000000000001180235387157383257"}
check("pair (2,1), mu = 1/3, rows at p = 2, 3, 5",
      all(plain(F(1, 3 ** (p ** n))) == v for (p, n), v in fig1.items()))
check("pair (10,9), mu = 1/19, row n = 3",
      plain(F(1, 19) ** (2 ** 3)) == "0.000000000058880459747221542992")

# ---------------------------------------------------------------- Table 1
head("Table 1  changing gear in the middle of a run")
mp.dps = 60
m0 = mpf(1) / 3
D0 = digits(m0)
sched = [mpf(3), mpf(1) / 2, mpf(4), mpf(1), mpf(5) / 2]
want = ["1.4314", "0.7157", "2.8627", "2.8627", "7.1568"]
m, got = m0, []
for q in sched:
    m = m ** q
    got.append(digits(m))
check("five gear steps", all(abs(g - mpf(w)) < mpf("5e-5") for g, w in zip(got, want)),
      " ".join(nstr(g, 6) for g in got))
prod = mpf(1)
for q in sched:
    prod *= q
check("one step at p = 15 gives the same", abs(digits(m0 ** prod) - mpf("7.1568")) < mpf("5e-5"))

# ---------------------------------------------------------------- Table 2
head("Table 2  naming the outcome and solving for the order")
for D, N, ps, W in [(10**2, 3, "5.94005425917", 100), (10**3, 5, "4.61608978261", 1000),
                    (10**6, 10, "4.28683851151", 10**6), (10**9, 12, "5.98109824762", 10**9)]:
    p = order_for(D, N, D0)
    check(f"D = 1e{len(str(D))-1}, N = {N}",
          nstr(p, 12).startswith(ps[:11]) and abs(digits(m0 ** (p ** N)) - W) / W < mpf("1e-20"),
          f"p = {nstr(p, 12)}")

# ---------------------------------------------------------------- Table 3
head("Table 3  tetration in the schedule, and what it costs")
logP = mpf(1 + 2 + 4 + 16 + 65536)
check("total exponent is 2^65559", int(logP) == 65559)
check("cost 2 log2(total) = 131,118", int(2 * logP) == 131118)
check("steps at p = 2 is 65,559", int(logP / log(mpf(2), 2)) == 65559)
check("steps at p = 5 is 28,234", int(logP / log(mpf(5), 2)) == 28234)
check("digit count just under 10^19735", 19734 < float(log(D0, 10) + logP * log(2, 10)) < 19735)

# ---------------------------------------------------------------- Table 5
head("Table 5  coefficients of V_p")
check("p = 9 row", coefficients(9) == ([1, 36, 126, 84, 9], [9, 84, 126, 36, 1]))
check("p = 5 row", coefficients(5) == ([1, 10, 5], [5, 10, 1]))

# ---------------------------------------------------------------- Sec 1.1 + Table 6
head("Section 1.1 and Table 6  the order is a real number")
mp.dps = 4000
e0 = 3 - 2 * sqrt(mpf(2))
E0 = digits(e0)
pirow = [2.405, 7.556, 23.74, 74.57, 234.3, 736.0]
m, ok_law, ok_page = e0, True, True
for n in range(1, 7):
    m = m ** pi
    ok_law &= abs(digits(m) - E0 * pi ** n) < mpf("1e-500")
    ok_page &= abs(float(digits(m)) - pirow[n - 1]) / pirow[n - 1] < 2e-4
check("p = pi is exactly pi^n D_0 at 4000 digits", ok_law)
check("matches the six values as printed", ok_page)

# ---------------------------------------------------------------- Table 8
head("Table 8  which chart makes the AGM a squaring")
mp.dps = 60
kp = mpf("0.3"); k = sqrt(1 - kp ** 2)
eps = (1 - sqrt(kp)) / (1 + sqrt(kp)) / 2
Q = qfrom(k=k)
k1 = (1 - kp) / (1 + kp); kp1 = sqrt(1 - k1 ** 2)
eps1 = (1 - sqrt(kp1)) / (1 + sqrt(kp1)) / 2
check("the nome squares exactly", abs(qfrom(k=k1) - Q ** 2) < mpf("1e-40"))
check("eps = mu/2 does not (3.91e-5)", abs(float(eps1 - eps ** 2) - 3.91e-5) < 1e-7)
check("Q - eps = 1.34e-4", abs(float(Q - eps) - 1.34e-4) < 1e-6)

# ---------------------------------------------------------------- Table 9
head("Table 9  order 3/2 for the elliptic machinery")
mp.dps = 50
def ascend(kk):
    return sqrt(1 - ((1 - kk) / (1 + kk)) ** 2)
def degree3(kk):
    f = lambda L: sqrt(kk * L) + sqrt(sqrt(1 - kk ** 2) * sqrt(1 - L ** 2)) - 1
    return findroot(f, kfrom(q=qfrom(k=kk) ** 3))
kk = mpf(1) / sqrt(2)
rows, want9 = [], ["2.04656", "3.06985", "4.60477", "6.90716"]
for _ in range(4):
    kk = ascend(degree3(kk))
    rows.append(digits(qfrom(k=kk)))
check("four steps, every operation algebraic",
      all(abs(g - mpf(w)) < mpf("1e-5") for g, w in zip(rows, want9)),
      " ".join(nstr(r, 6) for r in rows))

# ---------------------------------------------------------------- Table 10
head("Table 10  complex exponents oscillate")
mp.dps = 60
q0 = mpf(1) / 10
D = digits(q0)
for P, w in [(mpc(2, 1), [2.0, 3.0, 2.0, -7.0, -38.0]),
             (mpc(mpf(1) / 2, 2), [0.5, -3.75, -5.875, 10.06, 35.03])]:
    g = [float((P ** n).real * D) for n in range(1, 6)]
    check(f"p = {P}", all(abs(a - b) < 0.01 for a, b in zip(g, w)))

# ---------------------------------------------------------------- efficiency
head("Remark 13  Ostrowski's efficiency index")
check("p^(1/2log2 p) = sqrt(2) for every p",
      all(abs(efficiency_index(q) - sqrt(mpf(2))) < mpf("1e-40")
          for q in [2, 3, 5, 9, 64, 1024, 10 ** 6]))

# ---------------------------------------------------------------- summary
print()
if FAILS:
    print(f"{len(FAILS)} MISMATCH(ES): " + "; ".join(FAILS))
    sys.exit(1)
print("Every number printed in the paper reproduces.")
