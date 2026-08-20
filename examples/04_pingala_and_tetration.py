"""Growing gears beat any fixed order -- and cost exactly the same."""
from mpmath import mp, mpf, log, nstr
from vimana import mu, digits
mp.dps = 4000

m0 = mpf(1) / 3
D0 = digits(m0)
print("Pingala gears 2,3,5,8,13,21,34,55 on mu = 1/3")
print("  (the matrameru of Pingala's Chandahsastra, stated as a recurrence by")
print("   Virahanka, some fourteen centuries before Fibonacci)\n")
m, P = m0, 1
lg = []
for q in [2, 3, 5, 8, 13, 21, 34, 55]:
    m, P = m ** q, P * q
    lg.append(float(log(digits(m), 10)))
    print(f"  p = {q:>2}   digits = {float(digits(m)):>18,.2f}")
d = [round(lg[i+1] - lg[i], 3) for i in range(len(lg) - 1)]
print(f"\n  first differences of log10(digits): {d}")
print("  they keep growing, so log(digits) is quadratic in n; a constant p makes it linear.\n")

logP = mpf(1 + 2 + 4 + 16 + 65536)          # gears 2, 4, 16, 65536, 2^65536
print("tetration gears 2, 4, 16, 65536, 2^65536")
print(f"  after 5 steps: total exponent 2^{int(logP)}, "
      f"just under 10^{int(log(D0,10)+logP*log(2,10))+1} digits")
print(f"  cost = 2 log2(total exponent) = {int(2*logP):,} multiplications")
print(f"  same target at p = 2 needs {int(logP/log(mpf(2),2)):,} steps, "
      f"and costs the same {int(2*logP):,}")
print("  tetration buys fewer steps, not less work.")
