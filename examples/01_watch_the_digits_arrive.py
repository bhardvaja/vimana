"""The coordinate raised to p, over and over.  The leading zeros are the digits."""
from decimal import Decimal, getcontext
from fractions import Fraction as F
getcontext().prec = 120

def show(fr, N=30):
    d = Decimal(fr.numerator) / Decimal(fr.denominator)
    return format(d.quantize(Decimal(1).scaleb(-N)), "f")

print("pair (2,1) -> mu = 1/3 exactly, so mu_n = 3^(-p^n)\n")
for p in (2, 3, 5):
    print(f"  p = {p}")
    for n in range(5):
        e = p ** n
        s = show(F(1, 3 ** e))
        z = len(s.split(".")[1]) - len(s.split(".")[1].lstrip("0"))
        if z > 26:
            break
        print(f"    n={n}  {s}   {z} leading zeros")
    print()
