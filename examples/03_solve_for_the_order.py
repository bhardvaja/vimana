"""Name the accuracy and the step count; the order you need falls out."""
from mpmath import mp, mpf, nstr
from vimana import mu, digits, order_for
mp.dps = 60

m0 = mu(mpf(2), mpf(1))
D0 = digits(m0)
print(f"start D0 = log10(3) = {nstr(D0, 8)}\n")
print(f"  {'digits wanted':>14} {'steps':>6} {'p required':>18} {'delivered':>16}")
for D, N in [(10**2, 3), (10**3, 5), (10**6, 10), (10**9, 12)]:
    p = order_for(D, N, D0)
    print(f"  {D:>14} {N:>6} {nstr(p, 12):>18} {nstr(digits(m0 ** (p ** N)), 10):>16}")
print("\np grows only like the N-th root of the target: a billion digits in twelve steps wants p ~ 6.")
