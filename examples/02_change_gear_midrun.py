"""Speed up, slow down, pause.  Only the product of the exponents survives."""
from mpmath import mp, mpf
from vimana import mu, V, digits, run_schedule
mp.dps = 40

x, y = mpf(2), mpf(1)
print(f"start: mu = 1/3, digits = {float(digits(mu(x, y))):.4f}\n")

schedule = [mpf(3), mpf(1)/2, mpf(4), mpf(1), mpf(5)/2]
labels = ["speed up", "slow down (digits given back)", "speed up", "pause", "speed up"]
a, b = x, y
for p, lab in zip(schedule, labels):
    a, b = V(p, a, b)
    print(f"  p = {str(p):<4}  digits = {float(digits(mu(a, b))):8.4f}   {lab}")

prod = mpf(1)
for p in schedule:
    prod *= p
c, d = V(prod, x, y)
print(f"\none step at p = {prod}: digits = {float(digits(mu(c, d))):.4f}  <- identical")

import itertools
vals = {float(digits(mu(*run_schedule(list(q), x, y)))) for q in itertools.permutations(schedule)}
print(f"all {len(list(itertools.permutations(schedule)))} reorderings give {len(vals)} distinct value(s)")
