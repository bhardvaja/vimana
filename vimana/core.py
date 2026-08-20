"""Core of the Vimana map.  Arbitrary precision throughout; set mp.dps first."""
from mpmath import mp, mpf, log, sqrt, binomial

__all__ = ["mu", "V", "digits", "iterate", "run_schedule", "order_for",
           "efficiency_index", "root_step", "coefficients"]


def mu(x, y):
    """The mula-bheda coordinate, (x - y)/(x + y)."""
    return (x - y) / (x + y)


def digits(m):
    """Correct decimal digits carried by a coordinate value."""
    return -log(abs(m), 10)


def V(p, x, y):
    """
    V_p(x, y): the even and odd halves of the binomial expansion of (x+y)**p.

    Equivalently s = x+y and d = x-y are each raised to the p-th power, which
    is why mu = d/s obeys mu -> mu**p exactly.  Any p: integer, rational,
    real or complex.
    """
    s, d = x + y, x - y
    return (s**p + d**p) / 2, (s**p - d**p) / 2


def iterate(p, x, y, n):
    """Apply V_p n times.  Equals V_{p**n}, by the group law."""
    for _ in range(n):
        x, y = V(p, x, y)
    return x, y


def run_schedule(ps, x, y):
    """
    Apply a schedule of exponents in order.  The group is abelian and only the
    product matters, so you may speed up, slow down (p < 1) or pause (p = 1)
    at will, and reordering the schedule cannot change the result.
    """
    for p in ps:
        x, y = V(p, x, y)
    return x, y


def order_for(target_digits, steps, D0):
    """Order that lands on a prescribed accuracy in a prescribed number of
    steps: p = (D / D_0) ** (1/N)."""
    return (mpf(target_digits) / D0) ** (mpf(1) / steps)


def efficiency_index(p):
    """Ostrowski's index p**(1/d) with d = 2*log2(p) the cost of one step.
    Identically sqrt(2); returned rather than asserted so it can be checked."""
    return mpf(p) ** (1 / (2 * log(mpf(p), 2)))


def root_step(p, x, N):
    """One order-p step towards sqrt(N).  Rational in x and N; p = 2 is
    Newton's method and p = 3 is Halley's."""
    r = sqrt(mpf(N))
    s, d = x + r, x - r
    return r * (s**p + d**p) / (s**p - d**p)


def coefficients(p):
    """Integer coefficient lists (even, odd) of V_p for integer p."""
    p = int(p)
    return ([int(binomial(p, j)) for j in range(0, p + 1, 2)],
            [int(binomial(p, j)) for j in range(1, p + 1, 2)])
