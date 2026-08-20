"""
The Bharadvaja Vimana: V_p, the p-th power map on the split-complex numbers.

In the mula-bheda coordinate mu = (x-y)/(x+y) the map is the exact monomial
mu -> mu**p, so correct digits obey D_n = p**n * D_0 with no error constant
and no restriction on p.

    >>> from mpmath import mp, mpf, sqrt
    >>> from vimana import mu, V
    >>> mp.dps = 50
    >>> x, y = sqrt(mpf(2)), mpf(1)
    >>> X, Y = V(3, x, y)
    >>> abs(mu(X, Y) - mu(x, y)**3) < mpf(10)**-45
    True

Companion code for R. Bharadvaja, "Bharadvaja Vimana: The Axis of
Convergence" (2026), doi:10.5281/zenodo.22027073.
"""
from .core import (mu, V, digits, iterate, run_schedule, order_for,
                   efficiency_index, root_step, coefficients)

__version__ = "1.0.0"
__all__ = ["mu", "V", "digits", "iterate", "run_schedule", "order_for",
           "efficiency_index", "root_step", "coefficients"]
