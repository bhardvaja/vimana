"""The laws of V_p, exactly as stated in the paper."""
import pytest
from mpmath import mp, mpf, mpc, sqrt, cos, sin, cosh, sinh, pi, j, chebyt, chebyu
from vimana import mu, V, iterate, run_schedule, order_for, efficiency_index, root_step

mp.dps = 60
TOL = mpf(10) ** -40
X0, Y0 = mpf("1.3"), mpf("0.4")


@pytest.mark.parametrize("p", [2, 3, 5, 9, mpf("2.5"), mpf(-2), mpc(2, 1), pi])
def test_exactness(p):
    """Theorem 1: mu(V_p(x,y)) = mu(x,y)**p."""
    assert abs(mu(*V(p, X0, Y0)) - mu(X0, Y0) ** p) < TOL


@pytest.mark.parametrize("p", [2, 3, 5, 9, mpf("2.5"), mpf(-2), mpc(2, 1)])
def test_norm_law(p):
    """Theorem 3: X^2 - Y^2 = (x^2 - y^2)^p."""
    A, B = V(p, X0, Y0)
    assert abs((A * A - B * B) - (X0 * X0 - Y0 * Y0) ** p) < TOL


@pytest.mark.parametrize("p", [2, 3, 4, 5, 7, 9])
def test_identification_split_complex(p):
    """Theorem 4: V_p is the p-th power map of R[j], j^2 = +1."""
    def smul(a, b):
        return (a[0] * b[0] + a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    w = (X0, Y0)
    for _ in range(p - 1):
        w = smul(w, (X0, Y0))
    assert all(abs(u - v) < TOL for u, v in zip(w, V(p, X0, Y0)))


@pytest.mark.parametrize("p,r", [(2, 3), (3, 5), (mpf("2.5"), 4)])
def test_group_law(p, r):
    """Theorem 5: V_p o V_r = V_{pr}."""
    a, b = V(p, *V(r, X0, Y0))
    c, d = V(p * r, X0, Y0)
    assert abs(mu(a, b) - mu(c, d)) < TOL


@pytest.mark.parametrize("p", [2, 3, 5])
def test_inverse(p):
    """V_{1/p} o V_p = identity; Q^x is a group of exponents."""
    A, B = V(p, X0, Y0)
    assert all(abs(u - v) < TOL for u, v in zip(V(mpf(1) / p, A, B), (X0, Y0)))


def test_schedule_collapses_to_product():
    """Corollary 8: only the product of the exponents matters."""
    sched = [3, mpf(1) / 2, 4, 1, mpf(5) / 2]
    prod = mpf(1)
    for q in sched:
        prod *= q
    a, b = run_schedule(sched, X0, Y0)
    c, d = V(prod, X0, Y0)
    assert abs(mu(a, b) - mu(c, d)) < TOL


def test_schedule_is_order_free():
    """The group is abelian: reordering a schedule changes nothing."""
    import itertools
    sched = [3, 4, mpf("2.5"), mpf("0.7")]
    vals = [mu(*run_schedule(list(q), X0, Y0)) for q in itertools.permutations(sched)]
    assert max(abs(v - vals[0]) for v in vals) < TOL


def test_involution():
    """Proposition 10: V_{-1} is (x, -y)/sd, i.e. mu -> 1/mu."""
    s, d = X0 + Y0, X0 - Y0
    A, B = V(-1, X0, Y0)
    assert abs(A - X0 / (s * d)) < TOL and abs(B + Y0 / (s * d)) < TOL


@pytest.mark.parametrize("p", [2, 3, 4, 5, 7, 9, 12])
def test_chebyshev(p):
    """On x^2 - y^2 = 1, V_p is (T_p, y U_{p-1})."""
    u = mpf("1.6"); v = sqrt(u * u - 1)
    A, B = V(p, u, v)
    assert abs(A - chebyt(p, u)) < TOL and abs(B - v * chebyu(p - 1, u)) < TOL


@pytest.mark.parametrize("p", [2, 3, mpf("2.5"), mpf("7.3")])
def test_de_moivre(p):
    """With y -> iy on the unit circle, V_p is de Moivre."""
    t = mpf("0.7")
    A, B = V(p, cos(t), j * sin(t))
    assert abs(A - cos(p * t)) < TOL and abs(B - j * sin(p * t)) < TOL


@pytest.mark.parametrize("p", [2, 3, mpf("2.5"), pi, mpf(-3)])
def test_hyperbolic_angle_map(p):
    """V_p(cosh t, sinh t) = (cosh pt, sinh pt): the p-fold angle map."""
    t = mpf("0.83")
    A, B = V(p, cosh(t), sinh(t))
    assert abs(A - cosh(p * t)) < TOL and abs(B - sinh(p * t)) < TOL


@pytest.mark.parametrize("p", [2, 3, 4, 5, 6])
def test_brahmagupta_bhavana(p):
    """V_p over Z[sqrt D] is Brahmagupta's composition, iterated p times."""
    D, x0, y0 = 2, mpf(3), mpf(2)
    A, Bs = V(p, x0, y0 * sqrt(D)); B = Bs / sqrt(D)
    a, b = x0, y0
    for _ in range(p - 1):
        a, b = a * x0 + D * b * y0, a * y0 + b * x0
    assert abs(A * A - D * B * B - 1) < TOL
    assert abs(a - A) < TOL and abs(b - B) < TOL


def test_bhavana_off_the_lattice():
    """At p = pi the image still lies on x^2 - Dy^2 = 1, exactly."""
    D, x0, y0 = 2, mpf(3), mpf(2)
    A, Bs = V(pi, x0, y0 * sqrt(D)); B = Bs / sqrt(D)
    assert abs(A * A - D * B * B - 1) < TOL


def test_newton_and_halley():
    """root_step is Newton at p = 2 and Halley at p = 3."""
    N = mpf(7)
    assert abs(root_step(2, mpf(3), N) - (mpf(3) + N / 3) / 2) < TOL
    assert abs(root_step(3, mpf(3), N) - 3 * (9 + 3 * N) / (27 + N)) < TOL


@pytest.mark.parametrize("u,v", [(mpc(1, 2), mpc(2, -1)), (mpc(3, 0), mpc(0, 5))])
def test_julia_set_is_orthogonality_locus(u, v):
    """Theorem 19: |mu| = 1 exactly when Re(x conj y) = 0."""
    assert abs(abs(mu(u, v)) - 1) < TOL


@pytest.mark.parametrize("p", [2, 3, 9, 64, 10**6])
def test_efficiency_index_is_sqrt2(p):
    """Ostrowski's index is sqrt(2) along the whole axis."""
    assert abs(efficiency_index(p) - sqrt(mpf(2))) < TOL


def test_order_for_lands_exactly():
    """Naming D and N fixes p, and the run lands on D."""
    from vimana import digits
    m0 = mu(mpf(2), mpf(1))
    D0 = digits(m0)
    for D, N in [(100, 3), (1000, 5), (10**6, 10)]:
        p = order_for(D, N, D0)
        assert abs(digits(m0 ** (p ** N)) - D) / D < mpf(10) ** -20
