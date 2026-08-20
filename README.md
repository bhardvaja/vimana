# Bhāradvāja Vimāna

**The Axis of Convergence**

[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22027073.svg)](https://doi.org/10.5281/zenodo.22027073)

**For every `p`, integer or not, an elementary map multiplies the number of correct digits by
exactly `p` at every step.**

Send the pair `(x, y)` to the even and the odd half of the binomial expansion of `(x+y)^p`. In
the mūla-bheda coordinate `μ = (x−y)/(x+y)` that map is nothing but

```
μ  ↦  μ^p
```

so the digit count obeys `D_n = p^n · D_0`. That is an equality, not an estimate. There is no
error constant, no asymptotic regime, and no restriction on `p`. The map also composes:
`V_p ∘ V_r = V_pr`. So the order of convergence stops being a label you work out after the
method is built, and becomes **a coordinate on a group**.

This is the companion code for the paper *Bhāradvāja Vimāna: The Axis of Convergence*
([doi:10.5281/zenodo.22027073](https://doi.org/10.5281/zenodo.22027073)). It continues the
mūla-bheda coordinate of
[*Bhāradvāja Mūla-Bheda Jyāmiti*](https://doi.org/10.5281/zenodo.20044554).

## Install

```bash
pip install -e ".[dev]"
```

## Sixty seconds

```python
from mpmath import mp, mpf, sqrt, pi
from vimana import mu, V, digits, run_schedule, order_for

mp.dps = 60
x, y = sqrt(mpf(2)), mpf(1)

X, Y = V(pi, x, y)                              # p need not be an integer
abs(mu(X, Y) - mu(x, y)**pi) < mpf(10)**-58     # exact, to working precision

run_schedule([3, mpf(1)/2, 4, 1, mpf(5)/2], x, y)  # change gear mid-run;
                                                   # only the product matters
order_for(1000, 5, digits(mu(x, y)))            # 1000 digits in 5 steps? p = 4.1995...
```

## What is here

| | |
|---|---|
| `vimana/` | the map, the coordinate, schedules, and the design helpers |
| `tests/` | the laws of `V_p` as pytest: exactness, norm, group, inverses, Chebyshev, de Moivre, Brahmagupta, Julia set |
| `reproduce/` | **regenerates every number printed in the paper and checks it against the page** |
| `examples/` | four short scripts: watching digits arrive, changing gear, solving for the order, Piṅgala and tetration gears |

Run the checks:

```bash
pytest tests -q
```

and the reproduction:

```bash
python reproduce/reproduce_all.py
```

It recomputes Figure 1, Tables 1, 2, 3, 5, 6, 8, 9, 10 and the efficiency index from scratch,
then compares each against the typeset value. Nothing here trusts the paper. The paper is the
thing being tested.

## Some things it will tell you

- **`p` is a real number.** At `p = π` the first six steps carry `2.405, 7.556, 23.74, 74.57,
  234.3, 736.0` digits, each exactly `π^n·D_0`. A family indexed by its order has no member
  at `p = π`.
- **Order is not merit.** Reaching `D` digits costs `2 log₂ D` multiplications at any integer
  order, and Ostrowski's efficiency index is identically `√2` along the whole axis. The quartic
  algorithm *is* the quadratic composed with itself.
- **Gears.** Speed up, slow down (`p < 1` gives digits back), pause (`p = 1`). All 120
  reorderings of a five-step schedule give one value.
- **Tetration buys nothing.** From `μ₀ = 1/3`, gears `2, 4, 16, 65536, 2^65536` reach just under
  `10^19735` digits in five steps, for the same 131,118 multiplications that `p = 2` spends
  over 65,559 steps.
- **Where the map already lives.** Brahmagupta's *bhāvanā* (628 CE), Newton at `p = 2`, Halley
  at `p = 3`, Chebyshev on the hyperbola, de Moivre on the circle, the AGM and Landen in the
  nome chart.

## Citing

Cite the paper, not the repository:

```bibtex
@article{bhardwaj2026vimana,
  author = {Ram Bhardwaj},
  title  = {Bhāradvāja Vimāna: The Axis of Convergence},
  year   = {2026},
  doi    = {10.5281/zenodo.22027073},
  url    = {https://doi.org/10.5281/zenodo.22027073},
}
```

`CITATION.cff` carries the same metadata for GitHub's "Cite this repository" button.

## Licence

The code is MIT, in `LICENSE`. This repository holds code only. The paper is not distributed
here. It is licensed separately under CC BY-NC-ND 4.0, matching the Zenodo record.
