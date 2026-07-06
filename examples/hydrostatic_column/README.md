# `hydrostatic_column` — hydrostatic pressure and tensile-stability check

A granular-SPH column settles under gravity onto frozen boundary layers. At rest,
the pressure should rise linearly with depth, with `dp/dz = -rho*g`, and the
column should stay compressive and smooth rather than developing the tensile
instability described by Bui et al. (2008).

## Reference

- Bui, Fukagawa, Sako & Ohno (2008), *IJNAMG* **32**:1537,
  DOI [10.1002/nag.688](https://doi.org/10.1002/nag.688), for the tensile
  instability failure mode: spurious tension and particle clustering.
- Hydrostatic equilibrium gives the pressure-gradient reference
  `dp/dz = -rho*g` with `rho = 1500 kg/m^3` and `g = 9.81 m/s^2`.

## Run and plot

```sh
source ~/projects/.build-env
$BENCH_PYTHON examples/hydrostatic_column/sweep.py
```

The script runs the example, parses its own PASS output, and regenerates the
checked-in figure below.

![hydrostatic column validation](plots/hydrostatic_column_validation.svg)

Pressure slab means are compared with the hydrostatic reference, the accepted
`0.7-1.3 x -rho*g` gradient band, and the regression floor; the same run also
shows the tensile-instability pressure and density-spread gates passing.

Latest measured result: `dp/dz / (-rho*g) = 0.8208`, `p_min = 8.21 Pa`
against the `-3.31 Pa` limit, density spread `0.016%` against the `2%` limit.
Exit 0 = PASS, nonzero = FAIL.
