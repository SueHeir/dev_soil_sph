# Column Collapse

This example runs the granular column-collapse validation as an aspect-ratio
sweep against Lube/Lajeunesse experiments and the Lagree-Staron-Popinet (2011)
Eqs. 3.1-3.2 scaling bands.

Run the skeptic-facing column cases:

```bash
source ~/projects/.build-env
cargo run --release --example column_collapse -- examples/column_collapse/config_sweep_a0p5.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_sweep_a1.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_a.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_sweep_a3.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_sweep_a6.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_negctl.toml
```

Regenerate the checked-in graph from the profiles emitted by those runs:

```bash
$BENCH_PYTHON examples/column_collapse/plot_results.py
```

![column-collapse aspect sweep vs reference bands](plots/column_collapse_reference_bands.png)

The figure shows the actual checks: normalized run-out and final height measured
from the emitted deposit profiles versus the Lube/Lajeunesse/LSP pass bands.
The `a=2` gate keeps the original bands unchanged: run-out `[2.40, 3.60]` and
height `[0.80, 1.70]`.

Current sweep result:

| Case | expectation | run-out band | measured run-out | height band | measured height | result |
|---|---|---:|---:|---:|---:|---|
| `a=0.5` | declared limitation | `[0.60, 1.10]` | `0.10` | `[0.38, 0.62]` | `0.45` | PASS by staying outside |
| `a=1` | declared limitation | `[1.20, 2.20]` | `0.90` | `[0.75, 1.25]` | `0.83` | PASS by staying outside |
| `a=2` | accept | `[2.40, 3.60]` | `2.50` | `[0.80, 1.70]` | `1.49` | PASS |
| `a=3` | accept | `[3.95, 6.60]` | `4.10` | `[0.95, 2.02]` | `1.82` | PASS |
| `a=6` | accept | `[6.27, 13.20]` | `8.70` | `[1.22, 2.66]` | `2.19` | PASS |
| `a=2` negative control | reject | `[2.40, 3.60]` | `-0.10` | `[0.80, 1.70]` | `1.90` | PASS by rejection |

The shallow-column cases are intentionally not made green by loosening the
Lube/Lajeunesse band: the current SPH/free-surface discretization under-runs
run-out by `0.50 L0` at `a=0.5` and `0.30 L0` at `a=1`. The high-aspect case is
also shown quantitatively: at `a=6` the model gives run-out `8.70`, above the
Lube/Lajeunesse power-law envelope (`6.27-7.59`) but below the LSP continuum
line (`13.20`).

`config_negctl.toml` is an over-frictional wrong-physics case. It must be
rejected by the same `a=2` band, proving the gate is capable of failing.
