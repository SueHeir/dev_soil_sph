# Column Collapse

This example runs the granular column-collapse validation at aspect ratio
`a = H0/L0 = 2`. Configs `a`, `b`, and `c` form a resolution study against the
Lube/Lajeunesse experimental run-out band as summarized by Lagree, Staron &
Popinet (2011); `config_negctl.toml` is an over-frictional negative control that
must be rejected by the same band.

Run the validation cases:

```bash
source ~/projects/.build-env
cargo run --release --example column_collapse -- examples/column_collapse/config_a.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_b.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_c.toml
cargo run --release --example column_collapse -- examples/column_collapse/config_negctl.toml
```

Regenerate the checked-in graph from the profiles emitted by those runs:

```bash
$BENCH_PYTHON examples/column_collapse/plot_results.py
```

![column-collapse measured-vs-reference bands](plots/column_collapse_reference_bands.png)

Measured normalized run-out and deposit height from the current example output
are plotted against the Lube/Lajeunesse/LSP pass bands; configs `a`, `b`, and
`c` PASS inside the bands, while the negative control PASSes by being rejected
outside them.
