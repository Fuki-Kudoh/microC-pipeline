# Configuration directory

This directory contains example configuration files for the lightweight config-driven single-sample Micro-C runner.

- `example.single-sample.yaml` documents the v0.4.0 single-sample config format.

Validate a config before launching a full run:

```bash
bin/wulf3c validate-config --config config/example.single-sample.yaml
```

Print the planned command sequence without running external tools:

```bash
bin/wulf3c run --config config/example.single-sample.yaml --dry-run
```

The default v0.4.0 workflow is Micro-C-oriented and does not require restriction enzyme information.
