# grainsim
Monte Carlo grain coarsening simulation with support for transformed grain boundaries.

To compile on Windows, simply run "compile_windows.bat". To compile on Linux, run "make" to compile via the Makefile.

# optional config override
Simulation parameters are read from a config file, but key values may be overridden with:
```bash
grainsim [--config <file>] [--initial <seed.ph>] [--output <folder>] [--transition-count <int>]
```
example slurm config with all overrides:
```bash
srun grainsim.out \
    --config grainsim_config.txt \
    --initial "${INITIAL_STATE_FILE}" \
    --output "${OUTPUT_FOLDER}" \
    --transition-count "$TRANSITION_COUNT"
```
