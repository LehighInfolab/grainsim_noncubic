# grainsim
Monte Carlo grain coarsening simulation with support for transformed grain boundaries.

## Build
### Linux
Compile using the provided Makefile:
```bash
make
```
This produces the executable:
```bash
grainsim.out
```

### Windows
Run:
```bash
compile_windows.bat
```
## Usage
By default, the simulation reads parameters from a configuration file.

Basic usage with optional parameters:
```bash
./grainsim.out [--config <file>] [--initial <seed.ph>] [--output <folder>] [--transition-count <int>]
```
If ```--config ``` is not provided, the program uses the default configuration file (if defined in the source).

## Configuration
Simulation parameters are defined in a configuration file.

### Custom Z-Plane Propagation
Support for custom z-plane propagation can be enabled by specifying:
```ini
Z_PROP_PLANE = 1024
```
Add this line to your configuration file to control propagation behavior.

### Ramp-up nucleation rate (Discrete)
If enabled, the transition count starts at 0, increases by a fixed increment every TRANSITION_INTERVAL timesteps, capped at a maximum limit.

```ini
# Enable/disable ramp-up
RAMP_UP_ENABLED = true

# How often to transition boundaries (in timesteps / MCS)
TRANSITION_INTERVAL = 200

# How much to increase the transition count each interval
TRANSITION_TIME_COUNT_INCREASE = 5

# Maximum transition count per interval during ramp-up
HIGH_NUCLEATION_RATE_LIMIT = 200
```

## Example: SLURM Execution
Example SLURM command with all overrides:
```bash
srun ./grainsim.out \
    --config grainsim_config.txt \
    --initial "${INITIAL_STATE_FILE}" \
    --output "${OUTPUT_FOLDER}" \
    --transition-count "$TRANSITION_COUNT"

```
## Changelog
### 2026-02-21
### Added
- Configurable discrete ramp-up nucleation rate.
- New configuration parameters:
  - RAMP_UP_ENABLED
  - TRANSITION_TIME_COUNT_INCREASE
  - HIGH_NUCLEATION_RATE_LIMIT
- Example SLURM array script for ramp-up testing.

### 2026-02-11
- Added support for custom z-plane propagation via Z_PROP_PLANE in the configuration file.