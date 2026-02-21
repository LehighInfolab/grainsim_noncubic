#!/bin/bash
#
# SLURM job array for ramp-up nucleation testing.
# Each array index runs one seed with a corresponding transition count.
#

#SBATCH --job-name=grainsim_array
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=0-1                # 1 seed: 0-0, 2 seeds → 0-1, 3 seeds -> 0-2...
#SBATCH --output=/share/ceph/hawk/nhi_122121/suh222/logs/rampup/grainsim_%A_%a.out
#SBATCH --error=/share/ceph/hawk/nhi_122121/suh222/logs/rampup/grainsim_%A_%a.err
#SBATCH -p hawkcpu

# --- list of seed filenames ---
# use # to comment out seeds
SEEDS=(
  "0Seed256-256-256_T1_0000_0000000.ph"
  #"0Seed256-256-512_T1_0000_0000000.ph"
  #"0Seed256-256-1024_T1_0000_0000000.ph"
  #"0Seed256-256-2048_T1_0000_0000000.ph"
  # "0Seed256-256-4096_T1_0000_0000000.ph"
  # "0Seed256-256-8192_T1_0000_0000000.ph"
  # "0Seed512-512-512_T1_0000_0000000.ph"
  # "0Seed512-512-1024_T1_0000_0000000.ph"
  # "0Seed512-512-2048_T1_0000_0000000.ph"
)

# Transition count for each seeds
# use # to comment out transition count
TRANSITION_COUNTS=(
  200
  #200
  #42
  # 83
  #167
  #333
  # 666
  # 1332
  # 333
  #666
  #1332
)

# Seed path and output path

SEED_DIR="/share/ceph/hawk/nhi_122121/suh222/data/mc_seeds"
OUT_ROOT="/share/ceph/hawk/nhi_122121/suh222/results"

# Select seed and transition count based on SLURM array index
SEED_NAME=${SEEDS[$SLURM_ARRAY_TASK_ID]}
TRANSITION_COUNT=${TRANSITION_COUNTS[$SLURM_ARRAY_TASK_ID]}

# one output directory per seed, e.g.
INITIAL_STATE_FILE="${SEED_DIR}/${SEED_NAME}"
OUTPUT_FOLDER="${OUT_ROOT}/${SEED_NAME%.ph}/"

mkdir -p "${OUTPUT_FOLDER}"


# Logging simulation metadata
echo "Job ID:        ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Seed file:     ${INITIAL_STATE_FILE}"
echo "Transition count:  $TRANSITION_COUNT"
echo "Output folder: ${OUTPUT_FOLDER}"
echo "Hostname:           $(hostname)"
echo "Start time:         $(date)"

# Run simulation
srun grainsim.out \
    --config grainsim_config.txt \
    --initial "${INITIAL_STATE_FILE}" \
    --output "${OUTPUT_FOLDER}" \
    --transition-count "$TRANSITION_COUNT"

echo "End time:           $(date)"