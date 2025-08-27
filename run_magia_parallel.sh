#!/bin/bash

#SBATCH --job-name=magia
#SBATCH --gres=lscratch:20
#SBATCH --time=4:00:00
#SBATCH --partition=quick
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Load necessary modules and environment
export FS_LICENSE=/data/yanw4/license.txt
module load freesurfer
module load connectome-workbench
module load matlab

# Check for exactly 7 arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: sbatch run_magia_parallel.sh <data_dir> <sbj> <pet_tr> <workfolder> <sessions> <modeling_options.mat> <specs.mat>"
    exit 1
fi

# Assign arguments to variables
DATA_DIR="$1"
SB="$2"
PET_TR="$3"
WORKFOLDER="$4"
SESSIONS="$5"
MODELING_OPTIONS_FILE="$6"
SPECS_FILE="$7"

echo "========== Starting MAGIA run =========="
echo "Subject:           $SB"
echo "Tracer:            $PET_TR"
echo "Data directory:    $DATA_DIR"
echo "Output folder:     $WORKFOLDER"
echo "Sessions:          $SESSIONS"
echo "Modeling options:  $MODELING_OPTIONS_FILE"
echo "Specs file:        $SPECS_FILE"
echo "========================================"

# Launch MATLAB job
matlab -nodisplay -nosplash -r "\
data_parent_directory = '$DATA_DIR'; \
sb = '$SB'; \
pet_tr = '$PET_TR'; \
workfolder = '$WORKFOLDER'; \
sessions = str2num('$SESSIONS'); \
modeling_options_file = '$MODELING_OPTIONS_FILE'; \
specs_file = '$SPECS_FILE'; \
run_magia_with_params(data_parent_directory, sb, pet_tr, workfolder, sessions, modeling_options_file, specs_file); \
exit;"
