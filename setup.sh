#!/bin/bash

set -euo pipefail

# Create folers 
mkdir -p /input
mkdir -p /output
mkdir -p /slurm-log

# Activate Conda


# Create environment for YAML file
conda env create -f environment.yml

# Activate PhyLoad environment
conda activate phyload


