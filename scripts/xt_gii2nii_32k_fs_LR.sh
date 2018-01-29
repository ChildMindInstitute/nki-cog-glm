#!/usr/bin/env bash

export PATH=/home2/txu/lfcd/workbench/v1.2.3/bin_linux64:$PATH

fgii=$1
fnii=$2

wb_command -metric-convert -to-nifti ${fgii} ${fnii}
