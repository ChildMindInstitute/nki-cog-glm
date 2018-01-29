#!/usr/bin/env bash

export PATH=/home2/txu/lfcd/workbench/v1.2.3/bin_linux64:$PATH
hcptemplateDir=/home/txu/lfcd/hcpPipelines/xt/templates/32k_fs_LR

fnii=$1
fgii=$2
Hemi=$3

wb_command -metric-convert -from-nifti ${fnii} $hcptemplateDir/Conte69.${Hemi}.midthickness.32k_fs_LR.surf.gii ${fgii}
