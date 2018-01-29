#!/usr/bin/env bash
Hemi=$1
prefix=$2
outdir=$3

if [[ ${Hemi} == "L" ]]; then
    hemisphere=lh
elif [[ ${Hemi} == "R" ]]; then
    hemisphere=rh
else
    echo >&2 "Must enter Hemi as L or R"; 
    exit 1;
fi

export PATH=/home2/txu/lfcd/workbench/v1.2.3/bin_linux64:$PATH

templateDir=/home/txu/lfcd/hcpPipelines/xt/templates/32k_fs_LR
groupDir=/projects/txu/NKI_lifespan/group

cwd=$( pwd )
workdir=/projects/txu/NKI_lifespan/group/boundary/GLM/gradient_age_age2_gm_cov/residual_cognitive/
cd ${outdir}

pwd
mkdir clusters

p_uncorr=0.05
area=344

wb_command -metric-find-clusters ${templateDir}/Conte69.${Hemi}.midthickness.32k_fs_LR.surf.gii ${prefix}.${Hemi}.32k_fs_LR.func.gii ${p_uncorr} ${area} clusters/${Hemi}_${prefix}_cluster.func.gii -less-than -roi ${groupDir}/masks/${Hemi}.brain.NKI323.wb.32k_fs_LR.shape.gii

wb_command -metric-math 'pmap*(mask>0)' clusters/${Hemi}_${prefix}_cluster_marked.func.gii -var pmap ${prefix}.${Hemi}.32k_fs_LR.func.gii -var mask clusters/${Hemi}_${prefix}_cluster.func.gii
wb_command -metric-convert -to-nifti clusters/${Hemi}_${prefix}_cluster.func.gii clusters/${hemisphere}.${prefix}_cluster.nii.gz
wb_command -metric-convert -to-nifti clusters/${Hemi}_${prefix}_cluster_marked.func.gii clusters/${hemisphere}.${prefix}.nii.gz


#for invar in age age2 meanFD sex gm; do
#  for Hemi in L R; do
#    echo ${prefix}
#    wb_command -metric-math 'pmap*(mask>0)' clusters/${Hemi}_${invar}_t_value_cluster_marked.func.gii -var pmap ${invar}_t_value.${Hemi}.32k_fs_LR.func.gii -var mask clusters/${Hemi}_${invar}_log_p_cluster.func.gii 
#  done
#done

