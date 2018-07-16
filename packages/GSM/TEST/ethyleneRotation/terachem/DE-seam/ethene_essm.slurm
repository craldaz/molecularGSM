#!/bin/bash
#SBATCH -p guest --job-name=ethene_essm
#SBATCH --output=ethene_essm.o
#SBATCH --error=ethene_essm.e
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=14:00:00

if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi
if [ -z  ]; then
	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID"
else
	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_JOB_ID"
fi

#run job
/export/zimmerman/craldaz/ESSM/map/trymap6/packages/GSM/src/gsm.molpro.exe 1 1 > output
exit
