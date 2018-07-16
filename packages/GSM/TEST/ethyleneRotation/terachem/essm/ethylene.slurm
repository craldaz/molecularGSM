#!/bin/bash
#SBATCH -p guest --job-name=ethylene
#SBATCH --output=ethylene.o
#SBATCH --error=ethylene.e
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
/home/caldaz/my_programs/molecularGSM/build_tc/packages/GSM/src/gsm.mopac.exe 1 1 > output
exit
