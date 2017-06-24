#!/bin/bash
#This file was created on 02/02/2017


folders=`ls -d */ | egrep -v 'alanineDipeptideIsomerization|ammoniaBorane|dielsAlder|methanolFormaldehydeHydTransfer' `
for folder in $folders;do 
echo $folder
	echo '#!/bin/bash' > $folder/molpro/essm/guest.slurm.in
	echo '#SBATCH -p guest --job-name=${TEST_NAME}'>>$folder/molpro/essm/guest.slurm.in
	echo '#SBATCH --output=${TEST_NAME}.o'>>$folder/molpro/essm/guest.slurm.in
	echo '#SBATCH --error=${TEST_NAME}.e'>>$folder/molpro/essm/guest.slurm.in
	echo '#SBATCH --array=1'>>$folder/molpro/essm/guest.slurm.in
	echo '#SBATCH --nodes=1'>>$folder/molpro/essm/guest.slurm.in
	echo '#SBATCH --ntasks=1'>>$folder/molpro/essm/guest.slurm.in
	echo '#SBATCH --time=14:00:00'>>$folder/molpro/essm/guest.slurm.in
	echo ''>>$folder/molpro/essm/guest.slurm.in
	echo 'if [ -f /etc/bashrc ]; then'>>$folder/molpro/essm/guest.slurm.in
	echo '        . /etc/bashrc'>>$folder/molpro/essm/guest.slurm.in
	echo 'fi'>>$folder/molpro/essm/guest.slurm.in
	echo 'if [ -z ${SLURM_ARRAY_TASK_ID} ]; then'>>$folder/molpro/essm/guest.slurm.in
	echo '	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID"'>>$folder/molpro/essm/guest.slurm.in
	echo 'else' >> $folder/molpro/essm/guest.slurm.in
	echo '	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_JOB_ID"'>>$folder/molpro/essm/guest.slurm.in
	echo 'fi'>>$folder/molpro/essm/guest.slurm.in
	echo '' >> $folder/molpro/essm/guest.slurm.in
	echo '#run job'>>$folder/molpro/essm/guest.slurm.in
	echo '${PROJECT_BINARY_DIR}/packages/GSM/src/gsm.molpro.exe 1 1 > output'>>$folder/molpro/essm/guest.slurm.in
	echo 'exit'>>$folder/molpro/essm/guest.slurm.in
done


folders=`ls -d */ | egrep -v '2h-azir|benzene|CHD|ethylidene|nbor|ncpi|ny14|ny16|2plus2|pbuchi'`
for folder in $folders;do 
echo $folder
	echo '#!/bin/bash' > $folder/mopac/de-gsm/guest.slurm.in
	echo '#SBATCH -p guest --job-name=${TEST_NAME}'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '#SBATCH --output=${TEST_NAME}.o'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '#SBATCH --error=${TEST_NAME}.e'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '#SBATCH --array=1'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '#SBATCH --nodes=1'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '#SBATCH --ntasks=1'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '#SBATCH --time=14:00:00'>>$folder/mopac/de-gsm/guest.slurm.in
	echo ''>>$folder/mopac/de-gsm/guest.slurm.in
	echo 'if [ -f /etc/bashrc ]; then'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '        . /etc/bashrc'>>$folder/mopac/de-gsm/guest.slurm.in
	echo 'fi'>>$folder/mopac/de-gsm/guest.slurm.in
	echo 'if [ -z ${SLURM_ARRAY_TASK_ID} ]; then'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID"'>>$folder/mopac/de-gsm/guest.slurm.in
	echo 'else' >> $folder/mopac/de-gsm/guest.slurm.in
	echo '	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_JOB_ID"'>>$folder/mopac/de-gsm/guest.slurm.in
	echo 'fi'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '' >> $folder/mopac/de-gsm/guest.slurm.in
	echo '#run job'>>$folder/mopac/de-gsm/guest.slurm.in
	echo '${PROJECT_BINARY_DIR}/packages/GSM/src/gsm.mopac.exe 1 1 > output'>>$folder/mopac/de-gsm/guest.slurm.in
	echo 'exit'>>$folder/mopac/de-gsm/guest.slurm.in
done


folders=`ls -d */ | egrep -v '2h-azir|benzene|CHD|ethylidene|nbor|ncpi|ny14|ny16|2plus2|pbuchi'`
for folder in $folders;do 
echo $folder
	echo '#!/bin/bash' > $folder/qchem/de-gsm/guest.slurm.in
	echo '#SBATCH -p guest --job-name=${TEST_NAME}'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '#SBATCH --output=${TEST_NAME}.o'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '#SBATCH --error=${TEST_NAME}.e'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '#SBATCH --array=1'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '#SBATCH --nodes=1'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '#SBATCH --ntasks=1'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '#SBATCH --time=14:00:00'>>$folder/qchem/de-gsm/guest.slurm.in
	echo ''>>$folder/qchem/de-gsm/guest.slurm.in
	echo 'if [ -f /etc/bashrc ]; then'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '        . /etc/bashrc'>>$folder/qchem/de-gsm/guest.slurm.in
	echo 'fi'>>$folder/qchem/de-gsm/guest.slurm.in
	echo 'if [ -z ${SLURM_ARRAY_TASK_ID} ]; then'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID"'>>$folder/qchem/de-gsm/guest.slurm.in
	echo 'else' >> $folder/qchem/de-gsm/guest.slurm.in
	echo '	export OUTPUT_ID="$SLURM_JOB_NAME/$SLURM_JOB_ID"'>>$folder/qchem/de-gsm/guest.slurm.in
	echo 'fi'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '' >> $folder/qchem/de-gsm/guest.slurm.in
	echo '#run job'>>$folder/qchem/de-gsm/guest.slurm.in
	echo '${PROJECT_BINARY_DIR}/packages/GSM/src/gsm.qchem.exe 1 1 > output'>>$folder/qchem/de-gsm/guest.slurm.in
	echo 'exit'>>$folder/qchem/de-gsm/guest.slurm.in
done


