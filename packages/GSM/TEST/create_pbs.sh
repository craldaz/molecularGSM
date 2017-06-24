#!/bin/bash
#This file was created on 02/02/2017

echo "Making ciopt sub files"
folders=`ls -d */ | egrep -v 'alanineDipeptideIsomerization|ammoniaBorane|dielsAlder|methanolFormaldehydeHydTransfer|2plus2|pbuchi'`
for folder in $folders;do 
	echo $folder
	echo "#!/bin/bash" > $folder/molpro/ciopt_only/athena.pbs.in
	echo "#PBS -l nodes=1:ppn=1 -l walltime=18:00:00">> $folder/molpro/ciopt_only/athena.pbs.in
	echo '#PBS -N ${TEST_NAME}'>> $folder/molpro/ciopt_only/athena.pbs.in
	echo '#PBS -o ${TEST_NAME}'>> $folder/molpro/ciopt_only/athena.pbs.in
	echo '#PBS -e ${TEST_NAME}'>> $folder/molpro/ciopt_only/athena.pbs.in
	echo 'cd $PBS_O_WORKDIR'>> $folder/molpro/ciopt_only/athena.pbs.in
	echo '${PROJECT_BINARY_DIR}/packages/GSM/src/gsm.molpro.exe 1 1 > output'>> $folder/molpro/ciopt_only/athena.pbs.in
done

echo "Making map sub files"
folders=`ls -d */ | egrep -v '2h-azir|benzene|CHD|ethylidene|nbor|ncpi|ny14|ny16|alanineDipeptideIsomerization|ammoniaBorane|dielsAlder|methanolFormaldehydeHydTransfer|2plus2'`
for folder in $folders;do 
echo $folder
	echo "#!/bin/bash" > $folder/molpro/map/athena.pbs.in
	echo "#PBS -l nodes=1:ppn=1 -l walltime=18:00:00">> $folder/molpro/map/athena.pbs.in
	echo '#PBS -N ${TEST_NAME}'>> $folder/molpro/map/athena.pbs.in
	echo '#PBS -o ${TEST_NAME}'>> $folder/molpro/map/athena.pbs.in
	echo '#PBS -e ${TEST_NAME}'>> $folder/molpro/map/athena.pbs.in
	echo 'cd $PBS_O_WORKDIR'>> $folder/molpro/map/athena.pbs.in
	echo '${PROJECT_BINARY_DIR}/packages/GSM/src/gsm.molpro.exe 1 1 > output'>> $folder/molpro/map/athena.pbs.in
done

echo "Making qchem de-gsm sub files"
folders=`ls -d */ | egrep -v '2h-azir|benzene|CHD|ethylidene|nbor|ncpi|ny14|ny16|2plus2|pbuchi'`
for folder in $folders;do 
echo $folder
	echo "#!/bin/bash" > $folder/qchem/de-gsm/athena.pbs.in
	echo "#PBS -l nodes=1:ppn=1 -l walltime=18:00:00">> $folder/qchem/de-gsm/athena.pbs.in
	echo '#PBS -N ${TEST_NAME}'>> $folder/qchem/de-gsm/athena.pbs.in
	echo '#PBS -o ${TEST_NAME}'>> $folder/qchem/de-gsm/athena.pbs.in
	echo '#PBS -e ${TEST_NAME}'>> $folder/qchem/de-gsm/athena.pbs.in
	echo 'cd $PBS_O_WORKDIR'>> $folder/qchem/de-gsm/athena.pbs.in
	echo '${PROJECT_BINARY_DIR}/packages/GSM/src/gsm.qchem.exe 1 1 > output'>> $folder/qchem/de-gsm/athena.pbs.in
done

echo "Making mopac de-gsm sub files"
folders=`ls -d */ | egrep -v '2h-azir|benzene|CHD|ethylidene|nbor|ncpi|ny14|ny16|2plus2|pbuchi'`
for folder in $folders;do 
echo $folder
	echo "#!/bin/bash" > $folder/mopac/de-gsm/athena.pbs.in
	echo "#PBS -l nodes=1:ppn=1 -l walltime=18:00:00">> $folder/mopac/de-gsm/athena.pbs.in
	echo '#PBS -N ${TEST_NAME}'>> $folder/mopac/de-gsm/athena.pbs.in
	echo '#PBS -o ${TEST_NAME}'>> $folder/mopac/de-gsm/athena.pbs.in
	echo '#PBS -e ${TEST_NAME}'>> $folder/mopac/de-gsm/athena.pbs.in
	echo 'cd $PBS_O_WORKDIR'>> $folder/mopac/de-gsm/athena.pbs.in
	echo '${PROJECT_BINARY_DIR}/packages/GSM/src/gsm.mopac.exe 1 1 > output'>> $folder/mopac/de-gsm/athena.pbs.in
done


