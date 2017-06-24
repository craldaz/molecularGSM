# TEST_DIR is passed as an argument to GFSTRINGQ_TEST function
# TEST Dir is where the input files are kept
FUNCTION(GFSTRINGQ_TEST TEST_DIR TEST_NAME)

 		#copy input files to the testing folder
		file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/ DESTINATION
${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})

    #configure the PBS script to know where the executable is.
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/athena.pbs.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/${TEST_NAME}.pbs")
    #configure the SLURM script to know where the executable is.
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/guest.slurm.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/${TEST_NAME}.slurm")


		#create a bash script to run all test cases simultaneusly
    execute_process(COMMAND bash -c "echo 'cd ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/' >>${PROJECT_BINARY_DIR}/gsm.sh")
    #execute_process(COMMAND bash -c "echo 'qsub ${TEST_NAME}.pbs'>>${PROJECT_BINARY_DIR}/gsm.sh")
    #execute_process(COMMAND bash -c "echo 'cd ${PROJECT_BINARY_DIR}'>>${PROJECT_BINARY_DIR}/gsm.sh")
    execute_process(COMMAND bash -c "echo 'sbatch ${TEST_NAME}.slurm'>>${PROJECT_BINARY_DIR}/gsm.sh")
    execute_process(COMMAND bash -c "echo 'cd ${PROJECT_BINARY_DIR}'>>${PROJECT_BINARY_DIR}/gsm.sh")


	#configure the inpfileq. CAN ADD EXTRA ARGUMENTS TO GFSTRINGQ_TEST FUNCTION
	#AND PLACE THEM IN INPFILEQ TO BE CONFIGURED.
	CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/inpfileq.in" "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/inpfileq")	

	#CTEST tests works but better to run bash script
	TRIBITS_ADD_ADVANCED_TEST( 
    ${TEST_NAME}
    TEST_0
        MESSAGE "-- Working directory is
        ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}"
        # the execution command
        EXEC gsm.${CALCULATOR} NOEXEPREFIX 
				DIRECTORY ${PROJECT_BINARY_DIR}/${PROJECT_NAME}
        # the args for gstringq, 1 is the file name and 4 is the number of processors
        ARGS 1 4
        # where it executes
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}
        # don't clean it because we copied things there that it needs
        SKIP_CLEAN_WORKING_DIRECTORY
        OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/scratch/test_output
        NO_ECHO_OUTPUT
        PASS_ANY
 				COMM serial
        TIMEOUT 100000

    TEST_1
        MESSAGE "-- Comparing output in ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}"
        # the execution command
				CMND python
        # the args for CMND
        ARGS compare.py
        # where it executes
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}
        SKIP_CLEAN_WORKING_DIRECTORY
        # where the output should go
        OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}/scratch/compare_output
        #NO_ECHO_OUTPUT
        # don't clean it because we copied things there that it needs
        SKIP_CLEAN_WORKING_DIRECTORY
        FAIL_FAST
        NUM_TOTAL_CORES_USED 1
        CATEGORIES BASIC
        COMM serial
	)
ENDFUNCTION()


#configure the top level PBS
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/gsm.sh.in"
"${CMAKE_BINARY_DIR}/gsm.sh")
execute_process(COMMAND bash -c "chmod +x ${CMAKE_BINARY_DIR}/gsm.sh")


#call the function
if(GSM_ENABLE_QCHEM)
    MESSAGE(STATUS "QCHEM is enabled.")
		#don't run because it crashes
    #GFSTRINGQ_TEST(alanineDipeptideIsomerization/qchem/de-gsm alanine)
    GFSTRINGQ_TEST(ammoniaBorane/qchem/de-gsm ammonia)
    GFSTRINGQ_TEST(dielsAlder/qchem/de-gsm diels)
		##don't run ethylene because it crashes
    GFSTRINGQ_TEST(ethyleneRotation/qchem/de-gsm ethylene)
    GFSTRINGQ_TEST(methanolFormaldehydeHydTransfer/qchem/de-gsm methanol)

elseif(GSM_ENABLE_GAUSSIAN)
    MESSAGE(STATUS "GAUSSIAN is enabled.")
    #GFSTRINGQ_TEST(alanineDipeptideIsomerization/gaussian/de-gsm alanine)
    #GFSTRINGQ_TEST(ammoniaBorane/gaussian/de-gsm ammonia)
    GFSTRINGQ_TEST(dielsAlder/gaussian/de-gsm diels)
    #GFSTRINGQ_TEST(ethyleneRotation/gaussian/de-gsm ethylene)
    #GFSTRINGQ_TEST(methanolFormaldehydeHydTransfer/gaussian/de-gsm methanol)

elseif(GSM_ENABLE_ORCA)
    MESSAGE(STATUS "ORCA is enabled.")
		#
		# Don't run alanine --causes GSM to stall
    #GFSTRINGQ_TEST(alanineDipeptideIsomerization/orca/de-gsm alanine)
   # GFSTRINGQ_TEST(ammoniaBorane/orca/de-gsm ammonia)
    GFSTRINGQ_TEST(dielsAlder/orca/de-gsm diels)
		#
		#Don't run ethene --causes GSM to stall
    #GFSTRINGQ_TEST(ethyleneRotation/orca/de-gsm ethylene)
    #GFSTRINGQ_TEST(methanolFormaldehydeHydTransfer/orca/de-gsm methanol)
elseif(GSM_ENABLE_MOLPRO)
    MESSAGE(STATUS "MOLPRO is enabled.")

		# DE-GSM
    #GFSTRINGQ_TEST(ethyleneRotation/molpro/de-gsm/ ethene_de-gsm)
    #GFSTRINGQ_TEST(ethylidene/molpro/de-gsm/ ethylidene_de-gsm)
    #GFSTRINGQ_TEST(benzene/molpro/de-gsm/ benzene_de-gsm 15 0.8)
    #GFSTRINGQ_TEST(CHD/molpro/de-gsm/ CHD_de-gsm 15 0.8)
    #GFSTRINGQ_TEST(nbor/molpro/de-gsm/ nbor_de-gsm 15 0.8)
    #GFSTRINGQ_TEST(ncpi/molpro/de-gsm/ ncpi_de-gsm 15 0.8)
    #GFSTRINGQ_TEST(ny14/molpro/de-gsm/ ny14_de-gsm 15 0.8)
    #GFSTRINGQ_TEST(ny16/molpro/de-gsm/ ny10_de-gsm 15 0.8)
    #GFSTRINGQ_TEST(2h-azir/molpro/de-gsm/ 2h-azir_de-gsm 15 0.8)

		#Conical intersection optimizer
    #GFSTRINGQ_TEST(ethyleneRotation/molpro/ciopt_only/ cethene)
    #GFSTRINGQ_TEST(2h-azir/molpro/ciopt_only/ 2h-azir)
    #GFSTRINGQ_TEST(benzene/molpro/ciopt_only/ benzene)
    #GFSTRINGQ_TEST(CHD/molpro/ciopt_only/ CHD)
    #GFSTRINGQ_TEST(ethylidene/molpro/ciopt_only/ ethylidene)
    #GFSTRINGQ_TEST(nbor/molpro/ciopt_only/ nbor)
    #GFSTRINGQ_TEST(ncpi/molpro/ciopt_only/ ncpi)
    #GFSTRINGQ_TEST(ny14/molpro/ciopt_only/ ny14)
    #GFSTRINGQ_TEST(ny16/molpro/ciopt_only/ ny16)

		#Single-Ended Excited State String Method
    #GFSTRINGQ_TEST(ethyleneRotation/molpro/essm/ ethene_essm)
    #GFSTRINGQ_TEST(ethylidene/molpro/essm/ ethylidene_essm)
    #GFSTRINGQ_TEST(benzene/molpro/essm/ benzene_essm)
    ###GFSTRINGQ_TEST(cyclobutene/molpro/essm/ cyclobutene)
    #GFSTRINGQ_TEST(CHD/molpro/essm/ CHD_essm)
    #GFSTRINGQ_TEST(nbor/molpro/essm/ nbor_essm)
    #GFSTRINGQ_TEST(ncpi/molpro/essm/ ncpi_essm)
    #GFSTRINGQ_TEST(ny14/molpro/essm/ ny14_essm)
    #GFSTRINGQ_TEST(ny16/molpro/essm/ ny16_essm)
    #GFSTRINGQ_TEST(2h-azir/molpro/essm/ 2h-azir_essm)
    #GFSTRINGQ_TEST(2plus2/molpro/essm/ 2plus2)

		#Double-Ended Excited State String Method
    #GFSTRINGQ_TEST(ethyleneRotation/molpro/map/ map_ethene)
    #GFSTRINGQ_TEST(pbuchi/molpro/map/ map_pbuchi)

else()
    MESSAGE(STATUS "MOPAC is enabled.")
		#
    GFSTRINGQ_TEST(alanineDipeptideIsomerization/mopac/de-gsm alanine)
    GFSTRINGQ_TEST(ammoniaBorane/mopac/de-gsm ammonia)
		GFSTRINGQ_TEST(dielsAlder/mopac/de-gsm diels)
    GFSTRINGQ_TEST(ethyleneRotation/mopac/de-gsm ethylene)
    GFSTRINGQ_TEST(methanolFormaldehydeHydTransfer/mopac/de-gsm methanol)
endif()