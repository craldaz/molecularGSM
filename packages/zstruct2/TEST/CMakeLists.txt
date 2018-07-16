# TEST_DIR is passed as an argument to GFSTRINGQ_TEST function
# TEST Dir is where the input files are kept
FUNCTION(ZStruct_TEST TEST_DIR TEST_NAME)

 		#copy input files to the testing folder
		file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})

ENDFUNCTION()

#call the function
ZStruct_TEST(ammoniaborane ammoniaborane)

#if(GSM_ENABLE_QCHEM)
#    MESSAGE(STATUS "QCHEM is enabled.")
#elseif(GSM_ENABLE_QCHEM_SF)
#    MESSAGE(STATUS "QCHEM_SF is enabled.")
#elseif(GSM_ENABLE_GAUSSIAN)
#    MESSAGE(STATUS "GAUSSIAN is enabled.")
#elseif(GSM_ENABLE_ORCA)
#    MESSAGE(STATUS "ORCA is enabled.")
#elseif(GSM_ENABLE_MOLPRO)
#    MESSAGE(STATUS "MOLPRO is enabled.")
#elseif (GSM_ENABLE_TC)
#    MESSAGE(STATUS "TC is enabled.")
#else()
#    MESSAGE(STATUS "MOPAC is enabled.")
#endif()
