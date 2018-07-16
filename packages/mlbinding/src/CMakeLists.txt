
# source files
file(GLOB MLB_SRC
    "*.h"
    "*.cpp"
)

# set up libraries
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MKL_LIBRARIES})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# set target name
TRIBITS_ADD_EXECUTABLE(${TARGET_NAME3} NOEXEPREFIX SOURCES ${MLB_SRC} IMPORTEDLIBS ${LIBS})
