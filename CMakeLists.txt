cmake_minimum_required(VERSION 3.0.2)
project(ACSE-4-project-3)
include(CTest)

set(CMAKE_CXX_STANDARD 17)

#find_package(MPI REQUIRED)

#include_directories(${MPI_INCLUDE_PATH})

if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Compling with debug flags")
else()
    message("Release mode => compling with assertions off and -O3 to improve performance")
    add_compile_options ( -DNDEBUG -O3 )
endif()

#message("Adding -fPIC compile option")
#add_compile_options ( -fPIC )

include_directories(includes)
add_subdirectory(src)
add_subdirectory(tests)

add_executable(SPH_main main.cpp)
target_link_libraries(SPH_main SPH file_writer)

#target_link_libraries(life_mpi life comms ${MPI_LIBRARIES})

# if(MPI_COMPILE_FLAGS)
#   set_target_properties(life_mpi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
# endif()

# if(MPI_LINK_FLAGS)
#   set_target_properties(life_mpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
# else()
#   message(FATAL_ERROR "MPI link flags not set by find_package(MPI)")
# endif()


