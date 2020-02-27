add_library(file_writer SHARED file_writer.cpp)

add_library(SPH SHARED SPH_2D.cpp)
target_link_libraries(SPH OpenMP::OpenMP_CXX)