add_executable(test_file_writer test_file_writer.cpp)
target_link_libraries(test_file_writer file_writer SPH)
add_test(NAME test_file_writer COMMAND test_file_writer)

add_executable(test_SPH test_SPH_2D.cpp)
target_link_libraries(test_SPH file_writer SPH)
add_test(NAME test_SPH COMMAND test_SPH)

configure_file(python_tests.py python_tests.py COPYONLY)
add_test(NAME python_tests COMMAND python3 python_tests.py)
