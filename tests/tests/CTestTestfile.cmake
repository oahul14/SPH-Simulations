# CMake generated Testfile for 
# Source directory: /mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests
# Build directory: /mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_file_writer "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/tests/test_file_writer")
add_test(test_SPH_forward_euler "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/tests/test_SPH" "10" "forward_euler")
add_test(test_SPH_improved_euler "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/tests/test_SPH" "10" "improved_euler")
add_test(test_SPH_AB2 "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/tests/test_SPH" "10" "AB2")
add_test(python_tests "python3" "python_tests.py")
