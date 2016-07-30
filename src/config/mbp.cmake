set(CMAKE_CXX_FLAGS "-DEIGEN_USE_MKL_ALL -L/opt/intel/composer_xe_2015.0.077/mkl/lib/ -mkl=sequential -I/opt/intel/composer_xe_2015.0.077/mkl/include -I/opt/local/include/eigen3/ ${CMAKE_CXX_FLAGS} -wd597 -wd1224")

#set(CMAKE_CXX_FLAGS "-I/opt/local/include/eigen3/ ${CMAKE_CXX_FLAGS}")
