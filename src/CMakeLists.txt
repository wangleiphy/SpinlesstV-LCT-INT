cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(spinlessctint NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_ROOT} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

if (USE_MACHINE)
  if (EXISTS ${USE_MACHINE})
    # given full-path case
    set(MACHINE_CONF ${USE_MACHINE})
    message("-- Loading machine configuration: ${MACHINE_CONF}")
    include(${MACHINE_CONF})
  else (EXISTS ${USE_MACHINE})
    # given machine-name case
    set(MACHINE_CONF "${CMAKE_CURRENT_SOURCE_DIR}/config/${USE_MACHINE}.cmake")
    if (EXISTS ${MACHINE_CONF})
      message("-- Loading machine configuration: ${MACHINE_CONF}")
      include(${MACHINE_CONF})
    else (EXISTS ${MACHINE_CONF})
      message(FATAL_ERROR "Can't find your machine-configuration file")
    endif (EXISTS ${MACHINE_CONF})
  endif (EXISTS ${USE_MACHINE})
else(USE_MACHINE)
    message(FATAL_ERROR "where am I?")
endif(USE_MACHINE)

# enable C++ compilers
enable_language(CXX)

set (SOURCE
     fastupdate.cpp
     solver.cpp 
     io.cpp 
     interaction_expansion.cpp 
     model.cpp
     observables.cpp 
     measure.cpp 
     )

add_executable(main ${SOURCE} main.cpp)
target_link_libraries(main ${ALPS_LIBRARIES})

install(TARGETS main DESTINATION bin)
