cmake_minimum_required(VERSION 2.8)
project(cmaes)

add_definitions(
  #-g
  -O3
  -pedantic
  -Wall
  -Wshadow
  -Wconversion
  -DNDEBUG
)

include_directories(cma-es)

# doxygen, source:
# http://permalink.gmane.org/gmane.comp.programming.tools.cmake.user/4930
FIND_PACKAGE(Doxygen)
if(DOXYGEN_FOUND)
set(DOXYGEN_INPUT doxygen.cfg)
add_custom_command(
  OUTPUT doxygen
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT} > /dev/null
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS ${DOXYGEN_INPUT}
)
add_custom_target(apidoc ALL DEPENDS doxygen)
add_custom_target(apidoc_forced
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
endif(DOXYGEN_FOUND)

add_subdirectory(cma-es)

add_executable(evo1 example1.cpp)
target_link_libraries(evo1 cmaes)
add_executable(evo2 example2.cpp)
target_link_libraries(evo2 cmaes)
