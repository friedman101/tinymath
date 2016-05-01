cmake_minimum_required(VERSION 2.8.9)
project(tinymath)
 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

file(GLOB SOURCES "src/*.c")

add_library(tinymath SHARED ${SOURCES})
target_link_libraries(tinymath m gsl gslcblas)
SET(CMAKE_C_FLAGS "-std=gnu99 -Wall -pedantic -Wextra -Werror -O3")

