cmake_minimum_required (VERSION 3.5)

project (fft)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

include_directories (include)

add_executable (fft src/fft.cpp)
