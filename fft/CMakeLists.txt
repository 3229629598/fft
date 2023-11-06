cmake_minimum_required (VERSION 3.5)

project (fft)

include_directories (include)

add_executable (fft src/fft.cpp)

# add_library (fft_shared SHARED src/fft.cpp)

# set_target_properties (fft_shared PROPERTIES OUTPUT_NAME "fft")

# set (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
