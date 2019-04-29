cmake_minimum_required(VERSION 3.9)
set(CMAKE_CXX_STANDARD 17)

project(RLANG_NATIVE)

set(ANACONDA_PATH /home/archbox/opt/anaconda3/)

#----------------------------------------------

add_library(rcxx SHARED rcxx.cpp)

set_target_properties(rcxx PROPERTIES
    PREFIX ""
    )

target_include_directories(rcxx PRIVATE
  ${ANACONDA_PATH}/lib/R/include/
  ${ANACONDA_PATH}/lib/R/library/Rcpp/include
)


add_custom_command(TARGET rcxx POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  $<TARGET_FILE:rcxx>                 # Source
  ${CMAKE_CURRENT_LIST_DIR}/          # Destination, file or directory
  # OR: ${CMAKE_CURRENT_SOURCE_DIR}/  # Directory destination.
)
