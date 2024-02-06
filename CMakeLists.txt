cmake_minimum_required(VERSION 3.23)
set(CMAKE_CXX_STANDARD 11)
project(scallop)

find_path(GMP_INCLUDE gmp.h REQUIRED)
find_library(GMP_LIB NAMES libgmp.a gmp REQUIRED)
message(STATUS "Using GMP library ${GMP_LIB}")
find_path(GLPK_INCLUDE glpk.h REQUIRED)
find_library(GLPK_LIB NAMES libglpk.a glpk REQUIRED)
message(STATUS "Using GLPK library ${GLPK_LIB}")

add_library(scalloplib lp.cc rational.cc word.cc)
target_sources(scalloplib PUBLIC FILE_SET scalloplib_headers TYPE HEADERS FILES lp.h rational.h word.h)
install(TARGETS scalloplib ARCHIVE FILE_SET scalloplib_headers)

target_include_directories(scalloplib PUBLIC ${GMP_INCLUDE})
target_link_libraries(scalloplib PUBLIC ${GMP_LIB})
target_include_directories(scalloplib PRIVATE ${GLPK_INCLUDE})
target_link_libraries(scalloplib PRIVATE ${GLPK_LIB})

if (USE_GUROBI)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
    find_package(GUROBI REQUIRED)
    target_compile_definitions(scalloplib PRIVATE GUROBI_INSTALLED)
    target_include_directories(scalloplib PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(scalloplib ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
endif()

set(subdirs exlp-package gallop hallop scabble scylla trollop)
foreach (subdir ${subdirs})
    add_subdirectory(${subdir})
endforeach()

add_executable(scallop scallop.cc)
target_link_libraries(scallop PUBLIC scalloplib)
install(TARGETS scallop RUNTIME)

