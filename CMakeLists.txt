cmake_minimum_required(VERSION 3.23)
project(scallop)

find_path(GMP_INCLUDE gmp.h REQUIRED)
find_library(GMP_LIB gmp REQUIRED)
find_path(GLPK_INCLUDE glpk.h REQUIRED)
find_library(GLPK_LIB glpk REQUIRED)

add_library(scalloplib lp.cc rational.cc word.cc)
target_sources(scalloplib PUBLIC FILE_SET scalloplib_headers TYPE HEADERS FILES lp.h rational.h word.h)
install(TARGETS scalloplib ARCHIVE FILE_SET scalloplib_headers)

target_include_directories(scalloplib PUBLIC ${GMP_INCLUDE})
target_link_libraries(scalloplib PUBLIC ${GMP_LIB})
target_include_directories(scalloplib PRIVATE ${GLPK_INCLUDE})
target_link_libraries(scalloplib PRIVATE ${GLPK_LIB})

set(subdirs exlp-package gallop hallop scabble scylla trollop)
foreach (subdir ${subdirs})
    add_subdirectory(${subdir})
endforeach()

add_executable(scallop scallop.cc)
target_link_libraries(scallop PUBLIC scalloplib)
install(TARGETS scallop RUNTIME)

