target_sources(scalloplib PRIVATE trollop_classes.cc trollop.cc)
target_sources(scalloplib PUBLIC FILE_SET scalloplib_headers TYPE HEADERS FILES trollop.h)
target_include_directories(scalloplib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
