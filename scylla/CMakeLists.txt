target_sources(scalloplib PRIVATE scylla_classes.cc scylla.cc)
target_sources(scalloplib PUBLIC FILE_SET scalloplib_headers TYPE HEADERS FILES scylla.h)
target_include_directories(scalloplib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
