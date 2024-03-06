target_sources(scalloplib PRIVATE gallop.cc graph.cc)
target_sources(scalloplib PUBLIC FILE_SET scalloplib_headers TYPE HEADERS FILES gallop.h graph.h)
