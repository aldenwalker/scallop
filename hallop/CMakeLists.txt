target_sources(scalloplib PRIVATE free_group_chain.cc hallop_lp.cc hallop.cc hyp.cc pieces.cc)
target_sources(scalloplib PUBLIC FILE_SET scalloplib_headers TYPE HEADERS FILES hallop.h)
target_include_directories(scalloplib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
