set(my_sources
    EXLPvector.c
    eta_file.c
    hash.c
    lpstruct.c
    lu.c
    main.c
    matrix.c
    mps.c
    mylib.c
    preprocess.c
    solve_ip.c
    solve_lp.c
)
target_compile_definitions(scalloplib PRIVATE HAVE_CONFIG_H NO_GMP_HASH)
target_sources(scalloplib PRIVATE ${my_sources})
target_include_directories(scalloplib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
