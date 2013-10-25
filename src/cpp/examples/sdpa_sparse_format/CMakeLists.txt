project(sdpa_sparse_format)
add_executable(sdpa_sparse_format "driver.cpp")
target_link_libraries(sdpa_sparse_format
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        driver.cpp
        sdpa_sparse_format.json
        sdpa_sparse_format_phase1.json
        lp.dat-s
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
    install(TARGETS sdpa_sparse_format
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
endif()

# Run some unit tests
add_sdpa_test("unit/*.dat-s" "unit/*.json" sdpa_sparse_format)
