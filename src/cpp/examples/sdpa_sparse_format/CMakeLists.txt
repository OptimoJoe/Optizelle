project(sdpa_sparse_format)
add_executable(sdpa_sparse_format "sdpa_sparse_format.cpp")
target_link_libraries(sdpa_sparse_format
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        sdpa_sparse_format.cpp
        lp.dat-s
        lp_phase1.json
        lp_phase2.json
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
    install(TARGETS sdpa_sparse_format
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
endif()

# Run some unit tests
add_sdpa_test("*.dat-s" "*.json" sdpa_sparse_format)
