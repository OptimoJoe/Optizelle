project(linear_algebra)

include_directories(${JSONCPP_INCLUDE_DIRS})
include_directories(../utility)

add_executable(gmres_full gmres_full.cpp)
add_executable(gmres_left_preconditioner gmres_left_preconditioner.cpp)
add_executable(gmres_restart gmres_restart.cpp)
add_executable(gmres_right_preconditioner gmres_right_preconditioner.cpp)
add_executable(tcd_basic tcd_basic.cpp)
add_executable(tcd_cp tcd_cp.cpp)
add_executable(tcd_nullspace_solve tcd_nullspace_solve.cpp)
add_executable(tcd_starting_solution tcd_starting_solution.cpp)
add_executable(tcd_tr_stopping tcd_tr_stopping.cpp)
add_executable(tcd_tr_stopping_moved_center tcd_tr_stopping_moved_center.cpp)
add_executable(tminres_basic tminres_basic.cpp)
add_executable(tminres_cp tminres_cp.cpp)
add_executable(tminres_nullspace_solve tminres_nullspace_solve.cpp)
add_executable(tminres_tr_stopping tminres_tr_stopping.cpp)
add_executable(tminres_tr_stopping_moved_center
    tminres_tr_stopping_moved_center.cpp)

target_link_libraries(gmres_full
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(gmres_left_preconditioner
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(gmres_restart
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(gmres_right_preconditioner
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tcd_basic
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tcd_cp
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tcd_nullspace_solve
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tcd_starting_solution
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tcd_tr_stopping
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tcd_tr_stopping_moved_center
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tminres_basic
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tminres_cp
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tminres_nullspace_solve
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tminres_tr_stopping
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(tminres_tr_stopping_moved_center
    peopt_static ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_test(linear_algebra/gmres_full gmres_full)
add_test(linear_algebra/gmres_left_preconditioner gmres_left_preconditioner)
add_test(linear_algebra/gmres_restart gmres_restart)
add_test(linear_algebra/gmres_right_preconditioner gmres_right_preconditioner)
add_test(linear_algebra/tcd_basic tcd_basic)
add_test(linear_algebra/tcd_cp tcd_cp)
add_test(linear_algebra/tcd_nullspace_solve tcd_nullspace_solve)
add_test(linear_algebra/tcd_starting_solution tcd_starting_solution)
add_test(linear_algebra/tcd_tr_stopping tcd_tr_stopping)
add_test(linear_algebra/tcd_tr_stopping_moved_center
    tcd_tr_stopping_moved_center)
add_test(linear_algebra/tminres_basic tminres_basic)
add_test(linear_algebra/tminres_cp tminres_cp)
add_test(linear_algebra/tminres_nullspace_solve tminres_nullspace_solve)
add_test(linear_algebra/tminres_tr_stopping tminres_tr_stopping)
add_test(linear_algebra/tminres_tr_stopping_moved_center
    tminres_tr_stopping_moved_center)
