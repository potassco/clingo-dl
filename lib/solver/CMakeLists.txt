# [[[header: include
set(ide_header_group "Header Files")
set(header-group-include-clingo-dl
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/config.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/graph.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/parsing.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/propagator.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/theory.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/util.hh")
source_group("${ide_header_group}\\include\\clingo-dl" FILES ${header-group-include-clingo-dl})
set(header-group-include-clingo-dl-impl
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/impl/util.hh")
source_group("${ide_header_group}\\include\\clingo-dl\\impl" FILES ${header-group-include-clingo-dl-impl})
set(header
    ${header-group-include-clingo-dl}
    ${header-group-include-clingo-dl-impl})
# ]]]
# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/graph.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/parsing.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/propagator.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/theory.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]
# [[[test: tests
set(ide_test_group "Test Files")
set(test-group-tests
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/parsing.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/propagator.cc")
source_group("${ide_test_group}\\tests" FILES ${test-group-tests})
set(test
    ${test-group-tests})
# ]]]

add_library(clingo-dl-solver STATIC)
clingodl_target_properties(TARGET clingo-dl-solver FOLDER "lib" TYPE default)
target_sources(clingo-dl-solver PRIVATE ${header} ${source})
target_include_directories(clingo-dl-solver PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(clingo-dl-solver PUBLIC clingo::clingo)
if (CLINGODL_PROFILE)
    target_compile_definitions(clingo-dl-solver PUBLIC CLINGODL_PROFILE)
    target_link_libraries(clingo-dl-solver PUBLIC Gperftools::Profiler)
endif()

if (CLINGODL_BUILD_TESTS)
    add_library(clingo-dl-solver_t STATIC)
    clingodl_target_properties(TARGET clingo-dl-solver_t FOLDER "lib" TYPE default)
    target_sources(clingo-dl-solver_t PRIVATE ${header} ${source})
    target_link_libraries(clingo-dl-solver_t PUBLIC clingo::clingo)
    target_include_directories(clingo-dl-solver_t PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
    target_compile_definitions(clingo-dl-solver_t PUBLIC CLINGODL_CROSSCHECK)

    add_executable(test_solver)
    clingodl_target_properties(TARGET test_solver FOLDER "test" SUBDIR "test")
    target_sources(test_solver PRIVATE ${test})
    target_link_libraries(test_solver PRIVATE clingo-dl-solver_t Catch2::Catch2WithMain)
    catch_discover_tests(test_solver)
endif()
