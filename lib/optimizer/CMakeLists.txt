# [[[header: include
set(ide_header_group "Header Files")
set(header-group-include-clingo-dl
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl/optimizer.hh")
source_group("${ide_header_group}\\include\\clingo-dl" FILES ${header-group-include-clingo-dl})
set(header
    ${header-group-include-clingo-dl})
# ]]]
# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/optimizer.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]
# [[[test: tests
set(ide_test_group "Test Files")
set(test-group-tests
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/optimize.cc")
source_group("${ide_test_group}\\tests" FILES ${test-group-tests})
set(test
    ${test-group-tests})
# ]]]

add_library(clingo-dl-optimizer STATIC)
clingodl_target_properties(TARGET clingo-dl-optimizer FOLDER "lib" TYPE extra)
target_sources(clingo-dl-optimizer PRIVATE ${header} ${source})
target_link_libraries(clingo-dl-optimizer PUBLIC clingo-dl)
target_include_directories(clingo-dl-optimizer PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")

if (CLINGODL_BUILD_TESTS)
    add_library(clingo-dl-optimizer_t STATIC)
    clingodl_target_properties(TARGET clingo-dl-optimizer_t FOLDER "lib")
    target_sources(clingo-dl-optimizer_t PRIVATE ${header} ${source})
    target_include_directories(clingo-dl-optimizer_t PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
    target_link_libraries(clingo-dl-optimizer_t PUBLIC clingo-dl_t)

    add_executable(test_optimizer)
    clingodl_target_properties(TARGET test_optimizer FOLDER "test" SUBDIR "tests")
    target_sources(test_optimizer PRIVATE ${test})
    target_link_libraries(test_optimizer PRIVATE clingo-dl-optimizer_t Catch2::Catch2WithMain)
    catch_discover_tests(test_optimizer)
endif()
