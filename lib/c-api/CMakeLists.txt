# [[[header: include
set(ide_header_group "Header Files")
set(header-group-include
    "${CMAKE_CURRENT_SOURCE_DIR}/include/clingo-dl.h")
source_group("${ide_header_group}\\include" FILES ${header-group-include})
set(header
    ${header-group-include})
# ]]]
# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/clingo-dl.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]
# [[[test: tests
set(ide_test_group "Test Files")
set(test-group-tests
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/solve.cc")
source_group("${ide_test_group}\\tests" FILES ${test-group-tests})
set(test
    ${test-group-tests})
# ]]]

if (CLINGODL_BUILD_SHARED)
    set(lib_type SHARED)
else()
    set(lib_type STATIC)
endif()

add_library(clingo-dl ${lib_type})
clingodl_target_properties(TARGET clingo-dl FOLDER "lib" TYPE default)
target_sources(clingo-dl PRIVATE ${header} ${source})
target_include_directories(clingo-dl PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(clingo-dl PUBLIC clingo-dl-solver)
if (CLINGODL_USE_VISIBILITY)
    target_compile_definitions(clingo-dl-solver PRIVATE "$<BUILD_INTERFACE:CLINGODL_BUILD_LIBRARY>")
endif()

if (CLINGODL_BUILD_TESTS)
    add_library(clingo-dl_t STATIC)
    clingodl_target_properties(TARGET clingo-dl_t FOLDER "lib" TYPE default)
    target_sources(clingo-dl_t PRIVATE ${header} ${source})
    target_link_libraries(clingo-dl_t PUBLIC clingo-dl-solver)
    target_include_directories(clingo-dl_t PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
    target_compile_definitions(clingo-dl_t PUBLIC CLINGODL_NO_VISIBILITY)

    add_executable(test_c-api)
    clingodl_target_properties(TARGET test_c-api FOLDER "test" SUBDIR "test")
    target_sources(test_c-api PRIVATE ${test})
    target_link_libraries(test_c-api PRIVATE clingo-dl_t Catch2::Catch2WithMain)
    catch_discover_tests(test_c-api)
endif()
