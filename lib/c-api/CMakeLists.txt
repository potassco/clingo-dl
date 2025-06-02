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
target_link_libraries(clingo-dl PUBLIC clingo::clingo PRIVATE clingo-dl-solver)
if (CLINGODL_USE_VISIBILITY)
    target_compile_definitions(clingo-dl PRIVATE "$<BUILD_INTERFACE:CLINGODL_BUILD_LIBRARY>")
else()
    target_compile_definitions(clingo-dl PUBLIC CLINGODL_NO_VISIBILITY)
endif()

if (CLINGODL_INSTALL_DEFAULT)
    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
elseif (CLINGODL_INSTALL_WHEEL AND CLINGODL_BUILD_SHARED)
    install(
        TARGETS clingo-dl
        RUNTIME DESTINATION "."
        LIBRARY DESTINATION ".")
    if (CLINGODL_BUILD_SHARED)
        if (APPLE)
            set_target_properties(clingo-dl
                PROPERTIES INSTALL_RPATH "@loader_path"
                BUILD_WITH_INSTALL_RPATH TRUE)
        elseif (UNIX)
            set_target_properties(clingo-dl
                PROPERTIES INSTALL_RPATH "$ORIGIN"
                BUILD_WITH_INSTALL_RPATH TRUE)
        endif()
    endif()
endif()

if (CLINGODL_BUILD_TESTS)
    if (WIN32)
        # NOTE: needed to get the test explorer working
        set(test_subdir "${CMAKE_BINARY_DIR}/bin/tests/")
        set(python_subdir "${CMAKE_BINARY_DIR}/bin/python/")
        get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
        if(is_multi_config)
            set(test_subdir "${test_subdir}$<CONFIG>/")
            set(python_subdir "${python_subdir}$<CONFIG>/")
        endif()
        if (CLINGODL_BUILD_SHARED)
            add_custom_command(
                TARGET clingo-dl POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory "${python_subdir}"
                COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:clingo-dl>" "${python_subdir}"
            )
        endif()
        add_custom_command(
            TARGET clingo-dl POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory "${test_subdir}"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${python_subdir}"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:clingo::clingo>" "${test_subdir}"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:clingo::clingo>" "${python_subdir}"
        )
    endif()

    add_library(clingo-dl_t STATIC)
    clingodl_target_properties(TARGET clingo-dl_t FOLDER "lib")
    target_sources(clingo-dl_t PRIVATE ${header} ${source})
    target_link_libraries(clingo-dl_t PUBLIC clingo-dl-solver_t)
    target_include_directories(clingo-dl_t PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
    target_compile_definitions(clingo-dl_t PUBLIC CLINGODL_NO_VISIBILITY)

    add_executable(test_clingo-dl)
    clingodl_target_properties(TARGET test_clingo-dl FOLDER "test" SUBDIR "tests")
    target_sources(test_clingo-dl PRIVATE ${test})
    target_link_libraries(test_clingo-dl PRIVATE clingo-dl_t Catch2::Catch2WithMain)
    catch_discover_tests(test_clingo-dl DL_PATH "$<TARGET_FILE_DIR:clingo::clingo>")
endif()
