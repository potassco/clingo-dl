# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/clingodl.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]

pybind11_add_module(clingo-dl-python-api)
target_sources(clingo-dl-python-api PRIVATE ${source})
clingodl_target_properties(TARGET clingo-dl-python-api FOLDER "lib" SUBDIR "python")
set_target_properties(clingo-dl-python-api PROPERTIES OUTPUT_NAME "clingodl")
target_link_libraries(clingo-dl-python-api PRIVATE clingo-dl Python::Module)

if (CLINGODL_INSTALL_WHEEL)
    install(TARGETS clingo-dl-python-api
        RUNTIME DESTINATION "."
        LIBRARY DESTINATION ".")
    if (CLINGODL_BUILD_SHARED)
        if (APPLE)
            set_target_properties(clingo-dl-python-api
                PROPERTIES INSTALL_RPATH "@loader_path"
                BUILD_WITH_INSTALL_RPATH TRUE)
        elseif (UNIX)
            set_target_properties(clingo-dl-python-api
                PROPERTIES INSTALL_RPATH "$ORIGIN"
                BUILD_WITH_INSTALL_RPATH TRUE)
        endif()
    endif()
elseif(CLINGODL_INSTALL_DEFAULT)
    if (PYCLINGODL_INSTALL_DIR STREQUAL "")
        if (PYCLINGODL_INSTALL STREQUAL "user")
            execute_process(
                COMMAND "${Python_EXECUTABLE}" "${CLINGODL_SOURCE_DIR}/cmake/python-site.py" target --user
                OUTPUT_VARIABLE _python_install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)
        elseif(PYCLINGODL_INSTALL STREQUAL "prefix")
            execute_process(
                COMMAND "${Python_EXECUTABLE}" "${CLINGODL_SOURCE_DIR}/cmake/python-site.py" target --prefix "${CMAKE_INSTALL_PREFIX}"
                OUTPUT_VARIABLE _python_install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)
        else()
            execute_process(
                COMMAND "${Python_EXECUTABLE}" "${CLINGODL_SOURCE_DIR}/cmake/python-site.py" target
                OUTPUT_VARIABLE _python_install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()
        file(TO_CMAKE_PATH "${_python_install_dir}" python_install_dir)
    else()
        set(python_install_dir "${PYCLINGODL_INSTALL_DIR}")
    endif()
    MESSAGE(STATUS "Python module path: ${python_install_dir}")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/INSTALLER" "pip\n")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/METADATA" "Metadata-Version: 1.1\nName: clingo\nVersion: ${CLINGODL_VERSION}\n")
    install(
        FILES "${CMAKE_CURRENT_BINARY_DIR}/INSTALLER" "${CMAKE_CURRENT_BINARY_DIR}/METADATA"
        DESTINATION "${python_install_dir}/clingodl-${CLINGODL_VERSION}.dist-info")
    install(
        TARGETS clingo-dl-python-api
        RUNTIME DESTINATION "${python_install_dir}"
        LIBRARY DESTINATION "${python_install_dir}")
endif()

if (CLINGODL_BUILD_TESTS AND PYTEST_EXECUTABLE)
    if (WIN32)
        # NOTE: python from conda-forge with more than one element in the path fails for whatever reason
        set(python_subdir "${CMAKE_BINARY_DIR}/bin/python/")
        get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
        if(is_multi_config)
            set(python_subdir "${python_subdir}$<CONFIG>/")
        endif()
        add_custom_command(
            TARGET clingo-dl-python-api POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory "${python_subdir}"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:clingo-dl-python-api>" "${python_subdir}"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:clingo::pyclingo>" "${python_subdir}"
        )
        pytest_discover_tests(
            test_pyclingo
            PYTHON_PATH_PREPEND "${python_subdir}"
            WORKING_DIRECTORY "tests"
            DEPENDS clingo-dl-python-api
        )
    else()
        pytest_discover_tests(
            test_pyclingo
            PYTHON_PATH_PREPEND "$<TARGET_FILE_DIR:clingo-dl-python-api>" "$<TARGET_FILE_DIR:clingo::pyclingo>"
            LIBRARY_PATH_PREPEND "$<TARGET_FILE_DIR:clingo-dl>" "$<TARGET_FILE_DIR:clingo::clingo>"
            WORKING_DIRECTORY "tests"
            DEPENDS clingo-dl-python-api
        )
    endif()
endif()
