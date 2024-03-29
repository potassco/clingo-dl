# [[[header: .
set(ide_header_group "Header Files")
set(header-group
    "${CMAKE_CURRENT_SOURCE_DIR}/clingo-dl.h")
source_group("${ide_header_group}" FILES ${header-group})
set(header-group-clingo-dl
    "${CMAKE_CURRENT_SOURCE_DIR}/clingo-dl/propagator.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/clingo-dl/util.hh")
source_group("${ide_header_group}\\clingo-dl" FILES ${header-group-clingo-dl})
set(header
    ${header-group}
    ${header-group-clingo-dl})
# ]]]
# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/clingo-dl.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/propagator.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]

if (PYCLINGODL_ENABLE STREQUAL "pip")
    set(clingodl_build_interface ON)
else()
    set(clingodl_build_interface OFF)
endif()

set(clingodl_public_scope_ PUBLIC)
set(clingodl_private_scope_ PRIVATE)
if (CLINGODL_BUILD_SHARED)
    set(clingodl_lib_type SHARED)
    if (clingodl_build_interface)
        set(clingodl_lib_type INTERFACE)
        set(clingodl_public_scope_ INTERFACE)
        set(clingodl_private_scope_ INTERFACE)
    endif()
else()
    set(clingodl_lib_type STATIC)
endif()

add_library(libclingo-dl ${clingodl_lib_type})
target_sources(libclingo-dl ${clingodl_private_scope_} ${header} ${source})
target_include_directories(libclingo-dl ${clingodl_public_scope_} "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_compile_definitions(libclingo-dl ${clingodl_private_scope_} CLINGODL_BUILD_LIBRARY)

if (NOT CLINGODL_BUILD_STATIC AND CLINGODL_BUILD_SHARED AND NOT clingodl_build_interface)
    set_target_properties(libclingo-dl PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(libclingo-dl ${clingodl_public_scope_} libclingo)
if (NOT CLINGODL_BUILD_SHARED)
    target_compile_definitions(libclingo-dl ${clingodl_public_scope_}  CLINGODL_NO_VISIBILITY)
elseif(NOT clingodl_build_interface)
    set_target_properties(libclingo-dl PROPERTIES PUBLIC_HEADER clingo-dl.h)
    install(TARGETS libclingo-dl
        EXPORT ClingoDLTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(EXPORT ClingoDLTargets
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ClingoDL)

    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        ${PROJECT_SOURCE_DIR}/cmake/ClingoDLConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/ClingoDLConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ClingoDL)
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/ClingoDLConfigVersion.cmake
        COMPATIBILITY SameMajorVersion)
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/ClingoDLConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/ClingoDLConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ClingoDL)
endif()
if (NOT clingodl_build_interface)
    set_target_properties(libclingo-dl PROPERTIES
        OUTPUT_NAME clingo-dl
        LIBRARY_OUTPUT_NAME clingo-dl
        FOLDER lib)
endif()

if (CLINGODL_BUILD_TESTS)
    # NOTE: we build statically without hiding symbols
    add_library(libclingo-dl_t STATIC ${header} ${source})
    target_include_directories(libclingo-dl_t PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
    target_compile_definitions(libclingo-dl_t PUBLIC CLINGODL_NO_VISIBILITY)
    target_link_libraries(libclingo-dl_t PUBLIC libclingo)
    set_target_properties(libclingo-dl_t PROPERTIES
        OUTPUT_NAME clingo-dl_t
        LIBRARY_OUTPUT_NAME clingo-dl_t
        FOLDER lib)
    add_subdirectory(tests)
endif()
