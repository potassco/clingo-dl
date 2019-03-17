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

if (CLINGODL_BUILD_SHARED)
    set(clingodl_lib_type SHARED)
else()
    set(clingodl_lib_type STATIC)
endif()

add_library(libclingo-dl ${clingodl_lib_type} ${header} ${source})
target_include_directories(libclingo-dl PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_compile_definitions(libclingo-dl PRIVATE CLINGODL_BUILD_LIBRARY)
if (NOT CLINGODL_BUILD_STATIC AND CLINGODL_BUILD_SHARED)
    set_target_properties(libclingo-dl PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(libclingo-dl PRIVATE libclingo)
if (NOT CLINGODL_BUILD_SHARED)
    target_compile_definitions(libclingo-dl PUBLIC CLINGODL_NO_VISIBILITY)
else()
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
set_target_properties(libclingo-dl PROPERTIES
    OUTPUT_NAME clingo-dl
    LIBRARY_OUTPUT_NAME clingo-dl
    FOLDER lib)

if (CLINGODL_BUILD_TESTS)
    add_subdirectory(tests)
endif()
