# [[[header: clingo-dl-app
set(ide_header_group "Header Files")
set(header-group-clingo-dl-app
    "${CMAKE_CURRENT_SOURCE_DIR}/clingo-dl-app/app.hh")
source_group("${ide_header_group}\\clingo-dl-app" FILES ${header-group-clingo-dl-app})
set(header
    ${header-group-clingo-dl-app})
# ]]]
# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/app.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]

add_library(libclingo-dl-app STATIC)
target_sources(libclingo-dl-app PRIVATE ${header} ${source})
target_include_directories(libclingo-dl-app PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")

if (NOT CLINGODL_BUILD_STATIC AND CLINGODL_BUILD_SHARED AND NOT clingodl_build_interface)
    set_target_properties(libclingo-dl-app PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(libclingo-dl-app PUBLIC libclingo-dl)

set_target_properties(libclingo-dl-app PROPERTIES
    OUTPUT_NAME clingo-dl-app
    LIBRARY_OUTPUT_NAME clingo-dl-app
    FOLDER lib)

if (CLINGODL_BUILD_TESTS)
    add_library(libclingo-dl-app_t STATIC ${header} ${source})
    target_include_directories(libclingo-dl-app_t PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
    target_link_libraries(libclingo-dl-app_t PUBLIC libclingo-dl_t)
    set_target_properties(libclingo-dl-app_t PROPERTIES
        OUTPUT_NAME clingo-dl-app_t
        LIBRARY_OUTPUT_NAME clingo-dl-app_t
        FOLDER lib)
endif()
