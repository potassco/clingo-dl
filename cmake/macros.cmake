find_package(Git QUIET)

function(clingodl_target_properties)
    set(options)
    set(single_values TARGET FOLDER TYPE SUBDIR)
    set(multi_values)
    cmake_parse_arguments(clingodl "${options}" "${single_values}" "${multi_values}" ${ARGV})

    set(binary_subdir "bin")
    set(library_subdir "lib")
    if(clingodl_SUBDIR)
        set(binary_subdir "bin/${clingodl_SUBDIR}")
        set(library_subdir "lib/${clingodl_SUBDIR}")
    endif()

    get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if(is_multi_config)
        set(binary_subdir "${binary_subdir}/$<CONFIG>")
        set(library_subdir "${library_subdir}/$<CONFIG>")
    endif()

    if (clingodl_FOLDER)
        set_target_properties(${clingodl_TARGET} PROPERTIES
            FOLDER "${clingodl_FOLDER}"
            POSITION_INDEPENDENT_CODE ON
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${binary_subdir}"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${binary_subdir}"
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${library_subdir}"
            PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${binary_subdir}")
    endif()

    if(clingodl_TYPE STREQUAL "extra" AND CLINGODL_INSTALL_EXTRA)
        install(
            TARGETS "${clingodl_TARGET}"
            EXPORT clingo-dl-targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
    elseif((clingodl_TYPE STREQUAL "default" OR clingodl_TYPE STREQUAL "binary") AND CLINGODL_INSTALL_DEFAULT)
        install(
            TARGETS "${clingodl_TARGET}"
            EXPORT clingo-dl-targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
    endif()
endfunction()

function(potassco_git_hash name)
    string(TOLOWER "${name}" name_lower)
    string(TOUPPER "${name}" name_upper)

    set(target "${name}-hash")
    add_library(${target} INTERFACE)

    set(header_directory
        "${CMAKE_CURRENT_BINARY_DIR}/generated"
    )
    set(header_file
        "${header_directory}/${name_lower}_git_hash.h"
    )
    set(git_hash_file
        "${PROJECT_SOURCE_DIR}/.git_hash"
    )
    set(git_hash_var
        "${name_upper}_GIT_HASH"
    )

    file(MAKE_DIRECTORY "${header_directory}")
    if(DEFINED ${git_hash_var} AND NOT "${${git_hash_var}}" STREQUAL "")
        file(WRITE "${header_file}"
"#pragma once

#define ${git_hash_var} \"${${git_hash_var}}\"
")
    elseif(EXISTS "${git_hash_file}")
        file(READ "${git_hash_file}" git_hash_from_file)
        string(STRIP "${git_hash_from_file}" git_hash_from_file)
        if(NOT "${git_hash_from_file}" STREQUAL "")
            file(WRITE "${header_file}"
"#pragma once

#define ${git_hash_var} \"${git_hash_from_file}\"
")
        else()
            file(WRITE "${header_file}" "#pragma once\n")
        endif()
    elseif(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        add_custom_target("${target}_git_hash" ALL
            COMMAND "${CMAKE_COMMAND}"
                "-DGIT_EXECUTABLE=${GIT_EXECUTABLE}"
                "-DSOURCE_DIR=${PROJECT_SOURCE_DIR}"
                "-DOUTPUT_FILE=${header_file}"
                "-DNAME=${git_hash_var}"
                -P "${PROJECT_SOURCE_DIR}/cmake/write-git-hash.cmake"
            BYPRODUCTS
                "${header_file}"
            VERBATIM)

        add_dependencies("${target}" "${target}_git_hash")

    else()
        file(WRITE "${header_file}" "#pragma once\n")
    endif()

    target_sources("${target}" INTERFACE "${header_file}")
    target_include_directories("${target}" INTERFACE "${header_directory}")
endfunction()
