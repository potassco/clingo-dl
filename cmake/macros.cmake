function(clingodl_target_properties)
    set(options)
    set(single_values TARGET FOLDER TYPE SUBDIR)
    set(multi_values)
    cmake_parse_arguments(clingodl "${options}" "${single_values}" "${multi_values}" ${ARGV})

    set(binary_subdir "bin")
    if(clingodl_SUBDIR)
        set(binary_subdir "bin/${clingodl_SUBDIR}")
    endif()

    get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    set(base_dir "${CMAKE_BINARY_DIR}")
    if(is_multi_config)
        set(base_dir "${base_dir}/$<CONFIG>")
    endif()

    set_target_properties("${clingodl_TARGET}" PROPERTIES
        FOLDER "${clingodl_FOLDER}"
        POSITION_INDEPENDENT_CODE ON
        RUNTIME_OUTPUT_DIRECTORY "${base_dir}/${binary_subdir}"
        LIBRARY_OUTPUT_DIRECTORY "${base_dir}/${binary_subdir}"
        ARCHIVE_OUTPUT_DIRECTORY "${base_dir}/lib"
        PDB_OUTPUT_DIRECTORY "${base_dir}/bin"
    )

    if(clingodl_TYPE STREQUAL "extra" AND CLINGODL_INSTALL_EXTRA)
        install(
            TARGETS "${clingodl_TARGET}"
            EXPORT clingo-targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
    elseif((clingodl_TYPE STREQUAL "default" OR clingodl_TYPE STREQUAL "binary") AND CLINGODL_INSTALL_DEFAULT)
        install(
            TARGETS "${clingodl_TARGET}"
            EXPORT clingo-targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
    endif()
endfunction()
