# [[[source: .
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/_clingodl.c")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]

if (NOT PYCLINGODL_INSTALL_DIR AND Python_EXECUTABLE)
    if (PYCLINGODL_INSTALL STREQUAL "user")
        execute_process(COMMAND ${Python_EXECUTABLE} ${CLINGODL_SOURCE_DIR}/cmake/python-site.py target --user OUTPUT_VARIABLE PYCLINGODL_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    elseif(PYCLINGODL_INSTALL STREQUAL "prefix")
        execute_process(COMMAND ${Python_EXECUTABLE} ${CLINGODL_SOURCE_DIR}/cmake/python-site.py target --prefix "${CMAKE_INSTALL_PREFIX}" OUTPUT_VARIABLE PYCLINGODL_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
        execute_process(COMMAND ${Python_EXECUTABLE} ${CLINGODL_SOURCE_DIR}/cmake/python-site.py target OUTPUT_VARIABLE PYCLINGODL_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
endif()

if (NOT PYCLINGODL_SUFFIX)
    if (Python_EXECUTABLE)
        execute_process(COMMAND ${Python_EXECUTABLE} ${CLINGODL_SOURCE_DIR}/cmake/python-site.py suffix OUTPUT_VARIABLE PYCLINGODL_SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
        if (CYGWIN)
            set(PYCLINGODL_SUFFIX ".dll")
        elseif (UNIX)
            set(PYCLINGODL_SUFFIX ".so")
        else()
            set(PYCLINGODL_SUFFIX ".pyd")
        endif()
    endif()
endif()

# Note: For the pip package, we actually not just want to load this library but
# also link agaist it. This is mainly for Windows to ensure that the import
# library is installed.
if (PYCLINGODL_ENABLE STREQUAL "pip")
    add_library(pyclingo-dl SHARED ${source})
else()
    add_library(pyclingo-dl MODULE ${source})
endif()
target_link_libraries(pyclingo-dl PRIVATE libclingo-dl Python::Module)
target_compile_definitions(pyclingo-dl PRIVATE _CFFI_NO_LIMITED_API)
set_target_properties(pyclingo-dl PROPERTIES
    FOLDER lib
    OUTPUT_NAME _clingodl
    SUFFIX ${PYCLINGODL_SUFFIX}
    PREFIX "")
if (CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set_target_properties(pyclingo-dl PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/python)
endif()
if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set_target_properties(pyclingo-dl PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/python)
endif()
if (CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set_target_properties(pyclingo-dl PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/python)
endif()

if (PYCLINGODL_INSTALL_DIR)
    file(TO_CMAKE_PATH "${PYCLINGODL_INSTALL_DIR}" _PYCLINGODL_INSTALL_DIR)
    install(TARGETS pyclingo-dl
        RUNTIME DESTINATION ${_PYCLINGODL_INSTALL_DIR}/clingodl
        LIBRARY DESTINATION ${_PYCLINGODL_INSTALL_DIR}/clingodl
        ARCHIVE DESTINATION ${_PYCLINGODL_INSTALL_DIR}/clingodl)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/libpyclingo-dl/clingodl/
        DESTINATION ${_PYCLINGODL_INSTALL_DIR}/clingodl
        FILES_MATCHING
        PATTERN "*.py"
        PATTERN "*.typed"
        PATTERN "tests" EXCLUDE)
    if (NOT PYCLINGODL_ENABLE STREQUAL "pip")
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/INSTALLER "pip\n")
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/METADATA "Metadata-Version: 1.1\nName: clingo-dl\nVersion: ${CLINGODL_VERSION}\n")
        install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/INSTALLER
            ${CMAKE_CURRENT_BINARY_DIR}/METADATA
            DESTINATION ${_PYCLINGODL_INSTALL_DIR}/clingodl-${CLINGODL_VERSION}.dist-info)
    else()
        install(FILES ${CMAKE_SOURCE_DIR}/libclingo-dl/clingo-dl.h
            DESTINATION ${_PYCLINGODL_INSTALL_DIR}/clingodl)
    endif()
endif()
