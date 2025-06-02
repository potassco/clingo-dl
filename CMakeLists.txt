cmake_minimum_required(VERSION 3.22.1)

file(STRINGS "lib/c-api/include/clingo-dl.h" CLINGODL_VERSION_LINE REGEX "#define CLINGODL_VERSION ")
string(REGEX REPLACE "#define CLINGODL_VERSION[ \t]+\"([^\"]+)\"" "\\1" CLINGODL_VERSION "${CLINGODL_VERSION_LINE}")

project(CLINGODL VERSION "${CLINGODL_VERSION}" LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS ON)

include(GNUInstallDirs)
include(CMakeDependentOption)
include(macros)

option(CLINGODL_BUILD_APP "build the app" ON)
option(CLINGODL_BUILD_TESTS "build tests" OFF)
option(CLINGODL_BUILD_SHARED "build clingo-dl library shared" ON)
option(CLINGODL_CROSSCHECK "enable additional state checking" OFF)
option(CLINGODL_USE_LOCAL_CATCH "Advanced option to build against bundled or installed catch." ON)
option(CLINGODL_PROFILE "profile clingo-dl" OFF)
option(CLINGODL_MANAGE_RPATH "manage rpath when installing" ON)

cmake_dependent_option(CLINGODL_BUILD_SHARED "build clingo as shared library" ON "NOT EMSCRIPTEN" OFF)
cmake_dependent_option(CLINGODL_USE_VERSION_SCRIPT "use a version script to hide C++ symbols" ON "CLINGODL_BUILD_SHARED" OFF)
cmake_dependent_option(CLINGODL_USE_VISIBILITY "enable visibility attributes" ON "CLINGODL_BUILD_SHARED" OFF)

set(PYCLINGODL_INSTALL "prefix" CACHE STRING "where to install python module: user, system, prefix")
set(PYCLINGODL_INSTALL_DIR "" CACHE STRING "path where to install python module")

if (MSVC)
    set(VC_RELEASE_LINK_OPTIONS /LTCG)
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE    "${CMAKE_EXE_LINKER_FLAGS_RELEASE}    ${VC_RELEASE_LINK_OPTIONS}")
    SET(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS_RELEASE} ${VC_RELEASE_LINK_OPTIONS}")
    SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} ${VC_RELEASE_LINK_OPTIONS}")
    SET(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} ${VC_RELEASE_LINK_OPTIONS}")
endif()

if(${CMAKE_VERSION} VERSION_LESS "3.26")
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

if (CLINGODL_USE_VISIBILITY)
    set(CMAKE_C_VISIBILITY_PRESET hidden)
    set(CMAKE_CXX_VISIBILITY_PRESET hidden)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
endif()

set(PYBIND11_NEWPYTHON ON)
set(PYBIND11_FINDPYTHON ON)
set(CLINGODL_INSTALL_DEFAULT OFF)
set(CLINGODL_INSTALL_EXTRA OFF)
set(CLINGODL_INSTALL_WHEEL OFF)

if (DEFINED SKBUILD)
    message(STATUS "Configuring for skbuild")
    set(CLINGODL_INSTALL_WHEEL ON)
    set(CLINGODL_BUILD_TESTS OFF)
    set(CLINGODL_BUILD_APP OFF)
    find_package(pybind11 CONFIG REQUIRED)

    add_library(_clingo INTERFACE)
    if (PYCLINGODL_PIP_PATH)
        set(pyclingodl_path_ "${PYCLINGODL_PIP_PATH}")
    else()
        execute_process(COMMAND ${Python_EXECUTABLE} ${CLINGODL_SOURCE_DIR}/cmake/python-site.py target OUTPUT_VARIABLE pyclingodl_path_ OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
    file(TO_CMAKE_PATH "${pyclingodl_path_}" pyclingodl_path)
    target_include_directories(_clingo INTERFACE "${pyclingodl_path}/clingo-include")
    target_link_libraries(_clingo INTERFACE clingo)
    target_link_directories(_clingo INTERFACE "${pyclingodl_path}")
    add_library(clingo::clingo ALIAS _clingo)
else()
    if(CLINGODL_MANAGE_RPATH)
        list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
        if(isSystemDir EQUAL -1)
            set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
        endif()
    endif()
    set(CLINGODL_INSTALL_DEFAULT ON)
    if (NOT CLINGODL_BUILD_SHARED)
        set(CLINGODL_INSTALL_EXTRA ON)
    endif()
    find_package(Python COMPONENTS Interpreter Development Development.Module)
    find_package(clingo 6.0.0 REQUIRED)
endif()

if (CLINGODL_BUILD_TESTS)
    enable_testing()
endif()

if (CLINGODL_PROFILE)
    find_package(Gperftools REQUIRED)
endif()

add_subdirectory(third_party)
add_subdirectory(lib libclingo-dl)
if (CLINGODL_BUILD_APP)
    add_subdirectory(app)
endif()

if (CLINGODL_INSTALL_DEFAULT)
    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/clingo-dl-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/clingo-dl-config.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/clingo-dl"
    )
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/clingo-dl-config-version.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
    )
    install(
        FILES
            "${CMAKE_CURRENT_BINARY_DIR}/clingo-dl-config.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/clingo-dl-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/clingo-dl"
    )
    install(EXPORT clingo-dl-targets
        FILE clingo-dl-targets.cmake
        NAMESPACE clingo-dl::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/clingo-dl"
    )
endif()
