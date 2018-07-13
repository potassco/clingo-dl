cmake_minimum_required(VERSION 3.1)
file(READ "clingoDL.h" main)
string(REGEX MATCH "#define CLINGODL_VERSION \"([^\"]*)\"" clingov ${main})
project(CLINGODL VERSION "${CMAKE_MATCH_1}" LANGUAGES CXX)
# Enable folders in IDEs like Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type selected - using 'Release'")
    set(CMAKE_BUILD_TYPE "Release")
endif()

include(GNUInstallDirs)

option(CLINGODL_BUILD_WITH_SYSTEM_CLINGO "if true the local clingo copy will not be used" OFF)

if (NOT MSVC)
    if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    endif()
    if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    endif()
    if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    endif()
else()
    set(VC_RELEASE_LINK_OPTIONS /LTCG)
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE    "${CMAKE_EXE_LINKER_FLAGS_RELEASE}    ${VC_RELEASE_LINK_OPTIONS}")
    SET(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS_RELEASE} ${VC_RELEASE_LINK_OPTIONS}")
    SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} ${VC_RELEASE_LINK_OPTIONS}")
    SET(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} ${VC_RELEASE_LINK_OPTIONS}")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/clingo/CMakeLists.txt OR CLINGODL_BUILD_WITH_SYSTEM_CLINGO)
    find_package(Clingo REQUIRED)
    message(STATUS "Local clingo copy not available - using system installation")
else()
    message(STATUS "Clingo is not installed - using local copy")
    add_subdirectory(clingo)
endif()

# TODO: folder properties for main.cpp
#       (in clingo/cmake there is a convienient script to generate those)
set(appfiles main.cpp propagator.cpp clingoDL.cpp)
set(libfiles propagator.cpp clingoDL.cpp)

add_executable(clingoDL ${appfiles})
target_link_libraries(clingoDL PRIVATE libclingo)
set_target_properties(clingoDL PROPERTIES FOLDER exe)
#install(TARGETS clingoDL RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
# LIBRARY DESTINATION DOES NOT SEEM TO WORK
add_library(libclingoDL SHARED ${libfiles})
target_compile_definitions(libclingoDL PRIVATE CLINGODL_BUILD_LIBRARY)
target_link_libraries(libclingoDL PRIVATE libclingo)
set_target_properties(libclingoDL PROPERTIES PUBLIC_HEADER clingoDL.h)
set_target_properties(libclingoDL PROPERTIES LIBRARY_OUTPUT_NAME clingoDL)

install(TARGETS clingoDL libclingoDL 
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
