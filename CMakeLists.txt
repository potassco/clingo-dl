cmake_minimum_required(VERSION 3.1)
file(READ "main.cpp" main)
string(REGEX MATCH "#define CLINGODL_VERSION \"([^\"]*)\"" clingov ${main})
project(CLINGODL VERSION "${CMAKE_MATCH_1}" LANGUAGES CXX)
# Enable folders in IDEs like Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
#if (POLICY CMP0063)
#	cmake_policy(SET CMP0063 NEW)
#endif()
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "No build type selected - using 'Release'")
	set(CMAKE_BUILD_TYPE "Release")
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#include(GNUInstallDirs)

# Configuration options
option(CLINGODL_BUILD_APP          "whether or not to build the clingoDL application" ON)
#option(CLINGODL_BUILD_STATIC       "whether or not to link statically (if supported)" OFF)
#option(CLINGODL_BUILD_TESTS        "whether or not to build clingoDL unit tests" OFF)
#option(CLINGODL_BUILD_EXAMPLES     "whether or not to build examples" OFF)
option(CLINGODL_BUILD_WITH_THREADS "whether or not to build clingoDL with threading support (requires C++11)" ON)
#option(CLINGODL_INSTALL_LIB        "whether or not to install libclingoDL" OFF)

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
#	if (CLINGODL_BUILD_STATIC)
#		# force static runtime
#		string(REGEX REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
#	endif()
endif()

set(clingoDL_include_dest "clingoDL-${CLINGODL_VERSION}")
set(clingoDL_library_dest "clingoDL-${CLINGODL_VERSION}")
set(cmake_dest            "clingoDL-${CLINGODL_VERSION}/cmake")

if (CLINGODL_INSTALL_LIB AND NOT CMAKE_INSTALL_LIBDIR)
	message(STATUS "LIBDIR no set - using lib")
	set(CMAKE_INSTALL_LIBDIR lib)
endif()


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CLASP_BUILD_WITH_THREADS CLINGODL_BUILD_WITH_THREADS CACHE BOOL "")

find_package(clingo)

set(files main.cpp)

# Build clingoDL
if(CLINGODL_BUILD_APP)
	add_executable(clingoDL ${files})
	target_link_libraries(clingoDL PRIVATE Clingo::Clingo)
	set_target_properties(clingoDL PROPERTIES FOLDER exe)
endif()



