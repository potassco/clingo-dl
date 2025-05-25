# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/clingodl.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]

pybind11_add_module(pyclingodl ${source})
set_target_properties(pyclingodl PROPERTIES OUTPUT_NAME "clingodl")
target_link_libraries(pyclingodl PRIVATE libclingo-dl Python::Module)
