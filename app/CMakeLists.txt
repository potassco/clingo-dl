# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]

add_executable(clingo-dl-bin)
target_sources(clingo-dl-bin PRIVATE ${source})
clingodl_target_properties(TARGET clingo-dl-bin FOLDER "bin" TYPE default)
target_link_libraries(clingo-dl-bin PRIVATE clingo-dl-optimizer)
set_target_properties(clingo-dl-bin PROPERTIES OUTPUT_NAME clingo-dl)
