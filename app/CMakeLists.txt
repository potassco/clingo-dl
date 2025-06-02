# [[[source: src
set(ide_source_group "Source Files")
set(source-group
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc")
source_group("${ide_source_group}" FILES ${source-group})
set(source
    ${source-group})
# ]]]

add_executable(clingo-dl-app)
target_sources(clingo-dl-app PRIVATE ${source})
clingodl_target_properties(TARGET clingo-dl-app FOLDER "bin" TYPE default)
target_link_libraries(clingo-dl-app PRIVATE clingo-dl-optimizer)
set_target_properties(clingo-dl-app PROPERTIES OUTPUT_NAME clingo-dl)
