set(KLFT_ANALYSIS_SCRIPTS
  plaquette_stats.py
  polyakov_stats.py
  stats_common.py
  wtemp_stats.py)

set(KLFT_ANALYSIS_OUTPUT_DIR ${PROJECT_BINARY_DIR}/analysis)
file(MAKE_DIRECTORY ${KLFT_ANALYSIS_OUTPUT_DIR})

foreach(script IN LISTS KLFT_ANALYSIS_SCRIPTS)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${script}
    ${KLFT_ANALYSIS_OUTPUT_DIR}/${script}
    COPYONLY)

  add_custom_command(
    OUTPUT ${KLFT_ANALYSIS_OUTPUT_DIR}/${script}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/${script}
            ${KLFT_ANALYSIS_OUTPUT_DIR}/${script}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${script}
    COMMENT "Copying analysis/${script}"
    VERBATIM)
  list(APPEND KLFT_ANALYSIS_OUTPUTS ${KLFT_ANALYSIS_OUTPUT_DIR}/${script})
endforeach()

add_custom_target(analysis_scripts ALL
  DEPENDS ${KLFT_ANALYSIS_OUTPUTS})
