# Runs the reference generator (generateReferenceData.cpp) and compares its output with the committed
# ReferenceData.hpp (MODE CHECK, the test "csLie reference data up to date") or writes it there (MODE WRITE, the target
# csLieRegenerateReferenceData).
#
#   cmake -DGENERATOR=<executable> -DREFERENCE=<ReferenceData.hpp> -DOUTPUT=<file> -DMODE=CHECK|WRITE
#     -P ReferenceData.cmake
#
# OUTPUT: where the output of the generator is written first (in the build folder).

foreach(var GENERATOR REFERENCE OUTPUT MODE)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "ReferenceData.cmake: ${var} is required")
  endif()
endforeach()

get_filename_component(output_dir "${OUTPUT}" DIRECTORY)
file(MAKE_DIRECTORY "${output_dir}")
execute_process(COMMAND "${GENERATOR}" OUTPUT_FILE "${OUTPUT}" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "ReferenceData.cmake: ${GENERATOR} failed (${result})")
endif()

if(MODE STREQUAL "CHECK")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${OUTPUT}" "${REFERENCE}"
    RESULT_VARIABLE different
  )
  if(different)
    message(
      FATAL_ERROR
      "${REFERENCE} differs from the output of the generator (${OUTPUT}): it was edited by hand, or the generator "
      "changed without regenerating it. Regenerate it with:\n"
      "  cmake --build <build dir> --config <config> --target csLieRegenerateReferenceData\n"
      "then review and commit the diff."
    )
  endif()
  message(STATUS "${REFERENCE} is up to date")
elseif(MODE STREQUAL "WRITE")
  file(COPY_FILE "${OUTPUT}" "${REFERENCE}" ONLY_IF_DIFFERENT)
  message(STATUS "${REFERENCE} regenerated: review and commit the diff")
else()
  message(FATAL_ERROR "ReferenceData.cmake: MODE is CHECK or WRITE, not ${MODE}")
endif()
