# Runs the reference generator (generateReferenceData.cpp) and compares its output with the committed data/ folder
# (MODE CHECK, the test "csLie reference data up to date") or writes it there (MODE WRITE, the target
# csLieRegenerateReferenceData).
#
#   cmake -DGENERATOR=<executable> -DREFERENCE=<data folder> -DOUTPUT=<folder> -DMODE=CHECK|WRITE
#     -P ReferenceData.cmake
#
# OUTPUT: the folder the generator writes first (in the build folder), emptied before.

foreach(var GENERATOR REFERENCE OUTPUT MODE)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "ReferenceData.cmake: ${var} is required")
  endif()
endforeach()

file(REMOVE_RECURSE "${OUTPUT}")
file(MAKE_DIRECTORY "${OUTPUT}")
execute_process(COMMAND "${GENERATOR}" "${OUTPUT}" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "ReferenceData.cmake: ${GENERATOR} failed (${result})")
endif()

file(GLOB generated RELATIVE "${OUTPUT}" "${OUTPUT}/*.hpp")
file(GLOB committed RELATIVE "${REFERENCE}" "${REFERENCE}/*.hpp")
list(SORT generated)
list(SORT committed)

if(MODE STREQUAL "CHECK")
  set(differences "")
  foreach(name IN LISTS generated)
    if(NOT name IN_LIST committed)
      list(APPEND differences "${name}: missing")
      continue()
    endif()
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E compare_files "${OUTPUT}/${name}" "${REFERENCE}/${name}"
      RESULT_VARIABLE different
    )
    if(different)
      list(APPEND differences "${name}: differs")
    endif()
  endforeach()
  foreach(name IN LISTS committed)
    if(NOT name IN_LIST generated)
      list(APPEND differences "${name}: not generated any more")
    endif()
  endforeach()
  if(differences)
    list(JOIN differences "\n  " differences)
    message(
      FATAL_ERROR
      "${REFERENCE} is not the output of the generator (${OUTPUT}): edited by hand, or the generator changed without "
      "regenerating it.\n  ${differences}\nRegenerate it with:\n"
      "  cmake --build <build dir> --config <config> --target csLieRegenerateReferenceData\n"
      "then review and commit the diff."
    )
  endif()
  message(STATUS "${REFERENCE} is up to date")
elseif(MODE STREQUAL "WRITE")
  file(MAKE_DIRECTORY "${REFERENCE}")
  foreach(name IN LISTS generated)
    file(COPY_FILE "${OUTPUT}/${name}" "${REFERENCE}/${name}" ONLY_IF_DIFFERENT)
  endforeach()
  foreach(name IN LISTS committed)
    if(NOT name IN_LIST generated)
      file(REMOVE "${REFERENCE}/${name}")
    endif()
  endforeach()
  message(STATUS "${REFERENCE} regenerated: review and commit the diff")
else()
  message(FATAL_ERROR "ReferenceData.cmake: MODE is CHECK or WRITE, not ${MODE}")
endif()
