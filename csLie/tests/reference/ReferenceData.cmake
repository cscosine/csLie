# Runs the reference generator (generateReferenceData.cpp) and compares its output with the committed data/ folder
# (MODE CHECK, the test "csLie reference data up to date") or writes it there (MODE WRITE, the target
# csLieRegenerateReferenceData). The comparison (compareReferenceData.cpp) allows one unit in the last place per
# number: libquadmath's last bit differs between GCC versions.
#
#   cmake -DGENERATOR=<executable> -DCOMPARATOR=<executable> -DREFERENCE=<data folder> -DOUTPUT=<folder>
#     -DMODE=CHECK|WRITE -P ReferenceData.cmake
#
# OUTPUT: the folder the generator writes first (in the build folder), emptied before.

# A -P script does not inherit the policies of the project that runs it: without this, CMake 3.x leaves CMP0057 unset
# and does not treat IN_LIST as an operator (CMake 4 always does). Same version as the top-level CMakeLists.txt.
cmake_minimum_required(VERSION 3.31)

foreach(
  var
  GENERATOR
  COMPARATOR
  REFERENCE
  OUTPUT
  MODE
)
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
      COMMAND "${COMPARATOR}" "${REFERENCE}/${name}" "${OUTPUT}/${name}"
      RESULT_VARIABLE different
      ERROR_VARIABLE comparison
    )
    if(different)
      string(STRIP "${comparison}" comparison)
      list(APPEND differences "${name}: differs (${comparison})")
      # the first lines that differ, where diff exists (it does on the CI runners and on Unix)
      find_program(REFERENCE_DATA_DIFF diff)
      if(REFERENCE_DATA_DIFF)
        execute_process(
          COMMAND "${REFERENCE_DATA_DIFF}" "${REFERENCE}/${name}" "${OUTPUT}/${name}"
          OUTPUT_VARIABLE diffOutput
        )
        string(SUBSTRING "${diffOutput}" 0 3000 diffOutput)
        message(STATUS "${name}, committed (<) and generated (>), the first differences:\n${diffOutput}")
      endif()
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
