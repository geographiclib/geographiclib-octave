project (GeographicLib-MATLAB NONE)

# Version information
set (PROJECT_VERSION_MAJOR 2)
set (PROJECT_VERSION_MINOR 0)
set (PROJECT_VERSION_PATCH 0)
set (PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
if (PROJECT_VERSION_PATCH GREATER 0)
  set (PROJECT_VERSION "${PROJECT_VERSION}.${PROJECT_VERSION_PATCH}")
endif ()
set (PROJECT_VERSION_SUFFIX "-alpha")
set (PROJECT_VERSION "${PROJECT_VERSION}${PROJECT_VERSION_SUFFIX}")

set (RELEASE_DATE "2022-mm-dd")

cmake_minimum_required (VERSION 3.13.0)

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

set (INSTDIR "share/octave/packages"
  CACHE PATH "Installation directory GeographicLib-Octave")
# directory for bundling up toolbox and packate
set (DISTDIR ${PROJECT_BINARY_DIR}/dist/geographiclib-octave-${PROJECT_VERSION})
file (REMOVE_RECURSE ${DISTDIR})

set (METAFILES COPYING INDEX NEWS)
set (PACKAGE_DEPENDS)
foreach (f ${METAFILES})
  configure_file (${f} ${DISTDIR}/${f} COPYONLY)
  set (PACKAGE_DEPENDS ${PACKAGE_DEPENDS} ${DISTDIR}/${f})
endforeach ()
configure_file (DESCRIPTION.in ${DISTDIR}/DESCRIPTION @ONLY)
set (PACKAGE_DEPENDS ${PACKAGE_DEPENDS} ${DISTDIR}/DESCRIPTION)
configure_file (DESCRIPTION.in DESCRIPTION @ONLY)
set (METAFILES ${METAFILES} DESCRIPTION)

enable_testing ()

add_subdirectory (inst)

if (EXISTS ${PROJECT_SOURCE_DIR}/.git AND NOT WIN32)
  add_custom_target (checktrailingspace
    COMMAND git ls-tree -r HEAD --name-only | xargs grep ' $$' || true
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Looking for trailing spaces")
  add_custom_target (checktabs
    COMMAND git ls-tree -r HEAD --name-only | grep -v makefile-admin |
    xargs grep '\t' || true
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Looking for tabs")
  add_custom_target (checkblanklines
    COMMAND git ls-tree -r HEAD --name-only |
    while read f\; do tr 'X\\n' 'YX' < $$f |
    egrep '\(^X|XXX|XX$$|[^X]$$\)' > /dev/null && echo $$f\; done || true
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Looking for extra blank lines")

  add_custom_target(sanitize)
  add_dependencies(sanitize checktrailingspace checktabs checkblanklines)
endif ()
