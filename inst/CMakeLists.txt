set (SOURCES
  geoddistance.m
  geodreckon.m
  geodarea.m
  tranmerc_fwd.m
  tranmerc_inv.m
  polarst_fwd.m
  polarst_inv.m
  eqdazim_fwd.m
  eqdazim_inv.m
  cassini_fwd.m
  cassini_inv.m
  gnomonic_fwd.m
  gnomonic_inv.m
  utmups_fwd.m
  utmups_inv.m
  mgrs_fwd.m
  mgrs_inv.m
  geoid_height.m
  geoid_load.m
  geocent_fwd.m
  geocent_inv.m
  loccart_fwd.m
  loccart_inv.m
  gedistance.m
  gereckon.m
  defaultellipsoid.m
  ecc2flat.m
  flat2ecc.m
  geographiclib_test.m
  geographiclib_signtest.m
  geoddoc.m
  projdoc.m
  gedoc.m)

set (PRIVATE_SOURCES
  private/A1m1f.m
  private/A2m1f.m
  private/A3coeff.m
  private/A3f.m
  private/AngDiff.m
  private/AngNormalize.m
  private/AngRound.m
  private/C1f.m
  private/C1pf.m
  private/C2f.m
  private/C3coeff.m
  private/C3f.m
  private/C4coeff.m
  private/C4f.m
  private/G4coeff.m
  private/GeoRotation.m
  private/LatFix.m
  private/SinCosSeries.m
  private/atan2dx.m
  private/cbrtx.m
  private/copysignx.m
  private/cvmgt.m
  private/eatanhe.m
  private/geoid_file.m
  private/geoid_load_file.m
  private/norm2.m
  private/remx.m
  private/sincosdx.m
  private/sincosde.m
  private/signbitx.m
  private/sumx.m
  private/tauf.m
  private/taupf.m)

configure_file (Contents.m.in ${DISTDIR}/inst/Contents.m @ONLY)
set (PACKAGE_DEPENDS ${PACKAGE_DEPENDS} ${DISTDIR}/inst/Contents.m)
foreach (f ${SOURCES} ${PRIVATE_SOURCES})
  configure_file (${f} ${DISTDIR}/inst/${f} COPYONLY)
  set (PACKAGE_DEPENDS ${PACKAGE_DEPENDS} ${DISTDIR}/inst/${f})
endforeach ()

install (FILES ${SOURCES} DESTINATION ${INSTDIR}/geographiclib)
install (FILES ${PRIVATE_SOURCES} DESTINATION ${INSTDIR}/geographiclib/private)
install (FILES ${DISTDIR}/inst/Contents.m DESTINATION ${INSTDIR}/geographiclib)

foreach (f ${METAFILES})
  install (FILES ${DISTDIR}/${f} DESTINATION ${INSTDIR}/geographiclib/packinfo)
endforeach ()

add_custom_command (OUTPUT ${TOOLBOX}
  WORKING_DIRECTORY ${DISTDIR}/inst
  COMMAND ${CMAKE_COMMAND} -E
  tar cf ${TOOLBOX} --format=zip -- Contents.m ${SOURCES} ${PRIVATE_SOURCES}
  COMMENT "Creating MATLAB toolbox dist/${TOOLBOX_NAME}"
  DEPENDS ${PACKAGE_DEPENDS})

add_custom_command (OUTPUT ${PACKAGE}
  WORKING_DIRECTORY ${DISTDIR}/..
  COMMAND ${CMAKE_COMMAND} -E
  tar cfz ${PACKAGE} --format=gnutar --
  geographiclib-octave-${PROJECT_VERSION}
  COMMENT "Creating Octave package dist/${PACKAGE_NAME}"
  DEPENDS ${PACKAGE_DEPENDS})

add_custom_target (package ALL DEPENDS ${TOOLBOX} ${PACKAGE})

find_program (OCTAVE octave)
if (OCTAVE)
  add_test (NAME testsuite
    COMMAND ${OCTAVE} --no-gui --no-window-system --eval geographiclib_test
    WORKING_DIRECTORY ${DISTDIR}/inst)
  add_test (NAME signtest
    COMMAND ${OCTAVE} --no-gui --no-window-system --eval geographiclib_signtest
    WORKING_DIRECTORY ${DISTDIR}/inst)
endif ()
