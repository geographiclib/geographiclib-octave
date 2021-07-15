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
  private/sumx.m
  private/swap.m
  private/tauf.m
  private/taupf.m)

configure_file (Contents.m.in Contents.m @ONLY)
foreach (f ${SOURCES} ${PRIVATE_SOURCES})
  configure_file (${f} ${f} COPYONLY)
endforeach ()

install (FILES ${SOURCES} DESTINATION ${INSTALL_PATH}/geographiclib)
install (FILES ${PRIVATE_SOURCES}
  DESTINATION ${INSTALL_PATH}/geographiclib/private)

set (PACKAGE
  ${PROJECT_BINARY_DIR}/dist/geographiclib_toolbox_${PROJECT_VERSION}.zip)

add_custom_command (OUTPUT ${PACKAGE}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/dist
  COMMAND ${CMAKE_COMMAND} -E
  tar cf ${PACKAGE} --format=zip -- Contents.m ${SOURCES} ${PRIVATE_SOURCES}
  COMMENT "Creating package")

add_custom_target (package ALL DEPENDS ${PACKAGE})

find_program (OCTAVE octave)
if (OCTAVE)
  add_test (NAME testsuite
    COMMAND ${OCTAVE} --no-gui --no-window-system --eval geographiclib_test)
endif ()
  
