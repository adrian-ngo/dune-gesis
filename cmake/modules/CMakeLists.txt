set(modules DuneGesisMacros.cmake FindFFTW3.cmake AddFFTW3Flags.cmake AddHDF5Flags.cmake)

install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
