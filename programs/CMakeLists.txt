include_directories(
  SYSTEM
  ../lib/cereal-1.2.1/include/
  )

# ~~ Programs

if (hdf5 AND mpi)

  # ~~ hdf5 programs

  set( program_sources
    #single_band_hubbard_holstein_MYFILEOLD.cpp
    test_boson_hilbert.cpp
    test_dense_hilbert.cpp
    test_tJ_hilbert.cpp
    test_tJ_twoband_hilbert.cpp
    single_band_bethe.cpp
    single_band_bose_hubbard_bethe.cpp
    single_band_hubbard_bethe.cpp
    single_band_hubbard_bethe_test_selfe.cpp
    single_band_hubbard_bethe_spin.cpp
    single_band_hubbard_bethe_sc.cpp
    single_band_hubbard_holstein_bethe.cpp
    #single_band_hubbard_holstein_bethe_lang_firsov.cpp
    single_band_hubbard_bosonbath_bethe.cpp
    two_band_hubbard_bethe.cpp
    two_band_hubbard_bethe_field.cpp
    #two_band_tJ_bethe_field_afm.cpp
    #two_band_tJ_bethe_field_afm_positive.cpp
    #two_band_tJ_bethe_field_afm_twostates.cpp
    two_band_hubbard_bethe_field_current_sk.cpp
    two_band_hubbard_bethe_field_scsusc_dense.cpp
    #three_band_hubbard_bethe.cpp
    #three_band_hubbard_bethe_field.cpp
    #LVO_hetero.cpp
    # -- NOT WORKING FIXME!
    # n_channel_siam.cpp
    # one_channel_siam_drive.cpp
    one_spinless_fermion.cpp
    one_spinful_fermion.cpp
    #xas_two_band_MYFILEOLD.cpp
    xas_two_band_atomic.cpp
    xas_two_band_corebath-expansion.cpp
    #xas_two_band_corebath-expansion_MYFILEOLD.cpp
    )

  foreach(source_file ${program_sources})

    STRING(REGEX REPLACE ".cpp" ".ex" executable_file "${source_file}")  
    add_executable(${executable_file} ${source_file})
    target_link_libraries(${executable_file} cntr nca oca hdf5 ${COMMON_LIBS})
    install(TARGETS ${executable_file} DESTINATION bin)

  endforeach(source_file)

endif (hdf5 AND mpi)
