# Install script for directory: /Users/admin/Documents/phd/codes/ppsc/programs

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/test_boson_hilbert.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_boson_hilbert.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_boson_hilbert.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_boson_hilbert.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_boson_hilbert.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/test_boson_hilbert.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/test_dense_hilbert.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_dense_hilbert.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_dense_hilbert.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_dense_hilbert.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_dense_hilbert.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/test_dense_hilbert.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/test_tJ_hilbert.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_hilbert.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_hilbert.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_hilbert.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_hilbert.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/test_tJ_hilbert.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/test_tJ_twoband_hilbert.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_twoband_hilbert.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_twoband_hilbert.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_twoband_hilbert.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_tJ_twoband_hilbert.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/test_tJ_twoband_hilbert.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_bethe.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bethe.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bethe.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bethe.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bethe.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_bethe.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_bose_hubbard_bethe.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bose_hubbard_bethe.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bose_hubbard_bethe.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bose_hubbard_bethe.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_bose_hubbard_bethe.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_bose_hubbard_bethe.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_hubbard_bethe.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_hubbard_bethe.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_hubbard_bethe_test_selfe.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_test_selfe.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_test_selfe.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_test_selfe.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_test_selfe.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_hubbard_bethe_test_selfe.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_hubbard_bethe_spin.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_spin.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_spin.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_spin.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_spin.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_hubbard_bethe_spin.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_hubbard_bethe_sc.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_sc.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_sc.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_sc.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_sc.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_hubbard_bethe_sc.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_hubbard_holstein_bethe.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_holstein_bethe.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_holstein_bethe.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_holstein_bethe.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_holstein_bethe.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_hubbard_holstein_bethe.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_hubbard_bosonbath_bethe.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bosonbath_bethe.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bosonbath_bethe.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bosonbath_bethe.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bosonbath_bethe.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_hubbard_bosonbath_bethe.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/two_band_hubbard_bethe.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/two_band_hubbard_bethe.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/two_band_hubbard_bethe_field.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/two_band_hubbard_bethe_field.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/single_band_hubbard_bethe_driven.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_driven.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_driven.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_driven.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/single_band_hubbard_bethe_driven.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/single_band_hubbard_bethe_driven.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/two_band_hubbard_bethe_field_current_sk.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_current_sk.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_current_sk.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_current_sk.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_current_sk.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/two_band_hubbard_bethe_field_current_sk.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/two_band_hubbard_bethe_field_scsusc_dense.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_scsusc_dense.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_scsusc_dense.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_scsusc_dense.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/two_band_hubbard_bethe_field_scsusc_dense.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/two_band_hubbard_bethe_field_scsusc_dense.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/one_spinless_fermion.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/one_spinless_fermion.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/one_spinless_fermion_bosonic-bath.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bath.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bath.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bath.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bath.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/one_spinless_fermion_bosonic-bath.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/one_spinless_fermion_bosonic-bathTEST.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bathTEST.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bathTEST.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bathTEST.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinless_fermion_bosonic-bathTEST.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/one_spinless_fermion_bosonic-bathTEST.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/one_spinful_fermion.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinful_fermion.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinful_fermion.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinful_fermion.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/one_spinful_fermion.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/one_spinful_fermion.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/xas_two_band_atomic.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_atomic.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_atomic.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_atomic.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_atomic.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/xas_two_band_atomic.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/xas_two_band_corebath-expansion.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/xas_two_band_corebath-expansion.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/phd/codes/ppsc/build/programs/xas_two_band_corebath-expansion2.ex")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion2.ex" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion2.ex")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/$"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/nca"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc/oca"
      -delete_rpath "/usr/local/lib"
      -delete_rpath "/Users/admin/Documents/phd/codes/ppsc/build/ppsc"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion2.ex")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xas_two_band_corebath-expansion2.ex")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/admin/Documents/phd/codes/ppsc/build/programs/CMakeFiles/xas_two_band_corebath-expansion2.ex.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/admin/Documents/phd/codes/ppsc/build/programs/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
