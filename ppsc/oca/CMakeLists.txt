include_directories(
  SYSTEM
  ../../lib/cereal-1.2.1/include/
  ../../cntr/
  ../
  )
  
set(OBJECT_FILES
   sigma.cpp
   sigma_disp.cpp
   gf.cpp
   gf_disp.cpp
   )

#   ../ppsc/cntr_types.cpp

# ~~ The actual target library ~~

# ~~ NEW OCA MATRIX OPERATOR LIB ~~

add_library(oca SHARED ${OBJECT_FILES})
set_property(TARGET oca APPEND_STRING PROPERTY COMPILE_FLAGS " -w")
target_link_libraries(oca cntr ppsc ${COMMON_LIBS})

if (hdf5)
#  target_link_libraries(oca cntr_hdf5)
   target_link_libraries(oca)
endif (hdf5)

# ~~ Installation targets ~~
install(TARGETS oca
        LIBRARY DESTINATION lib)
