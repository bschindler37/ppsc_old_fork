include_directories(
  SYSTEM
  ../../lib/cereal-1.2.1/include/
  ../../cntr/
  ../
  )
  
set(OBJECT_FILES
   sigma.cpp
   gf.cpp
   )

# ~~ The actual target library ~~

# ~~ NCA MATRIX OPERATOR LIB ~~

add_library(nca SHARED ${OBJECT_FILES})
set_property(TARGET nca APPEND_STRING PROPERTY COMPILE_FLAGS " -w")
target_link_libraries(nca cntr ppsc ${COMMON_LIBS})

if (hdf5)
#  target_link_libraries(nca cntr_hdf5)
  target_link_libraries(nca )
endif (hdf5)

# ~~ Installation targets ~~
install(TARGETS nca
        LIBRARY DESTINATION lib)
