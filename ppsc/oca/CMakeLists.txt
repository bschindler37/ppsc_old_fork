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

target_include_directories(oca PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CNTR_INCLUDE_DIR} ${CNTR_HDF5_INCLUDE_DIR} )

# Propagate usage requirements to dependents
target_link_libraries(oca 
  ppsc
  Eigen3::Eigen
  ${CNTR_LIBRARY}
  ${CNTR_HDF5_LIBRARY}
  OpenMP::OpenMP_CXX
)
# ~~ Installation targets ~~
install(TARGETS oca
        LIBRARY DESTINATION lib)
