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
target_include_directories(nca PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CNTR_INCLUDE_DIR} ${CNTR_HDF5_INCLUDE_DIR} )

# Propagate usage requirements to dependents
target_link_libraries(nca 
  ppsc
  Eigen3::Eigen
  ${CNTR_LIBRARY}
  ${CNTR_HDF5_LIBRARY}
  OpenMP::OpenMP_CXX
)
# ~~ Installation targets ~~
install(TARGETS nca
        LIBRARY DESTINATION lib)
