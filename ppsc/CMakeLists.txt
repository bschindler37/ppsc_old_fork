
add_subdirectory(nca)
add_subdirectory(oca)

include_directories(
  SYSTEM
  ../lib/cereal-1.2.1/include/
  ../
  )

set(OBJECT_FILES
   ppsc.cpp
   range.cpp
   )

# ~~ The actual target library ~~

# ~~ PPSC LIB ~~

add_library(ppsc SHARED ${OBJECT_FILES})
set_property(TARGET ppsc APPEND_STRING PROPERTY COMPILE_FLAGS " -w")
target_include_directories(ppsc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CNTR_INCLUDE_DIR} ${CNTR_HDF5_INCLUDE_DIR} )

# Propagate usage requirements to dependents
target_link_libraries(ppsc PUBLIC 
  Eigen3::Eigen
  ${CNTR_LIBRARY}
  ${CNTR_HDF5_LIBRARY}
  OpenMP::OpenMP_CXX
)
# ~~ Installation targets ~~
install(TARGETS ppsc
        LIBRARY DESTINATION lib)
