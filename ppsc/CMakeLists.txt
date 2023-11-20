
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
target_link_libraries(ppsc cntr)

if (hdf5)
   target_link_libraries(ppsc cntr_hdf5)
   target_link_libraries(ppsc cntr)
endif (hdf5)

# ~~ Installation targets ~~
install(TARGETS ppsc
        LIBRARY DESTINATION lib)
