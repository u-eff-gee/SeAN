find_package(ROOT REQUIRED)

add_library(lib_ExcitedState ExcitedState.cc)
target_link_libraries(lib_ExcitedState ${ROOT_LIBRARIES})

add_library(lib_Grid Grid.cc)
target_link_libraries(lib_Grid ${ROOT_LIBRARIES})

add_library(lib_Integrator Integrator.cc)
target_link_libraries(lib_Integrator ${ROOT_LIBRARIES})

add_library(lib_Nucleus Nucleus.cc)
target_link_libraries(lib_Nucleus ${ROOT_LIBRARIES})
target_link_libraries(lib_Nucleus lib_Grid)

add_library(lib_Target Target.cc)
target_link_libraries(lib_Target ${ROOT_LIBRARIES})

add_library(lib_VoigtProfile VoigtProfile.cc)
target_link_libraries(lib_VoigtProfile ${ROOT_LIBRARIES})