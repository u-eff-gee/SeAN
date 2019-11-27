find_package(ROOT REQUIRED)

add_library(lib_VoigtProfile VoigtProfile.cc)
target_link_libraries(lib_VoigtProfile ${ROOT_LIBRARIES})