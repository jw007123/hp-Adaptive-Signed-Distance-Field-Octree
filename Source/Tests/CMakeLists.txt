add_executable(Tests "Main.cpp")

# .obj resources
file(COPY "../../Resources" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(Tests LINK_PUBLIC Lib_HP)