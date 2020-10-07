FILE(GLOB_RECURSE headers *.hh *.h)
FILE(GLOB_RECURSE sources *.cc *.cpp)

set(CMAKE_CXX_FLAGS "-O3 -ggdb ") # no werror here!

add_library (Goldfarb_Optimizer SHARED ${headers} ${sources})
#add_library (Goldfarb_Optimizer STATIC ${headers} ${sources})

#install(TARGETS Goldfarb_Optimizer DESTINATION "${INSTALL_LIB_DIR}")
#install(FILES ${headers} DESTINATION
    #"${INSTALL_INCLUDE_DIR}/Optimizer/Goldfarb_Optimizer")
