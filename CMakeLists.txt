cmake_minimum_required(VERSION 2.8)

INCLUDE_DIRECTORIES(/usr/include/)

   SET(MAP_LIST ./src/main.cpp
        ./src/Types.h
        ./src/SplitString.h
        ./src/SplitString.cc
        ./src/map.h
        ./src/map.cpp
        ./src/mole.h
        ./src/mole.cpp
        ./src/gene.cpp
        ./src/gene.h
        ./src/constant.h
        ./src/constant.cpp
        )
ADD_EXECUTABLE(map ${MAP_LIST})

target_link_libraries(map boost_filesystem boost_system log4cxx)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wno-unused-variable -pthread -lboost_system -lboost_filesystem -llog4cxx")
INSTALL(TARGETS map RUNTIME DESTINATION bin)
