cmake_minimum_required(VERSION 3.14)
project(meshbConvert)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    libMeshb
    GIT_REPOSITORY https://github.com/LoicMarechal/libMeshb.git
    GIT_TAG        v7.80
)
FetchContent_MakeAvailable(libMeshb)

add_executable(stl2meshb src/stl2meshb.cpp)
target_link_libraries(stl2meshb PRIVATE Meshb.7)

add_executable(ugrid2meshb src/ugrid2meshb.cpp)
target_link_libraries(ugrid2meshb PRIVATE Meshb.7)
