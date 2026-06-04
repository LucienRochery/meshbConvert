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

add_executable(meshConvert src/meshConvert.cpp)
target_link_libraries(meshConvert PRIVATE Meshb.7)
