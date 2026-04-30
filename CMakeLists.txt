cmake_minimum_required(VERSION 3.10)
project(meshbConvert)

set(LIBMESHB_DIR "$ENV{HOME}/progs/libMeshb/sources")

add_executable(stl2meshb src/stl2meshb.cpp ${LIBMESHB_DIR}/libmeshb7.c)
target_include_directories(stl2meshb PRIVATE ${LIBMESHB_DIR})

add_executable(ugrid2meshb src/ugrid2meshb.cpp ${LIBMESHB_DIR}/libmeshb7.c)
target_include_directories(ugrid2meshb PRIVATE ${LIBMESHB_DIR})
