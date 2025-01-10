cmake_minimum_required(VERSION 3.10)

project(MyJacobiSolver LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall")

set(APPS jacobi_noft)

include_directories(${CMAKE_SOURCE_DIR}/include)
set(SOURCES src/main.c src/jacobi_noft.c)

find_package(MPI REQUIRED)
find_package(PNG REQUIRED)

add_executable(${APPS} ${SOURCES})

target_link_libraries(${APPS} PRIVATE MPI::MPI_C PNG::PNG m pthread)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")

add_custom_target(clean_all
  COMMAND ${CMAKE_COMMAND} -E remove -f ${APPS} ${APPS}.o main.o jacobi_noft.o
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Cleaning up build artifacts"
)
