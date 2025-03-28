cmake_minimum_required(VERSION 3.16)
project(Runner)

set(CMAKE_CXX_STANDARD 17)

find_package(Geant4 REQUIRED)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")

# set(LibPath ${CMAKE_CURRENT_SOURCE_DIR}/..)
# add_subdirectory(${LibPath} ${CMAKE_BINARY_DIR}/Geant4)
include(${Geant4_USE_FILE})

add_executable(Runner main.cpp)
target_link_libraries(Runner PUBLIC ${Geant4_LIBRARIES})
target_include_directories(Runner PUBLIC ${Geant4_INCLUDE_DIR})
