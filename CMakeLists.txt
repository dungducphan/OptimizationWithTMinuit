cmake_minimum_required(VERSION 3.13)
project(UnconstrainedOptimization)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

add_executable(UnconstrainedOptimization main.cpp)
target_link_libraries(UnconstrainedOptimization ${ROOT_LIBRARIES} Minuit)