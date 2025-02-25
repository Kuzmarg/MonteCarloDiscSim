cmake_minimum_required(VERSION 3.5)
project(MonteCarloParticles LANGUAGES C)

include_directories(include)

add_library(cellgen src/cell.c)
target_link_libraries(cellgen m)
set_property(TARGET cellgen PROPERTY C_STANDARD 11)

add_executable(cellgen_demo src/cellgen_test.c)
target_link_libraries(cellgen_demo cellgen)
set_property(TARGET cellgen_demo PROPERTY C_STANDARD 11)
