cmake_minimum_required(VERSION 3.5)
project(MonteCarloParticles LANGUAGES C)

add_executable(cellgen_demo src/cellgen_test.c
    src/cell.c src/utils.c src/particles.c)
set_property(TARGET cellgen_demo PROPERTY C_STANDARD 11)
target_link_libraries(cellgen_demo m)
target_include_directories(cellgen_demo PRIVATE include)
