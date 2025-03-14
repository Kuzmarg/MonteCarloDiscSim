cmake_minimum_required(VERSION 3.5)
project(MonteCarloParticles LANGUAGES C)

add_executable(simulation src/main.c
    src/cell.c src/utils.c src/particles.c src/simulation.c src/movement.c)
set_property(TARGET simulation PROPERTY C_STANDARD 11)
target_link_libraries(simulation m)
target_include_directories(simulation PRIVATE include)
