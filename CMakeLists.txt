cmake_minimum_required(VERSION 3.18)
project(MonteCarloParticles LANGUAGES CUDA)

add_executable(simulation src/main.cu
    src/cell.cu src/utils.cu src/particles.cu src/simulation.cu src/movement.cu)
target_link_libraries(simulation m)
target_include_directories(simulation PRIVATE include)
