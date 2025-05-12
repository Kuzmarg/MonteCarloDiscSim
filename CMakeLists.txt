cmake_minimum_required(VERSION 3.18)
project(MonteCarloParticles LANGUAGES CUDA)

find_package(libconfig REQUIRED)

add_executable(simulation src/main.cu
    src/cell.cu
    src/utils.cu
    src/particles.cu
    src/simulation.cu
    src/movement.cu)
set_property(TARGET simulation PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET simulation PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
target_link_libraries(simulation m ${LIBCONFIG_LIBRARIES})
target_include_directories(simulation PRIVATE include)
