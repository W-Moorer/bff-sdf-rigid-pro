cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{VCPKG_ROOT})
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "vcpkg toolchain file")
endif()

project(bff_sdf_contact_atlas LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BFF_SDF_BUILD_TESTS "Build unit tests" ON)
option(BFF_SDF_BUILD_EXPERIMENTS "Build benchmark executables" ON)

find_package(Eigen3 CONFIG REQUIRED)

add_library(bff_sdf_contact
    src/core/Mesh.cpp
    src/io/MeshIO.cpp
    src/sdf/AnalyticSDF.cpp
    src/atlas/PatchBuilder.cpp
    src/atlas/UVAcceleration.cpp
    src/atlas/BFFAdapter.cpp
    src/atlas/ContactAtlas.cpp
    src/atlas/ChartEvaluator.cpp
)

target_include_directories(bff_sdf_contact PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(bff_sdf_contact PUBLIC Eigen3::Eigen)

if(MSVC)
    target_compile_options(bff_sdf_contact PRIVATE /W4 /permissive-)
    target_compile_definitions(bff_sdf_contact PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
    target_compile_options(bff_sdf_contact PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)

function(add_contact_test test_name source_file)
    add_executable(${test_name} ${source_file})
    target_link_libraries(${test_name} PRIVATE bff_sdf_contact)
    target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

if(BFF_SDF_BUILD_TESTS)
    add_contact_test(test_mesh tests/test_mesh.cpp)
    add_contact_test(test_pose tests/test_pose.cpp)
    add_contact_test(test_sdf tests/test_sdf.cpp)
    add_contact_test(test_bff_adapter tests/test_bff_adapter.cpp)
    add_contact_test(test_chart_evaluator tests/test_chart_evaluator.cpp)
endif()
