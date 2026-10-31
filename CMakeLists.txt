# CMakeList.txt : CMake project for reservoir_project, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)
project(reservoir_simulator CXX)
set(CMAKE_CXX_STANDARD 20)

# Enable exception handling
if (MSVC)
  add_compile_options(/EHsc)
  # Static link MSVC runtime to avoid DLL dependencies
  string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
  string(REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif()

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

# Add source to this project's executable.
#add_executable (reservoir_simulator "reservoir_project.cpp" "reservoir_project.h")

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET reservoir_simulator PROPERTY CXX_STANDARD 20)
#endif()

# TODO: Add tests and install targets if needed.


# Add packages

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

include(FetchContent)
FetchContent_Declare(
	Eigen3
	GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
	GIT_TAG 5.0.1
	GIT_SHALLOW TRUE
)
FetchContent_Declare(
	pybind11
	GIT_REPOSITORY https://github.com/pybind/pybind11.git
	GIT_TAG v3.0.4
	GIT_SHALLOW TRUE
)
#FetchContent_Declare(
#	spdlog
#	GIT_REPOSITORY https://github.com/gabime/spdlog.git
#	GIT_TAG v1.13.0
#	GIT_SHALLOW TRUE
#)
FetchContent_MakeAvailable(Eigen3 pybind11)# spdlog)

add_library(simcore INTERFACE)
target_include_directories(simcore INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(simcore INTERFACE Eigen3::Eigen)# spdlog::spdlog)

pybind11_add_module(pyressim bindings/py_bindings.cpp)
target_link_libraries(pyressim PRIVATE simcore Eigen3::Eigen)

add_executable(test_grid test/test_grid.cpp)
target_link_libraries(test_grid PRIVATE simcore)
enable_testing()
add_test(NAME test_grid COMMAND test_grid)

add_executable(test_relperm test/test_relperm.cpp)
target_link_libraries(test_relperm PRIVATE simcore)
add_test(NAME test_relperm COMMAND test_relperm)

add_executable(test_fluid test/test_fluid.cpp)
target_link_libraries(test_fluid PRIVATE simcore)
add_test(NAME test_fluid COMMAND test_fluid)

add_executable(test_impes test/test_impes.cpp)
target_link_libraries(test_impes PRIVATE simcore)
add_test(NAME test_impes COMMAND test_impes)

add_executable(test_newton test/test_newton.cpp)
target_link_libraries(test_newton PRIVATE simcore)
add_test(NAME test_newton COMMAND test_newton)

add_executable(test_model test/test_model.cpp)
target_link_libraries(test_model PRIVATE simcore)
add_test(NAME test_model COMMAND test_model)

add_executable(test_simulator test/test_simulator.cpp)
target_link_libraries(test_simulator PRIVATE simcore)
add_test(NAME test_simulator COMMAND test_simulator)

add_executable(diag test/diag.cpp)
target_link_libraries(diag PRIVATE simcore)
add_test(NAME diag COMMAND diag)