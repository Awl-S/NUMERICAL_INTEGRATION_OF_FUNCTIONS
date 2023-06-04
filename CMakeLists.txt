cmake_minimum_required(VERSION 3.25)
project(NUMERICAL_INTEGRATION_OF_FUNCTIONS)

set(CMAKE_CXX_STANDARD 17)

add_executable(NUMERICAL_INTEGRATION_OF_FUNCTIONS main.cpp src/NumericalIntegrator.h)

find_package(GSL REQUIRED)
target_link_libraries(NUMERICAL_INTEGRATION_OF_FUNCTIONS PRIVATE GSL::gsl GSL::gslcblas)