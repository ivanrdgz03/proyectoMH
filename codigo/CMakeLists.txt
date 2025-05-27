cmake_minimum_required(VERSION 3.0)
project(cec17)

add_executable(woa woa.cpp)
add_executable(woa-hybrid woa-hibridado.cpp)

add_library(cec17_test_func SHARED cec17_test_func.c cec17.c)

target_link_libraries(woa cec17_test_func)
target_link_libraries(woa-hybrid cec17_test_func)

target_include_directories(woa-hybrid PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
