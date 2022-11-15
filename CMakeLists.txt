cmake_minimum_required(VERSION 3.18)

project("osc-daemon")

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES  
    ${CMAKE_CURRENT_SOURCE_DIR}/osc-daemon/main.cpp
    )
set(HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/osc-daemon/osctypes.hpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
