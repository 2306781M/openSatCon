cmake_minimum_required(VERSION 3.18)

project("osc-daemon")

set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_COMPILER /usr/bin/x86_64-w64-mingw32-g++)

find_package(Threads REQUIRED)

add_subdirectory(external)

set(SOURCE_FILES  
    ${CMAKE_CURRENT_SOURCE_DIR}/osc-daemon/main.cpp
    )
set(HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/osc-daemon/osctypes.hpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_compile_definitions(${PROJECT_NAME} PUBLIC _CONSOLE)
target_link_libraries(${PROJECT_NAME} TBB::tbb Threads::Threads)
target_compile_options(${PROJECT_NAME} PRIVATE ${TBB_CXX_STD_FLAG})
target_include_directories(${PROJECT_NAME} PUBLIC external/oneTBB/include)
target_include_directories(${PROJECT_NAME} PUBLIC external/gcem/include)