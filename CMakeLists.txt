cmake_minimum_required(VERSION 2.8.12)

project(splitcode)

include(GNUInstallDirs)

set(CMAKE_CXX_FLAGS_PROFILE "-g")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(${CMAKE_VERSION} VERSION_LESS 3.1)
    add_compile_options(-std=c++11)
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

if(LINK MATCHES static)
    message("static build")
ELSE(LINK MATCHES shared)
    message("shared build")
ENDIF(LINK MATCHES static)


add_subdirectory(src)

if (BUILD_FUNCTESTING)
    add_subdirectory(func_tests)
    message("Functional testing enabled.")
    add_custom_target(test
    COMMAND /bin/bash ./func_tests/runtests.sh
    DEPENDS ./src/splitcode
    )
endif(BUILD_FUNCTESTING)

