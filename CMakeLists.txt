cmake_minimum_required(VERSION 3.10)

project(splitcode)

include(GNUInstallDirs)

set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(CMAKE_CXX_FLAGS_PROFILE "-g")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(${CMAKE_VERSION} VERSION_LESS 3.1)
    add_compile_options(-std=c++11)
    message("using C++11")
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
    if (COMPILER_SUPPORTS_CXX17)
        set(CMAKE_CXX_STANDARD 17)
        message("using C++17")
    else()
        message("using C++11")
    endif()
endif()

if(LINK MATCHES static)
    message("static build")
ELSE(LINK MATCHES shared)
    message("shared build")
ENDIF(LINK MATCHES static)

set(USE_HTSLIB "ON")
if (NO_HTSLIB)
    set(USE_HTSLIB OFF)
endif()

if(EMSCRIPTEN_BUILD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -s USE_ZLIB=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_ZLIB=1")
    set(USE_HTSLIB OFF)
endif()

add_subdirectory(src)

if(BUILD_FUNCTESTING)
    add_subdirectory(func_tests)
    message("Functional testing enabled.")
    add_custom_target(test
        COMMAND /bin/bash ./func_tests/runtests.sh
        DEPENDS ./src/splitcode
    )
endif(BUILD_FUNCTESTING)

include(ExternalProject)


if(NOT ZLIBNG_DISABLE)
    message("zlib-ng enabled.")
    set(ZLIBNG "ON")
    ExternalProject_Add(zlib-ng
        PREFIX ${PROJECT_SOURCE_DIR}/ext/zlib-ng
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/zlib-ng
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND
            mkdir -p zlib-ng && cd zlib-ng &&
            cmake .. -DZLIB_COMPAT=ON -DZLIB_ENABLE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${PREFIX}
        BUILD_COMMAND cd zlib-ng && make
        INSTALL_COMMAND ""
    )
endif(NOT ZLIBNG_DISABLE)


if(USE_HTSLIB)
    message("using htslib")
    ExternalProject_Add(htslib
        PREFIX ${PROJECT_SOURCE_DIR}/ext/htslib
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/htslib
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make lib-static
        INSTALL_COMMAND ""
    )
else()
    set(NO_HTSLIB "ON")
    add_compile_definitions("NO_HTSLIB=ON")
    message("not using htslib")
endif()
