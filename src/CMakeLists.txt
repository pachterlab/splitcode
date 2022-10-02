file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

add_library(splitcode_core ${sources} ${headers})
target_include_directories(splitcode_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(splitcode main.cpp)

find_package( Threads REQUIRED )
target_link_libraries(splitcode splitcode_core pthread)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)

    if (UNIX AND NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libstdc++")
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

    SET_TARGET_PROPERTIES(splitcode splitcode_core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINK MATCHES static)

find_package( ZLIB REQUIRED )

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries(splitcode splitcode_core ${ZLIB_LIBRARIES})
else()
    message(FATAL_ERROR "zlib not found. Required for reading FASTQ files" )
endif( ZLIB_FOUND )

if(LINK MATCHES static)
    if (UNIX AND NOT APPLE)
        target_link_libraries(splitcode librt.a)
    endif()
else()
    if (UNIX AND NOT APPLE)
        target_link_libraries(splitcode rt)
    endif()
endif(LINK MATCHES static)


install(TARGETS splitcode DESTINATION "${CMAKE_INSTALL_BINDIR}")