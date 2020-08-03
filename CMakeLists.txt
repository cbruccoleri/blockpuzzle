cmake_minimum_required(VERSION 3.11...3.16)

project(blockpuzzle VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(${PROJECT_NAME} blockpuzzle.cpp olcPixelGameEngine.h)

# Handle different requirements on Linux and Windows
# On Windows with the Visual Studio generator, the code has pragma directives to add these libraries,
# so nothing explicit is needed in CMake.
if (NOT MSVC)
    # On Linux, the following libraries should be already installed in most distributions.
    find_package(X11 REQUIRED)
    find_package(OpenGL REQUIRED)
    find_package(PNG REQUIRED)
    find_package(Threads REQUIRED)
    set(CMAKE_THREAD_PREFER_PTHREAD ON)

    # The link targets are made PRIVATE because they should not necessarily affect dependent projects.
    # If you add dependencies to this CMake build, you may need to change this.
    target_link_libraries(${PROJECT_NAME} PRIVATE ${X11_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL OpenGL::OpenGL OpenGL::GLU OpenGL::GLX)
    target_link_libraries(${PROJECT_NAME} PRIVATE PNG::PNG)
    target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)

    # GNU g++ requires an explicit linking flag for std::filesystem. 
    # This is a workaround because neither GNU or Kitware have achieved consensus on
    # how to handle this standard library.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(STATUS "GNU Compiler detected: linking with std::filesystem explicitly.")        
        target_link_libraries(${PROJECT_NAME} PRIVATE stdc++fs)
    endif()
    
endif()