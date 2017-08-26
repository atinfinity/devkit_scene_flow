cmake_minimum_required(VERSION 2.8)
project(evaluate_scene_flow)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-${CMAKE_CXX_FLAGS} -std=c++11")
endif()

#set(OpenCV_DIR C:/dev/opencv-3.2.0/build/install)
find_package(OpenCV REQUIRED)

if(OpenCV_FOUND)
    include_directories(${OpenCV_INCLUDE_DIRS})
    add_executable(evaluate_scene_flow evaluate_scene_flow.cpp)    
    target_link_libraries(evaluate_scene_flow ${OpenCV_LIBS})
endif(OpenCV_FOUND)
