# Minimum version of CMake
cmake_minimum_required(VERSION 2.8)

# Set C++ Version
#SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
#SET(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS ON) #...w/o or w/ compiler extensions like gnu++11

# Project name
project(css-or)

# OpenCV
find_package(OpenCV REQUIRED)

# Build css_or_test
add_executable(css_or_test
  src/css_or_test.cpp
  src/engine/EdgeProcessor.cpp
  src/engine/EdgeAnalyzer.cpp
  src/util/Edge.cpp
  src/util/InterestPoint.cpp)
  
# Build css_or_test2
add_executable(css_or_test2
  src/css_or_test2.cpp
  src/engine/EdgeProcessor.cpp
  src/engine/EdgeAnalyzer.cpp
  src/util/Edge.cpp
  src/util/InterestPoint.cpp)

# Build css_or_test3
add_executable(css_or_test3
  src/css_or_test3.cpp
  src/engine/EdgeProcessor.cpp
  src/engine/CurvatureScaleSpace.cpp
  src/util/Edge.cpp
  src/util/InterestPoint.cpp
  src/util/Visualization.cpp
  src/util/ShapeFeatures.cpp)
  
 # Build tracing_test
add_executable(tracing_test
  src/tracing_test.cpp
  src/engine/EdgeProcessor.cpp
  src/util/Edge.cpp
  src/util/InterestPoint.cpp
  src/util/Visualization.cpp)

# Build css_appr_test
add_executable(css_appr_test
  src/css_appr_test.cpp
  src/engine/EdgeProcessor.cpp
  src/engine/CurvatureScaleSpace.cpp
  src/util/Edge.cpp
  src/util/InterestPoint.cpp
  src/util/Visualization.cpp
  src/util/ShapeFeatures.cpp
  src/engine/ApprCSS.cpp)

# Linking
target_link_libraries(css_or_test ${OpenCV_LIBS})
target_link_libraries(css_or_test2 ${OpenCV_LIBS})
target_link_libraries(css_or_test3 ${OpenCV_LIBS})
target_link_libraries(tracing_test ${OpenCV_LIBS})
target_link_libraries(css_appr_test ${OpenCV_LIBS})
