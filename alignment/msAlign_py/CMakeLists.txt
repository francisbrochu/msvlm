cmake_minimum_required(VERSION 3.0)
project(msAlignForPy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}) #for shared library
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}) #for the executable
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}) #for static library

find_package (Boost)
if (Boost_FOUND)
    message("Boost found!")
else ()
    message(FATAL_ERROR "Boost not found")
endif ()

find_package(PythonLibs)
if (PYTHONLIBS_FOUND)
    message("PythonLibs found!")
else ()
    message(FATAL_ERROR "PythonLibs not found")
endif ()

include_directories(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

set(SOURCE_FILES msAlignForPy.cpp ActiveSequence.cpp Heap.cpp)
add_library(msAlignForPy MODULE ${SOURCE_FILES})

target_link_libraries(msAlignForPy ${BOOST_LIBRARIES} boost_python3 ${PYTHON_LIBRARIES})
