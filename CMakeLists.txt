#tkz-> use: cmake -DBOOST_ROOT="/home/hzpnotoad/Software/boost_1_59_0"


cmake_minimum_required(VERSION 3.0)
include(CPack)

project(psp)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#include(CMakeBuildTypes)

find_package(CPLEX REQUIRED)
if(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_LIBRARY})
endif()

if(CPLEX_CONCERT_FOUND)
    include_directories(${CPLEX_CONCERT_INCLUDE_DIR})
    link_directories(${CPLEX_CONCERT_LIBRARY})
endif()

if(CPLEX_ILOCPLEX_FOUND)
    include_directories(${CPLEX_ILOCPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_ILOCPLEX_LIBRARY})
endif()

if(CPLEX_CP_FOUND)
    include_directories(${CPLEX_CP_INCLUDE_DIR})
    link_directories(${CPLEX_CP_LIBRARY})
endif()


set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -Wall -Wunused -Wno-unknown-pragmas -fPIC -DIL_STD -fexceptions  -std=c++0x ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -lm -lpthread -ldl -mtune=native")

#find_package(cplex)
#find_package(cairo)
#find_package(GLUT)
#find_package(OpenGL)
#set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost 1.49.0 COMPONENTS program_options timer system chrono cplex)

#include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
#link_directories(${Boost_LIBRARY_DIRS})

#include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
#link_directories(${Boost_LIBRARY_DIRS})

#set(single_targets main mainIG mainInstChars mainRace mainTester mainPrepBatcher)
#set(single_targets main mainTester mainIg mainTimedIg mainRace mainTimedIgLbOrdered mainParallel)
#set(single_targets mainTester mainIg mainTimedIg mainRace mainTimedIgLbOrdered)
#set(single_targets main)
#foreach(target ${single_targets})
#  add_executable(${target} ${target}.cpp )
#  target_link_libraries(${target} ${Boost_LIBRARIES} rt)
#endforeach(target)

add_executable(main main.cpp Instance.cpp State.cpp Tabu.cpp Utilities.cpp) 
target_link_libraries(main ${Boost_LIBRARIES} rt cplex-library cplex-concert ilocplex cplex-cp )

#add_executable(intershop intershop.cpp jv/lap.cpp)
#target_link_libraries(intershop ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${Boost_LIBRARIES} rt)
#add_executable(gantt vis/gantt.cpp)
#target_link_libraries(gantt ${Boost_LIBRARIES} cairo rt)

## not working on header; left to remember
#set_source_files_properties(btst.hpp PROPERTIES COMPILE_FLAGS "-Wno-unused-but-set-variable -Wno-unused-variable")

set(CPACK_PACKAGE_EXECUTABLES "mainIG")

# Add new build types
message("* Adding build type \"profile\"")
set(CMAKE_CXX_FLAGS_PROFILE
    "${GCC_DEBUG_FLAGS} -pg"
    CACHE STRING "Flags used by the C++ compiler during profile builds."
    FORCE )
set(CMAKE_C_FLAGS_PROFILE
    "${GCC_DEBUG_FLAGS} -pg"
    CACHE STRING "Flags used by the C compiler during profile builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_PROFILE
    ""
    CACHE STRING "Flags used for linking binaries during profile builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE
    ""
    CACHE STRING "Flags used by the shared libraries linker during profile builds."
    FORCE )
mark_as_advanced(
    CMAKE_CXX_FLAGS_PROFILE
    CMAKE_C_FLAGS_PROFILE
    CMAKE_EXE_LINKER_FLAGS_PROFILE
    CMAKE_SHARED_LINKER_FLAGS_PROFILE )

# Add new build types
message("* Adding build type \"optprofile\"")
set(CMAKE_CXX_FLAGS_OPTPROFILE
    "-O -g -pg -DNDEBUG"
    CACHE STRING "Flags used by the C++ compiler during optimized profile builds."
    FORCE )
set(CMAKE_C_FLAGS_OPTPROFILE
    "-O3 -g -pg -DNDEBUG"
    CACHE STRING "Flags used by the C compiler during optimized profile builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_OPTPROFILE
    ""
    CACHE STRING "Flags used for linking binaries during optimized profile builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_OPTPROFILE
    ""
    CACHE STRING "Flags used by the shared libraries linker during optimized profile builds."
    FORCE )
mark_as_advanced(
    CMAKE_CXX_FLAGS_OPTPROFILE
    CMAKE_C_FLAGS_OPTPROFILE
    CMAKE_EXE_LINKER_FLAGS_OPTPROFILE
    CMAKE_SHARED_LINKER_FLAGS_OPTPROFILE )
