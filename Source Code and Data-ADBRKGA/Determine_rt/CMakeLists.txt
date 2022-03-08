cmake_minimum_required(VERSION 3.9)
project(WFGA)

include_directories(${PROJECT_SOURCE_DIR}/include)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "-fopenmp")
add_subdirectory(${PROJECT_SOURCE_DIR}/lib)
add_executable(WFGA main.cpp include/config.h include/GenerateAChrom.h CGA.cpp CGA.h HGA.cpp HGA.h LWSGA.cpp LWSGA.h HEFT.cpp HEFT.h ADBRKGA.cpp ADBRKGA.h MOELS.cpp MOELS.h IHEFT3.cpp IHEFT3.h DHEFT.cpp DHEFT.h HPSO.cpp HPSO.h)
target_link_libraries(WFGA tools GenOperator pugixml classAndVarDefine config)