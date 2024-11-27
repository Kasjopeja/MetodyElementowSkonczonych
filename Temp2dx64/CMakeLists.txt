cmake_minimum_required(VERSION 3.29)
project(Temp2dx64)

set(CMAKE_CXX_STANDARD 20)

add_executable(Temp2dx64 main.cpp
        MyType.h
        IniEL4.cpp
        InpData.cpp
        GlobData.h
        IniEL4.h
        InpData.h
        GenGrid2d.cpp
        GenGrid2d.h
)
