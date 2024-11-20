cmake_minimum_required(VERSION 3.29)
project(Temp2dx64)

set(CMAKE_CXX_STANDARD 20)

add_executable(Temp2dx64 main.cpp
        GlobData.cpp
        MyType.cpp
        IniEL4.cpp
)
