cmake_minimum_required(VERSION 3.24)
project(pdp_project)

set(CMAKE_CXX_STANDARD 23)

add_executable(pdp_project main.cpp src/objects/Vector.h src/objects/Atom.h)
