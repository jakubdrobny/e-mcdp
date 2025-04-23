cmake_minimum_required(VERSION 3.31)
project(emcdp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")

add_executable(emcdp ${SOURCES})

target_compile_options(emcdp PRIVATE
  -Wall
  -O3
)

target_link_libraries(emcdp PRIVATE
  OpenMP::OpenMP_CXX
  GTest::GTest
)

install(TARGETS emcdp
  RUNTIME DESTINATION bin
  COMPONENT runtime
)
