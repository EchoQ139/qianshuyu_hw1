cmake_minimum_required(VERSION 3.10)
project(algebra)

# 设置C标准
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 添加可执行文件
add_executable(algebra
    src/algebra.c
    src/main.c
)

# 包含头文件目录
target_include_directories(algebra PRIVATE inc)

# 链接数学库
target_link_libraries(algebra m)