cmake_minimum_required(VERSION 3.5.1)

project(POPPYProject LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

option(BUILD_STATIC "Set to ON to include static versions of the library" OFF)
option(WITH_TCM "Activate tcmalloc by setting WITH_TCM to ON" ON)
option(METIS "Set to ON to include METIS" OFF)

# ===============================
# 查找 OpenFHE
# ===============================
find_package(OpenFHE CONFIG REQUIRED)
if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif()

string(REPLACE "-Werror " "" CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS})

include_directories(${OPENMP_INCLUDES})
include_directories(${OpenFHE_INCLUDE})
include_directories(${OpenFHE_INCLUDE}/third-party/include)
include_directories(${OpenFHE_INCLUDE}/core)
include_directories(${OpenFHE_INCLUDE}/pke)
include_directories(${OpenFHE_INCLUDE}/binfhe)

link_directories(${OpenFHE_LIBDIR})
link_directories(${OPENMP_LIBRARIES})
if(BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries(${OpenFHE_STATIC_LIBRARIES})
else()
    set(CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS})
    link_libraries(${OpenFHE_SHARED_LIBRARIES})
endif()

# ===============================
# 查找 CUDA Toolkit
# ===============================
find_package(CUDAToolkit REQUIRED)

# ===============================
# 加入 FIDESlib 头文件目录
# 假设 FIDESlib 和 poppy 在同级目录
# ===============================
#include_directories(${CMAKE_CURRENT_LIST_DIR}/../FIDESlib/FIDESlib/include)

# ===============================
# 加入 CUDA 头文件目录
# ===============================
include_directories(${CUDAToolkit_INCLUDE_DIRS})

# ===============================
# POPPY 源码
# ===============================

# 定义 CUDA 源文件列表
set(CUDA_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/src/main.cu           # 改为 .cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils.cu    # 改为 .cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils_naive.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils_poppy_uj.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils_poppy_m.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils_export.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils_crypto.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils_graph.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils_stats.cu
    ${CMAKE_CURRENT_LIST_DIR}/test/calcul_orbits.cu
    ${CMAKE_CURRENT_LIST_DIR}/test/test.cu
    ${CMAKE_CURRENT_LIST_DIR}/test/expes.cu
    ${CMAKE_CURRENT_LIST_DIR}/test/op_compar.cu
    ${CMAKE_CURRENT_LIST_DIR}/test/rcm.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/clear.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/naive.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/poppy_uj.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/poppy_m.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/multi_DC_naive.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/multi_DC_poppy_uj.cu
    ${CMAKE_CURRENT_LIST_DIR}/src/multi_DC_poppy_m.cu
)
foreach(source_file ${CUDA_SOURCES})
    if(NOT EXISTS ${source_file})
        message(SEND_ERROR "Source file ${source_file} does not exist.")
    endif()
endforeach()

add_executable(calculpagerank ${CUDA_SOURCES})

# 显式设置源文件的语言属性为CUDA
set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)

# 检查CUDA编译器
if(CMAKE_CUDA_COMPILER)
    message(STATUS "CUDA compiler: ${CMAKE_CUDA_COMPILER}")
else()
    message(FATAL_ERROR "CUDA compiler not found")
endif()
# 设置目标属性，确保正确的编译和链接
target_compile_features(calculpagerank PRIVATE cxx_std_17)
target_include_directories(calculpagerank PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/include
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/third-party/include
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
    ${CMAKE_CURRENT_LIST_DIR}/../FIDESlib/FIDESlib/include
    ${CUDAToolkit_INCLUDE_DIRS}
)

# ===============================
# 链接库 - 使用现代 CMake 目标链接方式
# ===============================
target_link_libraries(calculpagerank
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/../nauty/nauty.a
        CUDA::cudart
)


if(METIS)
    link_libraries(metis)
    
    # 定义 graphpartitionner 的源文件
    set(GRAPH_PARTITION_SOURCES
        ${CMAKE_CURRENT_LIST_DIR}/test/graph_partitionner.cu  # 改为 .cu
        ${CMAKE_CURRENT_LIST_DIR}/src/utils/utils.cu
    )
    
    add_executable(graphpartitionner ${GRAPH_PARTITION_SOURCES})
    
    target_include_directories(graphpartitionner PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/include
        ${CUDAToolkit_INCLUDE_DIRS}
    )
    
    target_link_libraries(graphpartitionner
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/../nauty/nauty.a
            CUDA::cudart
    )
endif()

