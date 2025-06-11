cmake_minimum_required(VERSION 3.5.1)

project(LSH_OpenFHE_Tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置构建类型
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Building in ${CMAKE_BUILD_TYPE} mode")

#--------------------------------------------------------------------
# 设置 OpenFHE 路径
#--------------------------------------------------------------------
set(OPENFHE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)

# 验证 OpenFHE 根目录存在
if(NOT EXISTS "${OPENFHE_ROOT}")
    message(FATAL_ERROR "OpenFHE root directory not found: ${OPENFHE_ROOT}")
endif()

# 验证 OpenFHE 已经构建
if(NOT EXISTS "${OPENFHE_ROOT}/build")
    message(FATAL_ERROR "OpenFHE build directory not found. Please build OpenFHE first by running 'make' in ${OPENFHE_ROOT}")
endif()

if(NOT EXISTS "${OPENFHE_ROOT}/build/lib")
    message(FATAL_ERROR "OpenFHE libraries not found. Please ensure OpenFHE is properly built.")
endif()

message(STATUS "OpenFHE root directory: ${OPENFHE_ROOT}")
message(STATUS "OpenFHE build directory: ${OPENFHE_ROOT}/build")

#--------------------------------------------------------------------
# 查找 OpenFHE 库
#--------------------------------------------------------------------
# 首先尝试使用 find_package
list(APPEND CMAKE_PREFIX_PATH "${OPENFHE_ROOT}/build")
find_package(OpenFHE QUIET)

if(OpenFHE_FOUND)
    message(STATUS "Found OpenFHE package")
    message(STATUS "OpenFHE Version: ${OpenFHE_VERSION}")
    message(STATUS "OpenFHE include directories: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE libraries: ${OpenFHE_SHARED_LIBRARIES}")
else()
    message(STATUS "OpenFHE package not found, using manual configuration")
    
    # 手动设置包含目录
    include_directories(
        ${OPENFHE_ROOT}/src/core/include
        ${OPENFHE_ROOT}/src/pke/include
        ${OPENFHE_ROOT}/src/binfhe/include
        ${OPENFHE_ROOT}/third-party/include
        ${OPENFHE_ROOT}/third-party/cereal/include
        ${OPENFHE_ROOT}/build/src/core
        ${OPENFHE_ROOT}/build/src/pke
    )
    
    # 设置库目录
    link_directories(${OPENFHE_ROOT}/build/lib)
    
    # 检测库文件并设置库名称
    if(EXISTS "${OPENFHE_ROOT}/build/lib/libOPENFHEcore.so" OR EXISTS "${OPENFHE_ROOT}/build/lib/libOPENFHEcore.a")
        set(OpenFHE_LIBRARIES OPENFHEcore OPENFHEpke OPENFHEbinfhe)
        message(STATUS "Using library names: OPENFHEcore, OPENFHEpke, OPENFHEbinfhe")
    elseif(EXISTS "${OPENFHE_ROOT}/build/lib/libopenfhe-core.so" OR EXISTS "${OPENFHE_ROOT}/build/lib/libopenfhe-core.a")
        set(OpenFHE_LIBRARIES openfhe-core openfhe-pke openfhe-binfhe)
        message(STATUS "Using library names: openfhe-core, openfhe-pke, openfhe-binfhe")
    else()
        # 列出实际存在的库文件
        file(GLOB LIB_FILES "${OPENFHE_ROOT}/build/lib/lib*.so" "${OPENFHE_ROOT}/build/lib/lib*.a")
        message(STATUS "Available library files:")
        foreach(lib ${LIB_FILES})
            message(STATUS "  ${lib}")
        endforeach()
        message(FATAL_ERROR "Cannot find OpenFHE libraries in ${OPENFHE_ROOT}/build/lib/")
    endif()
endif()

#--------------------------------------------------------------------
# 查找线程库
#--------------------------------------------------------------------
find_package(Threads REQUIRED)

#--------------------------------------------------------------------
# 自动检测并编译所有测试文件（包含子目录）
#--------------------------------------------------------------------
# 在当前目录和 test 子目录中查找 .cpp 文件
file(GLOB TEST_SOURCES_CURRENT CONFIGURE_DEPENDS "*.cpp")
file(GLOB TEST_SOURCES_TEST CONFIGURE_DEPENDS "test/*.cpp")

# 合并所有测试源文件
set(TEST_SOURCES ${TEST_SOURCES_CURRENT} ${TEST_SOURCES_TEST})

if(NOT TEST_SOURCES)
    message(WARNING "No .cpp files found in current directory or test/ subdirectory")
else()
    message(STATUS "Found test files:")
    foreach(test_source ${TEST_SOURCES})
        message(STATUS "  ${test_source}")
    endforeach()
endif()

foreach(test_source ${TEST_SOURCES})
    # 获取文件名（不含扩展名）作为可执行文件名
    get_filename_component(test_name ${test_source} NAME_WE)
    
    # 创建可执行文件
    add_executable(${test_name} ${test_source})
    
    # 设置输出目录
    set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    
    if(OpenFHE_FOUND)
        # 使用 find_package 找到的配置
        target_link_libraries(${test_name} ${OpenFHE_SHARED_LIBRARIES})
        target_include_directories(${test_name} PRIVATE ${OpenFHE_INCLUDE})
    else()
        # 使用手动配置
        target_link_libraries(${test_name} ${OpenFHE_LIBRARIES})
    endif()
    
    # 链接线程库
    target_link_libraries(${test_name} Threads::Threads)
    
    # 设置编译器标志
    target_compile_options(${test_name} PRIVATE -Wall -Wextra -O3)
    
    message(STATUS "Added executable: ${test_name}")
endforeach()

# 创建一个总目标来编译所有测试
if(TEST_SOURCES)
    set(ALL_TEST_TARGETS "")
    foreach(test_source ${TEST_SOURCES})
        get_filename_component(test_name ${test_source} NAME_WE)
        list(APPEND ALL_TEST_TARGETS ${test_name})
    endforeach()
    
    add_custom_target(all_tests DEPENDS ${ALL_TEST_TARGETS})
    message(STATUS "Created target 'all_tests' for all test executables")
    message(STATUS "Usage:")
    message(STATUS "  make <test_name>  - Build specific test")
    message(STATUS "  make all_tests    - Build all tests")
    message(STATUS "  make              - Build all tests")
endif()

#--------------------------------------------------------------------
# 显示配置信息
#--------------------------------------------------------------------
message(STATUS "")
message(STATUS "=== Build Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Output directory: ${CMAKE_BINARY_DIR}/bin")
if(OpenFHE_FOUND)
    message(STATUS "OpenFHE: Found via find_package")
else()
    message(STATUS "OpenFHE: Manual configuration")
    message(STATUS "OpenFHE libraries: ${OpenFHE_LIBRARIES}")
endif()
message(STATUS "===========================")
message(STATUS "")