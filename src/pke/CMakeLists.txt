#
# CMakeLists.txt for PKE library
#

# all files named *.cpp are compiled to form the library
file(GLOB_RECURSE PKE_SRC_FILES CONFIGURE_DEPENDS lib/*.cpp)

include_directories(${CORE_INCLUDE_DIRS})
list(APPEND PKE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list(APPEND PKE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib")
include_directories(${PKE_INCLUDE_DIRS})

set(PKE_VERSION_MAJOR ${OPENFHE_VERSION_MAJOR})
set(PKE_VERSION_MINOR ${OPENFHE_VERSION_MINOR})
set(PKE_VERSION_PATCH ${OPENFHE_VERSION_PATCH})
set(PKE_VERSION ${PKE_VERSION_MAJOR}.${PKE_VERSION_MINOR}.${PKE_VERSION_PATCH})

add_library(pkeobj OBJECT ${PKE_SRC_FILES})
set_property(TARGET pkeobj PROPERTY POSITION_INDEPENDENT_CODE 1)

if(BUILD_SHARED)
    add_dependencies(pkeobj OPENFHEcore OPENFHEbinfhe)
    add_library(OPENFHEpke SHARED $<TARGET_OBJECTS:pkeobj>)
    set_property(TARGET OPENFHEpke PROPERTY VERSION ${PKE_VERSION})
    set_property(TARGET OPENFHEpke PROPERTY SOVERSION ${PKE_VERSION_MAJOR})
    set_property(TARGET OPENFHEpke PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    install(TARGETS OPENFHEpke EXPORT OpenFHETargets DESTINATION lib)
endif()

if(BUILD_STATIC)
    add_dependencies(pkeobj OPENFHEcore_static OPENFHEbinfhe_static)
    add_library(OPENFHEpke_static STATIC $<TARGET_OBJECTS:pkeobj>)
    set_property(TARGET OPENFHEpke_static PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    install(TARGETS OPENFHEpke_static EXPORT OpenFHETargets DESTINATION lib)
endif()

install(DIRECTORY include/ DESTINATION include/openfhe/pke)
install(DIRECTORY unittest/utils/ DESTINATION include/openfhe/pke/unittest/utils FILES_MATCHING PATTERN "*.h")

add_custom_target(allpke)

if(BUILD_SHARED)
    set(PKELIBS PUBLIC OPENFHEpke PUBLIC OPENFHEcore PUBLIC OPENFHEbinfhe ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS})
    target_link_libraries(OPENFHEpke PUBLIC OPENFHEcore PUBLIC OPENFHEbinfhe ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS} ${ADDITIONAL_LIBS})
    add_dependencies(allpke OPENFHEpke)
endif()

if(BUILD_STATIC)
    set(PKELIBS ${PKELIBS} PUBLIC OPENFHEpke_static PUBLIC OPENFHEcore_static PUBLIC OPENFHEbinfhe_static ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS})
    target_link_libraries(OPENFHEpke_static PUBLIC OPENFHEcore_static PUBLIC OPENFHEbinfhe_static ${THIRDPARTYSTATICLIBS} ${OpenMP_CXX_FLAGS} ${ADDITIONAL_LIBS})
    add_dependencies(allpke OPENFHEpke_static)
endif()

if(BUILD_UNITTESTS)
    if("${NATIVE_SIZE}" EQUAL 32)
        message("**** pke_tests are not linked for NATIVE_SIZE=32")
    else()
        file(GLOB_RECURSE PKE_TEST_SRC_FILES CONFIGURE_DEPENDS unittest/*.cpp)
    endif()
    add_executable(pke_tests ${PKE_TEST_SRC_FILES} ${UNITTESTMAIN})
    set_property(TARGET pke_tests PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/unittest)
    target_include_directories(pke_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/unittest/utils")
    target_link_libraries(pke_tests ${PKELIBS} ${ADDITIONAL_LIBS})
    if(NOT ${WITH_OPENMP} AND NOT EMSCRIPTEN)
        target_link_libraries(pke_tests PRIVATE Threads::Threads)
    endif()
    add_dependencies(allpke pke_tests)
    add_custom_command(OUTPUT runpketests WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND ${CMAKE_BINARY_DIR}/unittest/pke_tests)
    add_custom_target(testpke DEPENDS pke_tests runpketests)
endif()

set(PKEAPPS "")
if(BUILD_EXAMPLES)
    file(GLOB PKE_EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
    foreach(app ${PKE_EXAMPLES_SRC_FILES})
        get_filename_component(exe ${app} NAME_WE)
        if(${exe} STREQUAL "scheme-switching-serial")
            # add schemeswitching-data-serializer.cpp to link line and
            # include the unittest/utils directory where schemeswitching-data-serializer.h is
            add_executable(${exe} ${app} "${CMAKE_CURRENT_SOURCE_DIR}/unittest/utils/schemeswitching-data-serializer.cpp")
            target_include_directories(${exe} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/unittest/utils")
        else()
            add_executable(${exe} ${app})
        endif()
        set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/examples/pke)
        set(PKEAPPS ${PKEAPPS} ${exe})
        target_link_libraries(${exe} ${PKELIBS} ${ADDITIONAL_LIBS})
    endforeach()
    add_custom_target(allpkeexamples)
    add_dependencies(allpkeexamples ${PKEAPPS})
    add_dependencies(allpke allpkeexamples)
endif()


###########################################################################################
# Add workplace executables -- By SeehowLi
set(PKEAPPS_WORKPLACE "")
if(BUILD_WORKPLACE)
    message(STATUS "Building PKE workplace programs")
    
    # 递归查找workplace目录及其子目录下的所有.cpp文件
    file(GLOB_RECURSE PKE_WORKPLACE_SRC_FILES CONFIGURE_DEPENDS 
        ${CMAKE_CURRENT_SOURCE_DIR}/workplace/*.cpp
    )
    
    message(STATUS "Found ${CMAKE_CURRENT_SOURCE_DIR}/workplace/*.cpp files:")
    foreach(file ${PKE_WORKPLACE_SRC_FILES})
        message(STATUS "  - ${file}")
    endforeach()
    
    # 为每个.cpp文件创建可执行文件
    foreach(app ${PKE_WORKPLACE_SRC_FILES})
        get_filename_component(exe ${app} NAME_WE)
        
        # 直接使用文件名，不加前缀
        set(exe_name ${exe})
        
        # 检查是否已经存在同名目标
        if(TARGET ${exe_name})
            message(WARNING "Target ${exe_name} already exists! Skipping ${app}")
            continue()
        endif()
        
        message(STATUS "Creating executable: ${exe_name}")
        add_executable(${exe_name} ${app})
        set_property(TARGET ${exe_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/workplace)
        set(PKEAPPS_WORKPLACE ${PKEAPPS_WORKPLACE} ${exe_name})
        target_link_libraries(${exe_name} ${PKELIBS} ${ADDITIONAL_LIBS})
        
        # 添加必要的包含路径
        get_filename_component(app_dir ${app} DIRECTORY)
        target_include_directories(${exe_name} PRIVATE 
            "${CMAKE_CURRENT_SOURCE_DIR}/workplace"
            "${app_dir}"
        )
    endforeach()
    
    # 创建 allpkeworkplace 目标
    add_custom_target(allpkeworkplace)
    if(PKEAPPS_WORKPLACE)
        add_dependencies(allpkeworkplace ${PKEAPPS_WORKPLACE})
    endif()
    add_dependencies(allpke allpkeworkplace)
endif()
###########################################################################################

set(PKEEXTRAS "")
if(BUILD_EXTRAS)
    file(GLOB PKE_EXTRAS_SRC_FILES CONFIGURE_DEPENDS extras/*.cpp)
    foreach(app ${PKE_EXTRAS_SRC_FILES})
        get_filename_component(exe ${app} NAME_WE)
        add_executable(${exe} ${app})
        set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/extras/pke)
        set(PKEEXTRAS ${PKEEXTRAS} ${exe})
        target_link_libraries(${exe} ${PKELIBS} ${ADDITIONAL_LIBS})
    endforeach()
    add_custom_target(allpkeextras)
    add_dependencies(allpkeextras ${PKEEXTRAS})
    add_dependencies(allpke allpkeextras)
endif()

add_custom_command(OUTPUT pkeinfocmd COMMAND echo Builds OPENFHEpke and these apps: ${PKEAPPS})
add_custom_target(pkeinfo DEPENDS pkeinfocmd)

# Collect compile definitions and pass them upward
if(BUILD_SHARED)
    get_target_property(_compile_defs OPENFHEpke COMPILE_DEFINITIONS)
    set(_pal_pke_compile_defs ${_compile_defs} PARENT_SCOPE)
endif()

if(BUILD_STATIC)
    get_target_property(_compile_defs_static OPENFHEpke_static COMPILE_DEFINITIONS)
    set(_pal_pke_compile_defs_static ${_compile_defs_static} PARENT_SCOPE)
endif()
