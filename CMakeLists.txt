cmake_minimum_required(VERSION 2.6)
project (CompactHash)
# The version number.
set (COMPACTHASH_VERSION_MAJOR 1)
set (COMPACTHASH_VERSION_MINOR 2.0)
set (COMPACTHASH_VERSION "${COMPACTHASH_VERSION_MAJOR}.${COMPACTHASH_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

if (NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE RelWithDebInfo)
   set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -std=c99")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -vec-report3 -std=c99")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
endif()

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELDEBINFO} -g -O3")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -vec-report3")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
endif()

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELDEBINFO} -g -O3")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELDEBINFO} -g -O3")

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

#message("CMAKE_C_FLAGS ${CMAKE_C_FLAGS}")

#if (NOT MIC_NATIVE)
   find_package(OpenCL)
   #  OPENCL_FOUND        : TRUE if found, FALSE otherwise
   #  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
   #  OPENCL_LIBRARIES    : The libraries to link against
   if (OPENCL_FOUND)
      include_directories(${OPENCL_INCLUDE_DIRS})
   endif (OPENCL_FOUND)
#endif (NOT MIC_NATIVE)

#find_package(MPI)
#if (MPI_FOUND)
#   include_directories(${MPI_INCLUDE_PATH})
#endif (MPI_FOUND)

add_subdirectory(HashFactory)
add_subdirectory(kdtree)
add_subdirectory(timer)
add_subdirectory(simplehash)

#### global settings #######

set (NEIGH2D_SRCS neigh2d.c gpu.c)
set (NEIGH2D_HDRS gpu.h)

########### embed source target ############

#add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/neigh2d_kern.inc
#                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl 
#                      ${CMAKE_CURRENT_SOURCE_DIR}/neigh2d_kern.cl > 
#                      ${CMAKE_CURRENT_BINARY_DIR}/neigh2d_kern.inc
#                   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/neigh2d_kern.cl ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl)
#add_custom_target(neigh2d_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/neigh2d_kern.inc)

########### neigh2d target ###############

add_executable(neigh2d ${NEIGH2D_SRCS} ${NEIGH2D_HDRS})
target_link_libraries(neigh2d HashFactory simplehash kdtree timer m)
if (OPENCL_FOUND)
   set_target_properties(neigh2d PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
   target_link_libraries(neigh2d ${OPENCL_LIBRARIES})
#   add_dependencies(neigh2d neigh2d_kernel_source)
endif(OPENCL_FOUND)

########### clean 
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "")
