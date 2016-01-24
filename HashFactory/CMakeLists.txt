set(libHashFactory_LIB_SRCS CLHash_Utilities.c  HashFactory.c CLHash_Utilities.h HashFactory.h)

set(INDENT indent)
set(CPREPROCESSOR cpp)

#COMMAND ${CPREPROCESSOR} -CC -x c -imacros HashFactory.hm HashFactory.hp | ${INDENT} -linux -brf -o HashFactory.h

########### generate sources target ###############
add_custom_target(HashFactory.h ALL 
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.hm ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.hp
                          ${CMAKE_CURRENT_SOURCE_DIR}/perfectHash.hm ${CMAKE_CURRENT_SOURCE_DIR}/compactHash.hm
                  COMMAND cpp -P -CC -x c -imacros HashFactory.hm HashFactory.hp | indent -linux -brf -o HashFactory.h)

add_custom_target(HashFactory.c ALL
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.cm ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.cp 
                          ${CMAKE_CURRENT_SOURCE_DIR}/perfectHash.cm ${CMAKE_CURRENT_SOURCE_DIR}/compactHash.cm
                          HashFactory.cl HashFactory.h
                  COMMAND cpp -P -CC -x c -imacros HashFactory.cm HashFactory.cp | indent -linux -brf |
                     sed -e '/OMP_PRAGMA_PARALLEL_FOR/s/[[:space:]]*OMP_PRAGMA_PARALLEL_FOR\;/\#pragma omp parallel for/' > HashFactory.c &&
                      ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl HashFactory.cl >> ${CMAKE_CURRENT_BINARY_DIR}/HashFactory.c)

add_custom_target(HashFactory.cl ALL
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.clp ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.cm
                          ${CMAKE_CURRENT_SOURCE_DIR}/perfectHash.cm ${CMAKE_CURRENT_SOURCE_DIR}/compactHash.cm
                  COMMAND cpp -P -CC -x c -imacros HashFactory.cm HashFactory.clp | indent -linux -brf -o HashFactory.cl)

########### embed source target ###############
add_custom_target(HashFactory_source.inc ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/HashFactory.cl
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl
                     ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.cl > ${CMAKE_CURRENT_BINARY_DIR}/HashFactory_source.inc
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/HashFactory.cl ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl)

########## HashFactory target ##############
add_library(HashFactory STATIC ${libHashFactory_LIB_SRCS})
if (OPENCL_FOUND)
   set_target_properties(HashFactory PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
   add_dependencies(HashFactory HashFactory_source.inc)
endif (OPENCL_FOUND)

########## HashFactory_openmp target ##############
if (OPENMP_FOUND)
   add_library(HashFactory_openmp STATIC ${libHashFactory_LIB_SRCS})
   set_target_properties(HashFactory_openmp PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}")
   if (OPENCL_FOUND)
      set_target_properties(HashFactory_openmp PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
   endif (OPENCL_FOUND)
endif (OPENMP_FOUND)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
"HashFactory_source.inc;")
