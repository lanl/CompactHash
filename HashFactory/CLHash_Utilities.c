/* Copyright 2013-14.  Los Alamos National Security, LLC. This material was produced
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC
 * for the U.S. Department of Energy. The U.S. Government has rights to use,
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS
 * ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 * to produce derivative works, such modified software should be clearly marked,
 * so as not to confuse it with the version available from LANL.   
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at 
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed
 * under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 *
 * Under this license, it is required to include a reference to this work. We
 * request that each derivative work contain a reference to LANL Copyright 
 * Disclosure C14043/LA-CC-14-003 so that this work's impact can be roughly
 * measured. In addition, it is requested that a modifier is included as in
 * the following example:
 *
 * //<Uses | improves on | modified from> LANL Copyright Disclosure C14043/LA-CC-14-003
 *
 * This is LANL Copyright Disclosure C14043/LA-CC-14-003
 */

/**
 * @file   CLHash_Utilities.c
 * @author Ahrens/Robey -- rewritten for distribution
 * @date   Sun Aug 4 2013 
 */

#define _XOPEN_SOURCE 500

#include <execinfo.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "CLHash_Utilities.h"


//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif

#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif

#else
int clGetPlatformIDs(uint nentries, int *platforms, uint *nplatforms) { return 0; }
int clGetDeviceIDs(int platform, int device_type, uint nentries, int *devices, uint *ndevices) { return 0; } 
int clGetPlatformInfo(int platform, int param_name, size_t size, void *value, size_t *size_ret) { return 0; } 
int clGetDeviceInfo(int device, int param_name, size_t size, void *value, size_t *size_ret) { return 0; } 
int clCreateContext(int *properties, uint ndevices, const int *devices, void *notify, void *user_data, int *errcode_ret) { return 0; }
int clCreateCommandQueue(int context, int device, int properties, int *errcode_ret) { return 0; }
int clCreateProgramWithSource(int context, uint count, const char **strings, const size_t *lengths, int *errcode_ret) { return 0; }
int clBuildProgram(int program, uint ndevices, const int *device_list, const char *options, void *notify, void *user_data) { return 0; }
int clGetProgramBuildInfo(int program, int devices, int param_name, size_t size, void *value, size_t *size_ret) { return 0; }

#define CL_SUCCESS                           0
#define CL_PROGRAM_LOG                       0
#define CL_PROGRAM_BUILD_LOG                 0
#define CL_DEVICE_TYPE_GPU                   0
#define CL_DEVICE_TYPE_CPU                   0
#define CL_DEVICE_TYPE_ACCELERATOR           0
#define CL_DEVICE_NAME                       0
#define CL_DEVICE_VENDOR                     0
#define CL_PLATFORM_NAME                     0
#define CL_PLATFORM_VENDOR                   0
#define CL_PLATFORM_VERSION                  0
#define CL_QUEUE_PROFILING_ENABLE            0

#define CL_DEVICE_NOT_FOUND                 -1
#define CL_DEVICE_NOT_AVAILABLE             -2
#define CL_COMPILER_NOT_AVAILABLE           -3
#define CL_MEM_OBJECT_ALLOCATION_FAILURE    -4
#define CL_OUT_OF_RESOURCES                 -5
#define CL_OUT_OF_HOST_MEMORY               -6
#define CL_PROFILING_INFO_NOT_AVAILABLE     -7
#define CL_MEM_COPY_OVERLAP                 -8
#define CL_IMAGE_FORMAT_MISMATCH            -9
#define CL_IMAGE_FORMAT_NOT_SUPPORTED       -10
#define CL_BUILD_PROGRAM_FAILURE            -11
#define CL_MAP_FAILURE                      -12
#define CL_INVALID_VALUE                    -30
#define CL_INVALID_DEVICE_TYPE              -31
#define CL_INVALID_PLATFORM                 -32
#define CL_INVALID_DEVICE                   -33
#define CL_INVALID_CONTEXT                  -34
#define CL_INVALID_QUEUE_PROPERTIES         -35
#define CL_INVALID_COMMAND_QUEUE            -36
#define CL_INVALID_HOST_PTR                 -37
#define CL_INVALID_MEM_OBJECT               -38
#define CL_INVALID_IMAGE_FORMAT_DESCRIPTOR  -39
#define CL_INVALID_IMAGE_SIZE               -40
#define CL_INVALID_SAMPLER                  -41
#define CL_INVALID_BINARY                   -42
#define CL_INVALID_BUILD_OPTIONS            -43
#define CL_INVALID_PROGRAM                  -44
#define CL_INVALID_PROGRAM_EXECUTABLE       -45
#define CL_INVALID_KERNEL_NAME              -46
#define CL_INVALID_KERNEL_DEFINITION        -47
#define CL_INVALID_KERNEL                   -48
#define CL_INVALID_ARG_INDEX                -49
#define CL_INVALID_ARG_VALUE                -50
#define CL_INVALID_ARG_SIZE                 -51
#define CL_INVALID_KERNEL_ARGS              -52
#define CL_INVALID_WORK_DIMENSION           -53
#define CL_INVALID_WORK_GROUP_SIZE          -54
#define CL_INVALID_WORK_ITEM_SIZE           -55
#define CL_INVALID_GLOBAL_OFFSET            -56
#define CL_INVALID_EVENT_WAIT_LIST          -57
#define CL_INVALID_EVENT                    -58
#define CL_INVALID_OPERATION                -59
#define CL_INVALID_GL_OBJECT                -60
#define CL_INVALID_BUFFER_SIZE              -61
#define CL_INVALID_MIP_LEVEL                -62
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

typedef unsigned int uint;

char *program_name = NULL;

void CLHash_Init(const char *program_name_in)
{
   program_name = strdup(program_name_in);
}

/* Find a GPU or CPU associated with the first available platform */
void CLHash_Utilities_CreateContext_p(cl_context *context, cl_command_queue *command_queue, const char *file , int line) {

   uint num_platforms;
   cl_platform_id *platforms;
   cl_device_id device;
   int err;

   /* Get all available platforms */
   err = clGetPlatformIDs(0, NULL, &num_platforms);
   CLHash_Utilities_HandleError(err, "CLHash_Utilities_CreateContext", "clGetPlatformIDs");
   platforms = (cl_platform_id *)malloc(num_platforms*sizeof(cl_platform_id));
   /* Identify a platform */
   err = clGetPlatformIDs(num_platforms, platforms, NULL);
   CLHash_Utilities_HandleError(err, "CLHash_Utilities_CreateContext", "clGetPlatformIDs");

   if (DEBUG == 1) {
     char info[1024];
     for (int iplatform=0; iplatform<num_platforms; iplatform++){
       printf("  Platform %d:\n",iplatform+1);

       //clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_PROFILE,   1024L,info,0);
       //printf("    CL_PLATFORM_PROFILE    : %s\n",info);

       clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_VERSION,   1024L,info,0);
       printf("    CL_PLATFORM_VERSION    : %s\n",info);

       clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_NAME,      1024L,info,0);
       printf("    CL_PLATFORM_NAME       : %s\n",info);

       clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_VENDOR,    1024L,info,0);
       printf("    CL_PLATFORM_VENDOR     : %s\n",info);

       //clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_EXTENSIONS,1024L,info,0);
       // printf("    CL_PLATFORM_EXTENSIONS : %s\n",info);
     }
   }

   /* Access a device */
   for (int iplatform=0; iplatform<num_platforms; iplatform++){
     err = clGetDeviceIDs(platforms[iplatform], CL_DEVICE_TYPE_GPU, 1, &device, NULL);
     if(err != CL_DEVICE_NOT_FOUND){
       break;
     }
   }
   if(err == CL_DEVICE_NOT_FOUND){
     for (int iplatform=0; iplatform<num_platforms; iplatform++){
       err = clGetDeviceIDs(platforms[iplatform], CL_DEVICE_TYPE_ACCELERATOR, 1, &device, NULL);
       if(err != CL_DEVICE_NOT_FOUND){
         break;
       }
     }
     if(err == CL_DEVICE_NOT_FOUND){
       for (int iplatform=0; iplatform<num_platforms; iplatform++){
         err = clGetDeviceIDs(platforms[iplatform], CL_DEVICE_TYPE_CPU, 1, &device, NULL);
           if(err != CL_DEVICE_NOT_FOUND){
           break;
         }
       }
     }
   }
   CLHash_Utilities_HandleError(err, "CLHash_Utilities_CreateContext", "clGetDeviceIDs");

   if (DEBUG == 1) {
     char info[1024];
   
     printf("\n\n");
     printf("  Device:\n");
     clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(info), info, NULL);
     printf("    CL_DEVICE_NAME         : %s\n",info);

     clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(info), info, NULL);
     printf("    CL_DEVICE_VENDOR       : %s\n",info);
   }

   *context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
   if(err != CL_SUCCESS) CLHash_Utilities_PrintError_p(err, "CLHash_Utilities_CreateContext", "clCreateContext", file, line);

   *command_queue = clCreateCommandQueue(*context, device, CL_QUEUE_PROFILING_ENABLE, &err);
   if(err != CL_SUCCESS) CLHash_Utilities_PrintError_p(err, "CLHash_Utilities_CreateContext", "clCreateCommandQueue", file, line);

   free(platforms);
}

cl_program CLHash_Utilities_BuildProgramString_p(cl_context ctx, cl_device_id device, const char *programString, const char *file , int line){
   /* Create program from string */
   cl_int err;
   size_t program_size, nReportSize;
   program_size = strlen(programString);
   cl_program program = clCreateProgramWithSource(ctx, 1, (const char **)&programString, &program_size, &err);
   if(err != CL_SUCCESS) CLHash_Utilities_PrintError(err, "CLHash_Utilities_CreateDevice", "clCreateProgram");

   /* Build program */
   char *BuildReport;
   err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
   if(err != CL_SUCCESS) {
      /* Find size of log and print to std output */
      clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &nReportSize);

      if (err != CL_SUCCESS) {
        switch (err){
          case CL_INVALID_DEVICE:
            printf("Invalid device in clProgramBuildInfo\n");
            break;
          case CL_INVALID_VALUE:
            printf("Invalid value in clProgramBuildInfo\n");
            break;
          case CL_INVALID_PROGRAM:
            printf("Invalid program in clProgramBuildInfo\n");
            break;
        }
      }

      BuildReport = (char*) malloc(nReportSize + 1);
      BuildReport[nReportSize] = '\0';
      clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, nReportSize + 1, BuildReport, NULL);
      printf("%s\n", BuildReport);
      free(BuildReport);
      exit(1);
   }

   return program;
}

/* Create program from a file and compile it */
cl_program CLHash_Utilities_BuildProgramFile_p(cl_context ctx, cl_device_id device, const char* fileName, const char *file , int line) {

   cl_program program;
   FILE *program_handle;
   char *program_buffer;
   size_t program_size;
   int err;

   /* Read program file and place content into buffer */
   program_handle = fopen(fileName, "r");
   if(program_handle == NULL) {
      perror("Couldn't find the program file");
      exit(1);
   }
   fseek(program_handle, 0, SEEK_END);
   program_size = ftell(program_handle);
   rewind(program_handle);
   program_buffer = (char*)malloc(program_size + 1);
   program_buffer[program_size] = '\0';
   fread(program_buffer, sizeof(char), program_size, program_handle);
   fclose(program_handle);
   program = CLHash_Utilities_BuildProgramString(ctx, device, program_buffer);
   free(program_buffer);
   return program;
}

void CLHash_Utilities_PrintError_p(cl_int err, const char *routine, const char *cl_routine, const char *file, int line){
  switch(err){
    case CL_DEVICE_NOT_FOUND:                //#define CL_DEVICE_NOT_FOUND                 -1
      printf("\nERROR: %s -- Device not found in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_DEVICE_NOT_AVAILABLE:            //#define CL_DEVICE_NOT_AVAILABLE             -2
      printf("\nERROR: %s -- Device not available in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_COMPILER_NOT_AVAILABLE:          //#define CL_COMPILER_NOT_AVAILABLE           -3
      printf("\nERROR: %s -- CL compiler not available failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:   //#define CL_MEM_OBJECT_ALLOCATION_FAILURE    -4
      printf("\nERROR: %s -- Mem object allocation failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_OUT_OF_RESOURCES:                //#define CL_OUT_OF_RESOURCES                 -5
      printf("\nERROR: %s -- Out of resources in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_OUT_OF_HOST_MEMORY:              //#define CL_OUT_OF_HOST_MEMORY               -6
      printf("\nERROR: %s -- Out of host memory in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_PROFILING_INFO_NOT_AVAILABLE:    //#define CL_PROFILING_INFO_NOT_AVAILABLE     -7
      printf("\nERROR: %s -- Profiling info not available in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_MEM_COPY_OVERLAP:                //#define CL_MEM_COPY_OVERLAP                 -8
      printf("\nERROR: %s -- Mem copy overlap in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_IMAGE_FORMAT_MISMATCH:           //#define CL_IMAGE_FORMAT_MISMATCH            -9
      printf("\nERROR: %s -- Image format mismatch in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:      //#define CL_IMAGE_FORMAT_NOT_SUPPORTED      -10
      printf("\nERROR: %s -- Image format not supported in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_BUILD_PROGRAM_FAILURE:           //#define CL_BUILD_PROGRAM_FAILURE           -11
      printf("\nERROR: %s -- Build program failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_MAP_FAILURE:                     //#define CL_MAP_FAILURE                     -12
      printf("\nERROR: %s -- Map failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_VALUE:                   //#define CL_INVALID_VALUE                   -30
      printf("\nERROR: %s -- Invalid value in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_DEVICE_TYPE:             //#define CL_INVALID_DEVICE_TYPE             -31
      printf("\nERROR: %s -- Invalid device type in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_PLATFORM:                //#define CL_INVALID_PLATFORM                -32
      printf("\nERROR: %s -- Invalid platform in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_DEVICE:                  //#define CL_INVALID_DEVICE                  -33
      printf("\nERROR: %s -- Invalid device in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_CONTEXT:                 //#define CL_INVALID_CONTEXT                 -34
      printf("\nERROR: %s -- Invalid context in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_QUEUE_PROPERTIES:        //#define CL_INVALID_QUEUE_PROPERTIES        -35
      printf("\nERROR: %s -- Invalid queue properties in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_COMMAND_QUEUE:           //#define CL_INVALID_COMMAND_QUEUE           -36
      printf("\nERROR: %s -- Invalid command queue in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_HOST_PTR:                //#define CL_INVALID_HOST_PTR                -37
      printf("\nERROR: %s -- Invalid host pointer in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_MEM_OBJECT:              //#define CL_INVALID_MEM_OBJECT              -38
      printf("\nERROR: %s -- Invalid memory object in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: //#define CL_INVALID_IMAGE_FORMAT_DESCRIPTOR -39
      printf("\nERROR: %s -- Invalid image format descriptor in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_IMAGE_SIZE:              //#define CL_INVALID_IMAGE_SIZE              -40
      printf("\nERROR: %s -- Invalid image size in %s at line %d in file %s\n", routine, cl_routine, line, file); 
      break;
    case CL_INVALID_SAMPLER:                 //#define CL_INVALID_SAMPLER                 -41
      printf("\nERROR: %s -- Invalid sampler in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_BINARY:                  //#define CL_INVALID_BINARY                  -42
      printf("\nERROR: %s -- Invalid binary in %s at line %d in file %s\n", routine, cl_routine, line, file);  
      break;
    case CL_INVALID_BUILD_OPTIONS:           //#define CL_INVALID_BUILD_OPTIONS           -43
      printf("\nERROR: %s -- Invalid build options in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_PROGRAM:                 //#define CL_INVALID_PROGRAM                 -44
      printf("\nERROR: %s -- Invalid program in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_PROGRAM_EXECUTABLE:      //#define CL_INVALID_PROGRAM_EXECUTABLE      -45
      printf("\nERROR: %s -- Invalid program executable in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_KERNEL_NAME:             //#define CL_INVALID_KERNEL_NAME             -46
      printf("\nERROR: %s -- Invalid kernel name in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_KERNEL_DEFINITION:       //#define CL_INVALID_KERNEL_DEFINITION       -47
      printf("\nERROR: %s -- Invalid kernel definition in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_KERNEL:                  //#define CL_INVALID_KERNEL                  -48
      printf("\nERROR: %s -- Invalid kernel in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_ARG_INDEX:               //#define CL_INVALID_ARG_INDEX               -49
      printf("\nERROR: %s -- Invalid arg index in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_ARG_VALUE:               //#define CL_INVALID_ARG_VALUE               -50
      printf("\nERROR: %s -- Invalid arg value in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_ARG_SIZE:                //#define CL_INVALID_ARG_SIZE                -51
      printf("\nERROR: %s -- Invalid arg size in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_KERNEL_ARGS:             //#define CL_INVALID_KERNEL_ARGS             -52
      printf("\nERROR: %s -- Invalid kernel args in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_WORK_DIMENSION:          //#define CL_INVALID_WORK_DIMENSION          -53
      printf("\nERROR: %s -- Invalid work dimension in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_WORK_GROUP_SIZE:         //#define CL_INVALID_WORK_GROUP_SIZE         -54
      printf("\nERROR: %s -- Invalid work group size in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_WORK_ITEM_SIZE:          //#define CL_INVALID_WORK_ITEM_SIZE          -55
      printf("\nERROR: %s -- Invalid work item size in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_GLOBAL_OFFSET:           //#define CL_INVALID_GLOBAL_OFFSET           -56
      printf("\nERROR: %s -- Invalid global offset in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_EVENT_WAIT_LIST:         //#define CL_INVALID_EVENT_WAIT_LIST         -57
      printf("\nERROR: %s -- Invalid event wait list in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_EVENT:                   //#define CL_INVALID_EVENT                   -58
      printf("\nERROR: %s -- Invalid event in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_OPERATION:               //#define CL_INVALID_OPERATION               -59
      printf("\nERROR: %s -- Invalid operation in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_GL_OBJECT:               //#define CL_INVALID_GL_OBJECT               -60
      printf("\nERROR: %s -- Invalid GL object in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_BUFFER_SIZE:             //#define CL_INVALID_BUFFER_SIZE             -61
      printf("\nERROR: %s -- Invalid buffer size in %s at line %d in file %s\n", routine, cl_routine, line, file);
      break;
    case CL_INVALID_MIP_LEVEL:               //#define CL_INVALID_MIP_LEVEL               -62
      break;

    default:
      printf("\nERROR: %s -- %d in %s at line %d in file %s\n", routine, err, cl_routine, line, file);
      break;
  }
  void* callstack[128];
  int frames = backtrace(callstack, 128);
  if (frames > 2) {
#ifdef HAVE_ADDR2LINE
     char hex_address[21];
     char command_string[80];
#endif
     char** strs = backtrace_symbols(callstack, frames);
     fprintf(stderr,"\n  =============== Backtrace ===============\n");
     for (int i = 1; i < frames-1; ++i) {
       fprintf(stderr,"   %s    \t", strs[i]);
#ifdef HAVE_ADDR2LINE
       if (program_name != NULL) {
          sscanf(strs[i],"%*s [%[^]]s]",hex_address);
       //printf("DEBUG addr2line -e clamr -f -s -i -p %s\n",hex_address);
          //sprintf(command_string,"addr2line -e %s -f -s -i %s | sed -e \'1,1s\?\\n\?\?' \n",program_name,hex_address);
          sprintf(command_string,"addr2line -e %s -f -s -i %s",program_name,hex_address);
          //printf("%s",command_string);
          system(command_string);
       }
       // on mac, need to install binutils using macports "port install binutils"
#endif
       fprintf(stderr,"\n");
    }
    fprintf(stderr,"  =============== Backtrace ===============\n\n");
    free(strs);
  }
  exit(-1);
}

const char* CLHash_Utilities_HandleError_p(cl_int err, const char *routine, const char *cl_routine, const char *file, int line){
  if(err != CL_SUCCESS) CLHash_Utilities_PrintError_p(err, routine, cl_routine, file, line);
  return(NULL);
}
