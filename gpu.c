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
 * measured. 
 *
 * This is LANL Copyright Disclosure C14043/LA-CC-14-003
 */

/*
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Peter Ahrens            peter.ahrens@lanl.gov, ptrahrens@gmail.com
 *           Sara Hartse             sara@lanl.gov, sara.hartse@gmail.com
 *           Rebecka Tumblin         rtumblin@lanl.gov, rebeckatumblin@gmail.com
 */

/* Modified from LANL Copyright Disclosure C13002/LA-CC-12-022 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include "gpu.h"

#ifdef HAVE_CL_DOUBLE
typedef double real;
typedef cl_double cl_real;
typedef cl_double4 cl_real4;
#else
typedef float real;
typedef cl_float cl_real;
typedef cl_float4 cl_real4;
#endif

void GPUInit(cl_context *context, cl_command_queue *queue, int *is_nvidia, cl_program *program, char *filename, char *addlibsource) {
  printf("============== Hardware detection report ===============\n");

  cl_platform_id* platform;
  cl_uint num_platforms;
  cl_uint num_devices;
  cl_device_id* device;
  //cl_program program;
  cl_int error = 0;
  
  error = clGetPlatformIDs(0, NULL, &num_platforms);
  if(error != CL_SUCCESS) {
    printf("Error getting number of platforms\n");
    exit(error);
  }
  platform = (cl_platform_id*)malloc(num_platforms*sizeof(cl_platform_id));
  error = clGetPlatformIDs(num_platforms, platform, 0);
  if(error != CL_SUCCESS) {
    printf("Error getting platform id\n");
    exit(error);
  }
  char info[1024];
  clGetPlatformInfo(platform[0],CL_PLATFORM_VENDOR, 1024, info, 0);

  error = clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
  if (error != CL_SUCCESS) error = clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_ACCELERATOR, 0, NULL, &num_devices);
  if(error != CL_SUCCESS) {
    printf("Error getting device ids1\n");
    exit(error);
  }
  device = (cl_device_id*)malloc(num_devices*sizeof(cl_device_id));
  
  error = clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_GPU, num_devices, device, 0);
  if(error != CL_SUCCESS) error = clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_ACCELERATOR, num_devices, device, 0);
  if(error != CL_SUCCESS) {
    printf("Error getting device ids\n");
    exit(error);
  }
  
 // char platform_info[1024];
  for (int iplatform=0; iplatform<num_platforms; iplatform++){
    printf("  Platform %d:\n",iplatform+1);

    //clGetPlatformInfo(platform[iplatform],CL_PLATFORM_PROFILE,   1024L,info,0);
    //printf("    CL_PLATFORM_PROFILE    : %s\n",info);

    clGetPlatformInfo(platform[iplatform],CL_PLATFORM_VERSION,   1024L,info,0);
    printf("    CL_PLATFORM_VERSION    : %s\n",info);
     
    clGetPlatformInfo(platform[iplatform],CL_PLATFORM_NAME,      1024L,info,0);
    printf("    CL_PLATFORM_NAME       : %s\n",info);
     
    clGetPlatformInfo(platform[iplatform],CL_PLATFORM_VENDOR,    1024L,info,0);
    printf("    CL_PLATFORM_VENDOR     : %s\n",info);
     
    //clGetPlatformInfo(platform[iplatform],CL_PLATFORM_EXTENSIONS,1024L,info,0);
    //printf("    CL_PLATFORM_EXTENSIONS : %s\n",info);
  }
   
  *is_nvidia = 0;
  printf("\n");
  printf("  Device:\n");
  clGetDeviceInfo(device[0], CL_DEVICE_NAME, sizeof(info), &info, NULL);
  printf("    CL_DEVICE_NAME         : %s\n",info);
  
  clGetDeviceInfo(device[0], CL_DEVICE_VENDOR, sizeof(info), &info, NULL);
  printf("    CL_DEVICE_VENDOR       : %s\n",info);
  
  if (! strncmp(info,"NVIDIA",6) ) *is_nvidia = 1;
  //printf("DEBUG -- device vendor is |%s|, is_nvidia %d\n",info,*is_nvidia);

  *context = clCreateContext(0, 1, device, NULL, NULL, &error);
  if(error != CL_SUCCESS) {
    printf("Error creating context\n");
    exit(error);
  }
  *queue = clCreateCommandQueue(*context, device[0], CL_QUEUE_PROFILING_ENABLE, &error);
  if(error != CL_SUCCESS) {
    printf("Error creating command queue\n");
    exit(error);
  }
  
  printf("============== Hardware detection report ===============\n");

  // Load the kernel source code into the array source
  struct stat statbuf;
  FILE *fh;
  char *source;
  
  fh = fopen(filename, "r");
  if (!fh) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(-1);
  }
  stat(filename, &statbuf);
  
  int string_pos = 0;
  //int line_count = 0;
  int addlibsource_length = strlen(addlibsource);
  //printf("DEBUG -- length is addlibsource_length %d\n",addlibsource_length);
  if (addlibsource != NULL && addlibsource_length > 0){
    source = (char*)malloc(statbuf.st_size + 1 + addlibsource_length+1);
    strncpy(source, addlibsource, addlibsource_length);
    string_pos = addlibsource_length;
    source[string_pos] = '\n';
    string_pos++;
    //const char *src_string = source;
    //for (uint sindex = 0; sindex < addlibsource_length; src_string++, sindex++){
    //   if (! strncmp(src_string, '\n', 1)) line_count++;
    //}
    //printf("Offset line error in compile by %d lines\n",line_count);
  } else {
    source = (char*)malloc(statbuf.st_size + 1);
  }
  if( fread(&source[string_pos], statbuf.st_size, 1, fh) != 1) {
    printf("Problem reading program source file\n");
  }
  source[string_pos+statbuf.st_size] = '\0';
  //printf("DEBUG string is %s\n",source);
  fclose( fh );
  
  *program = clCreateProgramWithSource(*context, 1, (const char**) &source, NULL, &error);
  if (error != CL_SUCCESS){
    printf("clCreateProgramWithSource returned an error %d at line %d in file %s\n", error,__LINE__,__FILE__);
  }
  //printf("%d %s\n", (int)statbuf.st_size, source);
  
  size_t nReportSize;
  char* BuildReport;
  
  #ifdef HAVE_CL_DOUBLE
    if (*is_nvidia) {
      error = clBuildProgram(*program, 1, device, "-DHAVE_CL_DOUBLE -DIS_NVIDIA", NULL, NULL);
    } else {
      error = clBuildProgram(*program, 1, device, "-DHAVE_CL_DOUBLE", NULL, NULL);
    }
  #else
    if (*is_nvidia) {
       error = clBuildProgram(*program, 1, device, "-DNO_CL_DOUBLE -DIS_NVIDIA -cl-single-precision-constant", NULL, NULL);
    } else {
       error = clBuildProgram(*program, 1, device, "-DNO_CL_DOUBLE -cl-single-precision-constant", NULL, NULL);
    }
  #endif
  if (error != CL_SUCCESS){
    printf("clBuildProgram returned an error %d at line %d in file %s\n", error,__LINE__,__FILE__);
    error = clGetProgramBuildInfo(*program, device[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &nReportSize);
    if (error != CL_SUCCESS) {
      switch (error){
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
    
    BuildReport = (char *)malloc(nReportSize);
    
    error = clGetProgramBuildInfo(*program, device[0], CL_PROGRAM_BUILD_LOG, nReportSize, BuildReport, NULL);
    if (error != CL_SUCCESS) {
      switch (error){
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
    printf("%s\n", BuildReport);
  }
}
