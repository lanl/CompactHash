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
 * @file   CLHash_Utilities.h
 * @author Ahrens/Robey -- rewritten for distribution
 * @date   Sun Aug 4 2013 
 */
#ifndef CLHASH_UTILITIES_H
#define CLHASH_UTILITIES_H

#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#else
typedef unsigned int uint;
typedef int cl_int;
typedef int cl_context;
typedef int cl_program;
typedef int cl_command_queue;
typedef int cl_device_id;
typedef int cl_platform_id;

int clGetPlatformIDs(uint nentries, int *platforms, uint *nplatforms);
int clGetDeviceIDs(int platform, int device_type, uint nentries, int *devices, uint *ndevices);
int clGetPlatformInfo(int platform, int param_name, size_t size, void *value, size_t *size_ret);
int clGetDeviceInfo(int device, int param_name, size_t size, void *value, size_t *size_ret);
int clCreateContext(int *properties, uint ndevices, const int *devices, void *notify, void *user_data, int *errcode_ret);
int clCreateCommandQueue(int context, int device, int properties, int *errcode_ret);
int clCreateProgramWithSource(int context, uint count, const char **strings, const size_t *lengths, int *errcode_ret);
int clBuildProgram(int program, uint ndevices, const int *device_list, const char *options, void *notify, void *user_data);
int clGetProgramBuildInfo(int program, int devices, int param_name, size_t size, void *value, size_t *size_ret);
#endif

#define CLHash_Utilities_CreateContext(   context, command_queue) \
      ( CLHash_Utilities_CreateContext_p( context, command_queue, __FILE__, __LINE__) )
#define CLHash_Utilities_BuildProgramString(   context, device, programString) \
      ( CLHash_Utilities_BuildProgramString_p( context, device, programString, __FILE__, __LINE__) )
#define CLHash_Utilities_BuildProgramFile(   context, device, fileName) \
      ( CLHash_Utilities_BuildProgramFile_p( context, device, fileName, __FILE__, __LINE__) )

#define CLHash_Utilities_HandleError(   ierr, routine, cl_routine) \
      ( CLHash_Utilities_HandleError_p( ierr, routine, cl_routine, __FILE__, __LINE__) )
#define CLHash_Utilities_PrintError(   ierr, routine, cl_routine) \
      ( CLHash_Utilities_PrintError_p( ierr, routine, cl_routine, __FILE__, __LINE__) )

/* Initialize Hash lib with program name for error reporting */
void CLHash_Init(const char *program_name_in);

/* Hardware detection */
void CLHash_Utilities_CreateContext_p(cl_context *context, cl_command_queue *command_queue, const char *file, int line);

/* Build Programs */
cl_program CLHash_Utilities_BuildProgramString_p(cl_context contex, cl_device_id device, const char *programString, const char *file, int line);
cl_program CLHash_Utilities_BuildProgramFile_p(cl_context contex, cl_device_id device, const char* fileName, const char *file, int line);

/* Error Handling */
const char* CLHash_Utilities_HandleError_p(cl_int err, const char *routine, const char *cl_routine, const char *file, int line);
void CLHash_Utilities_PrintError_p(cl_int err, const char *routine, const char *cl_routine, const char *file, int line);

#endif
