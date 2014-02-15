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

#include "hasholdlib.h"
#include "hasholdlib_kern.inc"
#include "hashold_kern.inc"

static int reportLevel = 0;

char *get_hash_kernel_source_string(void){
   return (char*) hashold_kern_source;
}

cl_kernel init_kernel;
void hash_lib_init(cl_context context){

   cl_int error;

   cl_program program =  clCreateProgramWithSource(context, 1, (const char **)&hasholdlib_kern_source, NULL, &error);
   if (error != CL_SUCCESS){
       printf("clCreateProgramWithSource returned an error %d at line %d in file %s\n", error,__LINE__,__FILE__);
   }

   size_t nReportSize;
   char* BuildReport;
   #ifdef HAVE_CL_DOUBLE
     error = clBuildProgram(program, 0, NULL, "-DHAVE_CL_DOUBLE", NULL, NULL);
   #else
     error = clBuildProgram(program, 0, NULL, "-DNO_CL_DOUBLE -cl-single-precision-constant", NULL, NULL);
   #endif
    if (error != CL_SUCCESS){
       printf("clBuildProgram returned an error %d at line %d in file %s\n", error,__LINE__,__FILE__);

       cl_device_id device;

       clGetContextInfo(context, CL_CONTEXT_DEVICES, 1, &device, NULL);
       if (context == NULL){
          printf("EZCL_DEVTYPE_INIT: Failed to find device and setup context in file %s at line %d\n", __FILE__, __LINE__);
          exit(-1); // No device is available, something is wrong 
       }

       error = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &nReportSize);
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

       error = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, nReportSize, BuildReport, NULL);
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


   init_kernel = clCreateKernel(program, "init_kern", &error);
   if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
}

cl_mem hash_init (int hash_size, cl_context context, cl_command_queue queue, long *gpu_time)
{
   cl_int error;
   const int TILE_SIZE = 128;

   long gpu_time_start, gpu_time_end;

   cl_mem hash_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, hash_size*sizeof(int), NULL, &error);
   if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  //init to -1

   error = clSetKernelArg(init_kernel, 0, sizeof(cl_uint), &hash_size);
   if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
   error = clSetKernelArg(init_kernel, 1, sizeof(cl_mem), (void*)&hash_buffer);
   if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

   size_t global_work_size;
   size_t local_work_size;

   local_work_size = TILE_SIZE;
   global_work_size = ((hash_size+local_work_size-1)/local_work_size)*local_work_size;

   cl_event hash_init_event;

   error = clEnqueueNDRangeKernel(queue, init_kernel, 1, 0, &global_work_size, &local_work_size, 0, NULL, &hash_init_event);
   if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

   clWaitForEvents(1,&hash_init_event);

   clGetEventProfilingInfo(hash_init_event, CL_PROFILING_COMMAND_START, sizeof(gpu_time_start), &gpu_time_start, NULL);
   clGetEventProfilingInfo(hash_init_event, CL_PROFILING_COMMAND_END, sizeof(gpu_time_end), &gpu_time_end, NULL);


   *gpu_time = gpu_time_end - gpu_time_start;

   clReleaseEvent(hash_init_event);

   //if (DETAILED_TIMING) printf("\n\tinit %.6lf,", (double)(gpu_time_end - gpu_time_start)*1.0e-9);

   return(hash_buffer);
}
