/* Copyright 2013-19.  Triad National Security, LLC. This material was produced
 * under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 * Laboratory (LANL), which is operated by Triad National Security, LLC
 * for the U.S. Department of Energy. The U.S. Government has rights to use,
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 * TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
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
 *  Authors: Bob Robey         XCP-2   brobey@lanl.gov
 *           Peter Ahrens      XCP-2   peter.ahrens@lanl.gov, ptrahrens@gmail.com
 *           Sara Hartse       XCP-2   sara@lanl.gov, 
 *           Rebecka Tumblin           rtumblin@lanl.gov, rebeckatumblin@gmail.com
 */

/* Modified from LANL Copyright Disclosure C13002/LA-CC-12-022 */

#ifndef __GPU_H__
#define __GPU_H__

//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif

#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#endif

enum device_vendor_type {
   VENDOR_UNKNOWN = 0,
   NVIDIA,
   TESLA_K20Xm,
   TESLA_M2090,
   TESLA_2050_2070,
   GeFORCE_GT_650M,
   ATI,
   MIC,
};

#ifdef HAVE_OPENCL
cl_kernel interpolate_kernel;

void GPUInit(cl_context *context, cl_command_queue *queue, int *device_type, cl_program *program, char *filename, char *addlibsource);
int device_double_support(cl_device_id device);
void device_info(cl_device_id device);
#endif

#endif // __GPU_H__
