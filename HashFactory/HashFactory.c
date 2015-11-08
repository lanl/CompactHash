
/* Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* We do support the IEC 559 math functionality, real and complex.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
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
 * @file   HashFactory.c
 * @author Peter Ahrens
 * @date   Thu Jun 6 2013 
 */
//
#ifdef _OPENMP
#include <omp.h>
#endif
#include "HashFactory.h"
#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#else
typedef int cl_kernel;
#define CL_CONTEXT_DEVICES 0
#define CL_MEM_READ_WRITE 0
#define CL_TRUE 1
#define CL_MEM_READ_ONLY 0
#define CL_MEM_COPY_HOST_PTR 0
#define CL_MEM_WRITE_ONLY 0
#define CL_KERNEL_WORK_GROUP_SIZE 128
int clRetainContext(int context) {
	return context;
}
int clRetainCommandQueue(int command_queue) {
	return command_queue;
}
int clGetContextInfo(int context, int param, size_t size, void *value,
		     size_t * size_ret) {
	return 0;
}
int clReleaseContext(int context) {
	return context;
}
int clReleaseCommandQueue(int command_queue) {
	return command_queue;
}
int clReleaseProgram(int program) {
	return program;
}
int clRetainKernel(int kernel) {
	return kernel;
}
int clRetainProgram(int program) {
	return program;
}
cl_mem clCreateBuffer(int context, int flags, size_t size, void *value,
		      int *size_ret) {
	return 0;
}
int clEnqueueWriteBuffer(int command_queue, void *buffer, int blocking_write,
			 size_t offset, size_t cb, const void *ptr,
			 uint nevents, const int *wait_list, int *event) {
	return 0;
}
int clEnqueueReadBuffer(int command_queue, void *buffer, int blocking_write,
			size_t offset, size_t cb, const void *ptr, uint nevents,
			const int *wait_list, int *event) {
	return 0;
}
int clCreateKernel(int program, const char *kernel_name, int *errcode_ret) {
	return 0;
}
int clReleaseKernel(int kernel) {
	return kernel;
}
int clReleaseMemObject(void *memobj) {
	return 0;
}
int clSetKernelArg(int kernel, uint arg_index, size_t arg_size,
		   const void *arg_value) {
	return 0;
}
int clGetKernelWorkGroupInfo(int kernel, int device, int param_name,
			     size_t size, void *value, size_t * size_ret) {
	return 0;
}
int clEnqueueNDRangeKernel(int command_queue, int kernel, uint work_dim,
			   const size_t * offset, const size_t * size,
			   const size_t * local_size, uint nevents,
			   const int *wait_list, int *event) {
	return 0;
}
int clFinish(int command_queue) {
	return 0;
}
#endif
#define PRIME_NUM_CHECKS 20
#include <math.h>
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
static int reportLevel = 0;
const char *HashFactory_source;
const char *Hash_GetKernelSourceString() {
	return HashFactory_source;
}
size_t roundUpToNearest(size_t x, size_t r) {
	return (((x - 1) / r) + 1) * r;
}
int modularPow(int base, int exponent, int modulus) {
	int result = 1;
	while (exponent) {
		if (exponent & 1)
			result = ((long long int)result * base) % modulus;
		exponent >>= 1;
		base = ((long long int)base * base) % modulus;
	}
	return result;
}
int largestProthPrimeUnder(int N) {
	if (N < 4) {
		return N;
	}
	//determine the nearest proth number
	int n;
	int m;
	frexp((double)N, &n);
	n /= 2;
	int s = 1 << n;
	int p = s * ((N - 1) / s) + 1;
	int i;
	int a;
	srand(p);
	while (p > 3) {
		//check if a proth number is prime
		for (i = 0; i < PRIME_NUM_CHECKS; i++) {
			a = rand();
			if (modularPow(a, (p - 1) / 2, p) == p - 1) {
				return p;
			}
		}
		//determine the next proth number
		if (p - 1 == s * s / 4) {
			s /= 2;
		}
		p -= s;
	}
	return 3;
}
int smallestProthPrimeAbove(int N) {
	if (N < 4) {
		return N;
	}
	//determine the nearest proth number
	int n;
	int m;
	frexp((double)N, &n);
	n /= 2;
	int s = 1 << n;
	int p = s * ((N - 1) / s) + 1;
	int i;
	int a;
	srand(p);
	while (1) {
		//determine the next proth number
		if (p - 1 == s * s) {
			s *= 2;
		}
		p += s;
		//check if a proth number is prime
		for (i = 0; i < PRIME_NUM_CHECKS; i++) {
			a = rand();
			if (modularPow(a, (p - 1) / 2, p) == p - 1) {
				return p;
			}
		}
	}
	return 3;
}
int intLog2(int n) {
	int result = 0;
	while (n >>= 1) {
		result++;
	}
	return result;
}
void Hash_SetReportLevel(int level) {
	reportLevel = level;
}
int Hash_GetReportLevel() {
	return reportLevel;
}
char *Hash_ExitCodeString(int exitCode) {
	switch (exitCode) {
	case HASH_EXIT_CODE_NORMAL:
		return "Normal";
	case HASH_EXIT_CODE_ERROR:
		return "Error";
	case HASH_EXIT_CODE_OVERWRITE:
		return "Overwrite";
	case HASH_EXIT_CODE_KEY_DNE:
		return "Key Does Not Exist";
	case HASH_EXIT_CODE_CYCLE:
		return "Cycle";
	case HASH_EXIT_CODE_MAX_ENTRIES_EXCEEDED:
		return "Maximum Number Of Entries Exceeded";
	default:
		return "Unknown";
	}
}
void Hash_ExitCodeDebug(int exitCode) {
	if (exitCode != HASH_EXIT_CODE_NORMAL) {
		printf("HashExitCode: %s\n", Hash_ExitCodeString(exitCode));
	}
}

struct intintHash_Table_ {
	char *tableData;
	cl_mem tableDataBuffer;
	int (*destroyFunc) (intintHash_Table *);
	int (*setupFunc) (intintHash_Table *);
	int (*emptyFunc) (intintHash_Table *);
	int (*queryFunc) (intintHash_Table *, size_t, int *, int *);
	int (*querySingleFunc) (intintHash_Table *, int, int *);
	int (*insertFunc) (intintHash_Table *, size_t, int *, int *);
	int (*insertSingleFunc) (intintHash_Table *, int, int);
	int (*insertNoOverwriteFunc) (intintHash_Table *, size_t, int *, int *);
	int (*insertSingleNoOverwriteFunc) (intintHash_Table *, int, int);
	int (*bufferQueryFunc) (intintHash_Table *, size_t, cl_mem, cl_mem);
	int (*bufferInsertFunc) (intintHash_Table *, size_t, cl_mem, cl_mem);
	int (*bufferInsertNoOverwriteFunc) (intintHash_Table *, size_t, cl_mem,
					    cl_mem);
	cl_context context;
	cl_command_queue queue;
	cl_program utilProgram;
	cl_kernel emptyKernel;
	size_t emptyKernelLocalWorkSize;
	cl_program program;
	cl_kernel querySingleKernel;
	cl_kernel insertSingleKernel;
	cl_kernel insertSingleNoOverwriteKernel;
	size_t localWorkSize;
};
struct intintHash_Factory_ {
	cl_context context;
	cl_program program;
	cl_command_queue queue;
	int hashTypesAvailable;
	cl_program utilProgram[HASH_NUM_CL_HASHES];
	cl_kernel emptyKernel[HASH_NUM_CL_HASHES];
	size_t emptyKernelLocalWorkSize[HASH_NUM_CL_HASHES];
	cl_kernel querySingleKernel[HASH_NUM_CL_HASHES];
	cl_kernel insertSingleKernel[HASH_NUM_CL_HASHES];
	cl_kernel insertSingleNoOverwriteKernel[HASH_NUM_CL_HASHES];
	int emptyValue;
	size_t localWorkSize;
	intintHash_Table *(*createFunc[HASH_NUM_HASHES]) (intintHash_Factory *,
							  int hashIndex,
							  size_t keyRange,
							  size_t numEntries,
							  float loadFactor);
	int (*destroyFunc[HASH_NUM_HASHES]) (intintHash_Factory *,
					     int hashIndex);
};
intintHash_Factory *intintHash_CreateFactory(int hashTypes, int *emptyValue,
					     size_t localWorkSize,
					     cl_context * context,
					     cl_command_queue * queue) {
	if (hashTypes == 0) {
		hashTypes = HASH_ALL_C_HASHES;
	}
	if (!(hashTypes & HASH_ALL_HASHES)) {
		printf("Please specify a valid hash type to create.\n");
		exit(1);
	}
	hashTypes &= HASH_ALL_HASHES;
	if ((hashTypes & HASH_SENTINEL_PERFECT_HASHES) == hashTypes
	    && emptyValue == NULL) {
		printf
		    ("emptyValue must be valid if a sentinel perfect hash is the only option available.\n");
		exit(1);
	}
	intintHash_Factory *factory =
	    (intintHash_Factory *) malloc(sizeof(intintHash_Factory));
	if (emptyValue == NULL) {
		hashTypes &= !HASH_SENTINEL_PERFECT_HASHES;
	} else {
		factory->emptyValue = *emptyValue;
	}
	factory->hashTypesAvailable = hashTypes;
	if (hashTypes & HASH_ALL_CL_HASHES) {
		if (localWorkSize == 0) {
			factory->localWorkSize = HASH_DEFAULT_LOCAL_WORK_SIZE;
		} else {
			factory->localWorkSize = 1 << intLog2(localWorkSize);
		}
		if (context == NULL) {
			CLHash_Utilities_CreateContext(&factory->context,
						       &factory->queue);
		} else {
			factory->context = *context;
			clRetainContext(*context);
			if (queue == NULL) {
				printf
				    ("Please specify a command queue for your context.\n");
				exit(-1);
			}
			factory->queue = *queue;
			clRetainCommandQueue(*queue);
		}
		cl_int error;
		cl_device_id device;
		error =
		    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
				     sizeof(device), &device, NULL);
		CLHash_Utilities_HandleError(error, "intintHash_CreateFactory",
					     "clGetContextInfo");
		factory->program =
		    CLHash_Utilities_BuildProgramString(factory->context,
							device,
							Hash_GetKernelSourceString
							());
	}
	int hashType = 1;
	for (int hashIndex = 0; hashIndex < HASH_NUM_HASHES; hashIndex++) {
		hashType = 1 << hashIndex;
		switch (hashType & hashTypes) {
		case IDENTITY_PERFECT_HASH_ID:
			intintIdentityPerfectHash_CreateFactory(factory,
								hashIndex);
			break;
		case IDENTITY_PERFECT_CL_HASH_ID:
			intintIdentityPerfectCLHash_CreateFactory(factory,
								  hashIndex);
			break;
		case IDENTITY_PERFECT_OPENMP_HASH_ID:
			intintIdentityPerfectOpenMPHash_CreateFactory(factory,
								      hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_HASH_ID:
			intintIdentitySentinelPerfectHash_CreateFactory(factory,
									hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
			intintIdentitySentinelPerfectCLHash_CreateFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID:
			intintIdentitySentinelPerfectOpenMPHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_HASH_ID:
			intintLCGLinearOpenCompactHash_CreateFactory(factory,
								     hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
			intintLCGLinearOpenCompactCLHash_CreateFactory(factory,
								       hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID:
			intintLCGLinearOpenCompactOpenMPHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_HASH_ID:
			intintLCGQuadraticOpenCompactHash_CreateFactory(factory,
									hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
			intintLCGQuadraticOpenCompactCLHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID:
			intintLCGQuadraticOpenCompactOpenMPHash_CreateFactory
			    (factory, hashIndex);
			break;
		}
	}
	return factory;
}
int intintHash_DestroyFactory(intintHash_Factory * factory) {
	int hashType = 1;
	for (int hashIndex = 0; hashIndex < HASH_NUM_HASHES; hashIndex++) {
		hashType = 1 << hashIndex;
		switch (hashType & factory->hashTypesAvailable) {
		case IDENTITY_PERFECT_HASH_ID:
			intintIdentityPerfectHash_DestroyFactory(factory,
								 hashIndex);
			break;
		case IDENTITY_PERFECT_CL_HASH_ID:
			intintIdentityPerfectCLHash_DestroyFactory(factory,
								   hashIndex);
			break;
		case IDENTITY_PERFECT_OPENMP_HASH_ID:
			intintIdentityPerfectOpenMPHash_DestroyFactory(factory,
								       hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_HASH_ID:
			intintIdentitySentinelPerfectHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
			intintIdentitySentinelPerfectCLHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID:
			intintIdentitySentinelPerfectOpenMPHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_HASH_ID:
			intintLCGLinearOpenCompactHash_DestroyFactory(factory,
								      hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
			intintLCGLinearOpenCompactCLHash_DestroyFactory(factory,
									hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID:
			intintLCGLinearOpenCompactOpenMPHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_HASH_ID:
			intintLCGQuadraticOpenCompactHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
			intintLCGQuadraticOpenCompactCLHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID:
			intintLCGQuadraticOpenCompactOpenMPHash_DestroyFactory
			    (factory, hashIndex);
			break;
		}
		hashIndex++;
	}
	if (factory->hashTypesAvailable & HASH_ALL_CL_HASHES) {
		clReleaseContext(factory->context);
		clReleaseCommandQueue(factory->queue);
		clReleaseProgram(factory->program);
	}
	free(factory);
	return (0);
}
intintHash_Table *intintHash_CreateTable(intintHash_Factory * factory,
					 int hashTypes, size_t keyRange,
					 size_t numEntries, float loadFactor) {
	if (loadFactor > 1.0 || loadFactor < HASH_MIN_LOAD_FACTOR) {
		loadFactor = HASH_DEFAULT_LOAD_FACTOR;
	}
	if (hashTypes == 0) {
		hashTypes = factory->hashTypesAvailable;
		if ((hashTypes & HASH_ALL_CL_HASHES)
		    && (hashTypes & HASH_ALL_OPENMP_HASHES)
		    && (hashTypes & HASH_ALL_C_HASHES)) {
			hashTypes &= HASH_ALL_C_HASHES;
		}
	}
	if (!(hashTypes & factory->hashTypesAvailable)) {
		printf
		    ("None of the selected hash types are supported by this factory.\n");
		exit(1);
	}
	hashTypes &= factory->hashTypesAvailable;
	if ((hashTypes & HASH_ALL_CL_HASHES)
	    && (hashTypes & HASH_ALL_OPENMP_HASHES)
	    && (hashTypes & HASH_ALL_C_HASHES)) {
		printf("Please decide between OpenCL, OpenMP or C hash.\n");
		exit(1);
	}
	if ((hashTypes & HASH_PERFECT_HASHES) == hashTypes && keyRange == 0) {
		printf
		    ("keyRange must be set if a perfect hash is the only option available.\n");
		exit(1);
	}
	if ((hashTypes & HASH_COMPACT_HASHES) == hashTypes && numEntries == 0) {
		printf
		    ("numEntries must be set if a compact hash is the only option available.\n");
		exit(1);
	}
	if (numEntries == 0 && keyRange == 0) {
		printf("either numEntries or keyRange must be set.\n");
		exit(1);
	}
	size_t perfectNumBuckets = keyRange;
	size_t compactNumBuckets = (size_t) (numEntries / loadFactor);
	int hashIndex;
	if ((hashTypes & HASH_SENTINEL_PERFECT_HASHES)
	    && ((hashTypes == (hashTypes & HASH_SENTINEL_PERFECT_HASHES))
		|| (compactNumBuckets == 0
		    || (perfectNumBuckets / compactNumBuckets <
			HASH_PERFECT_COMPACT_SWITCH_FACTOR)))) {
		hashIndex = intLog2(hashTypes & HASH_SENTINEL_PERFECT_HASHES);
	} else if ((hashTypes & HASH_NOSENTINEL_PERFECT_HASHES)
		   && ((hashTypes == (hashTypes & HASH_PERFECT_HASHES))
		       || (compactNumBuckets == 0
			   || (perfectNumBuckets / compactNumBuckets <
			       HASH_PERFECT_COMPACT_SWITCH_FACTOR)))) {
		hashIndex = intLog2(hashTypes & HASH_NOSENTINEL_PERFECT_HASHES);
	} else if ((hashTypes & HASH_LINEAR_COMPACT_HASHES)
		   &&
		   ((hashTypes ==
		     (hashTypes &
		      (HASH_PERFECT_HASHES | HASH_LINEAR_COMPACT_HASHES)))
		    ||
		    ((hashTypes & HASH_COMPACT_HASHES &
		      HASH_LINEAR_COMPACT_HASHES) ==
		     (hashTypes & HASH_COMPACT_HASHES) || loadFactor > 0.5))) {
		hashIndex = intLog2(hashTypes & HASH_LINEAR_COMPACT_HASHES);
	} else {
		hashIndex = intLog2(hashTypes & HASH_QUADRATIC_COMPACT_HASHES);
	}
	intintHash_Table *table =
	    factory->createFunc[hashIndex] (factory, hashIndex, keyRange,
					    numEntries, loadFactor);
	return table;
}
int intintHash_SetupTable(intintHash_Table * table) {
	table->setupFunc(table);
	return (0);
}
int intintHash_EmptyTable(intintHash_Table * table) {
	table->emptyFunc(table);
	return (0);
}
int intintHash_DestroyTable(intintHash_Table * table) {
	table->destroyFunc(table);
	return (0);
}
cl_mem intintHash_GetTableDataBuffer(intintHash_Table * table) {
	return table->tableDataBuffer;
}
cl_mem *intintHash_GetTableDataBufferPtr(intintHash_Table * table) {
	return &table->tableDataBuffer;
}
int intintHash_GetTableType(intintHash_Table * table) {
	return ((int *)table->tableData)[0];
} int intintHash_Query(intintHash_Table * table, size_t numKeys, int *keys,
		       int *valuesOutput) {
	table->queryFunc(table, numKeys, keys, valuesOutput);
	return (0);
}
int intintHash_QuerySingle(intintHash_Table * table, int key, int *valueOutput) {
	table->querySingleFunc(table, key, valueOutput);
	return (0);
}
int intintHash_Insert(intintHash_Table * table, size_t numEntries, int *keys,
		      int *values) {
	table->insertFunc(table, numEntries, keys, values);
	return (0);
}
int intintHash_InsertSingle(intintHash_Table * table, int key, int value) {
	table->insertSingleFunc(table, key, value);
	return (0);
}
int intintHash_InsertNoOverwrite(intintHash_Table * table, size_t numEntries,
				 int *keys, int *values) {
	table->insertNoOverwriteFunc(table, numEntries, keys, values);
	return (0);
}
int intintHash_InsertSingleNoOverwrite(intintHash_Table * table, int key,
				       int value) {
	table->insertSingleNoOverwriteFunc(table, key, value);
	return (0);
}
int intintHash_BufferQuery(intintHash_Table * table, size_t numKeys,
			   cl_mem keys, cl_mem valuesOutput) {
	table->bufferQueryFunc(table, numKeys, keys, valuesOutput);
	return (0);
}
int intintHash_BufferInsert(intintHash_Table * table, size_t numEntries,
			    cl_mem keys, cl_mem values) {
	table->bufferInsertFunc(table, numEntries, keys, values);
	return (0);
}
int intintHash_BufferInsertNoOverwrite(intintHash_Table * table,
				       size_t numEntries, cl_mem keys,
				       cl_mem values) {
	table->bufferInsertNoOverwriteFunc(table, numEntries, keys, values);
	return (0);
}

typedef struct intintIdentityPerfectHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
} intintIdentityPerfectHash_TableData;
typedef struct intintIdentityPerfectHash_Bucket {
	int key;
	int value;
} intintIdentityPerfectHash_Bucket;
intintHash_Table *intintIdentityPerfectHash_CreateTable(intintHash_Factory *
							factory, int hashIndex,
							size_t keyRange,
							size_t numEntries,
							float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintIdentityPerfectHash_DestroyTable;
	table->setupFunc = &intintIdentityPerfectHash_SetupTable;
	table->emptyFunc = &intintIdentityPerfectHash_EmptyTable;
	table->queryFunc = &intintIdentityPerfectHash_Query;
	table->querySingleFunc = &intintIdentityPerfectHash_QuerySingle;
	table->insertFunc = &intintIdentityPerfectHash_Insert;
	table->insertSingleFunc = &intintIdentityPerfectHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintIdentityPerfectHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintIdentityPerfectHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(intintIdentityPerfectHash_TableData));
	((intintIdentityPerfectHash_TableData *) table->tableData)->hashID =
	    IDENTITY_PERFECT_HASH_ID;
	((intintIdentityPerfectHash_TableData *) table->tableData)->numBuckets =
	    keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(intintIdentityPerfectHash_TableData) +
			   ((intintIdentityPerfectHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(intintIdentityPerfectHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintIdentityPerfectHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintIdentityPerfectHash_CreateFactory(intintHash_Factory * factory,
					    int hashIndex) {
	factory->createFunc[hashIndex] = &intintIdentityPerfectHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintIdentityPerfectHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectHash_DestroyFactory(intintHash_Factory * factory,
					     int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintIdentityPerfectHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) & table->
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (int index = 0;
		     index <
		     ((intintIdentityPerfectHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentityPerfectHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) & table->
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	for (int index = 0;
	     index <
	     ((intintIdentityPerfectHash_TableData *) table->tableData)->
	     numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentityPerfectHash_InnerQuerySingle(char *tableData, int key,
					       int *valueOutput) {
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentityPerfectHash_TableData *)
					 tableData)->compressFuncData, key);
	if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintIdentityPerfectHash_InnerQuery(char *tableData, unsigned int numKeys,
					 int *keys, int *valuesOutput) {
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    intintHash_CompressIdentity(((intintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
		if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentityPerfectHash_InnerInsertSingle(char *tableData, int key,
						int value) {
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentityPerfectHash_TableData *)
					 tableData)->compressFuncData, key);
	if (((buckets[index].key ==
	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
					   key,
					   HASH_BUCKET_STATUS_EMPTY) :
	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentityPerfectHash_InnerInsert(char *tableData,
					  unsigned int numEntries, int *keys,
					  int *values) {
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentityPerfectHash_InnerInsertSingleNoOverwrite(char *tableData,
							   int key, int value) {
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentityPerfectHash_TableData *)
					 tableData)->compressFuncData, key);
	if (((buckets[index].key ==
	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
					   key,
					   HASH_BUCKET_STATUS_EMPTY) :
	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentityPerfectHash_InnerInsertNoOverwrite(char *tableData,
						     unsigned int numEntries,
						     int *keys, int *values) {
	intintIdentityPerfectHash_Bucket *buckets =
	    (intintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentityPerfectHash_QuerySingle(intintHash_Table * table, int key,
					  int *valueOutput) {
	return intintIdentityPerfectHash_InnerQuerySingle(table->tableData, key,
							  valueOutput);
}
int intintIdentityPerfectHash_Query(intintHash_Table * table, size_t numKeys,
				    int *keys, int *valuesOutput) {
	return intintIdentityPerfectHash_InnerQuery(table->tableData, numKeys,
						    keys, valuesOutput);
}
int intintIdentityPerfectHash_InsertSingle(intintHash_Table * table, int key,
					   int value) {
	return intintIdentityPerfectHash_InnerInsertSingle(table->tableData,
							   key, value);
}
int intintIdentityPerfectHash_Insert(intintHash_Table * table,
				     size_t numEntries, int *keys,
				     int *values) {
	return intintIdentityPerfectHash_InnerInsert(table->tableData,
						     numEntries, keys, values);
}
int intintIdentityPerfectHash_InsertSingleNoOverwrite(intintHash_Table * table,
						      int key, int value) {
	return intintIdentityPerfectHash_InnerInsertSingleNoOverwrite(table->
								      tableData,
								      key,
								      value);
}
int intintIdentityPerfectHash_InsertNoOverwrite(intintHash_Table * table,
						size_t numEntries, int *keys,
						int *values) {
	return intintIdentityPerfectHash_InnerInsertNoOverwrite(table->
								tableData,
								numEntries,
								keys, values);
}

typedef struct intintIdentityPerfectCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
} intintIdentityPerfectCLHash_TableData;
typedef struct intintIdentityPerfectCLHash_Bucket {
	int key;
	int value;
} intintIdentityPerfectCLHash_Bucket;
intintHash_Table *intintIdentityPerfectCLHash_CreateTable(intintHash_Factory *
							  factory,
							  int hashIndex,
							  size_t keyRange,
							  size_t numEntries,
							  float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintIdentityPerfectCLHash_DestroyTable;
	table->setupFunc = &intintIdentityPerfectCLHash_SetupTable;
	table->emptyFunc = &intintIdentityPerfectCLHash_EmptyTable;
	table->queryFunc = &intintIdentityPerfectCLHash_Query;
	table->querySingleFunc = &intintIdentityPerfectCLHash_QuerySingle;
	table->insertFunc = &intintIdentityPerfectCLHash_Insert;
	table->insertSingleFunc = &intintIdentityPerfectCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintIdentityPerfectCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintIdentityPerfectCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(intintIdentityPerfectCLHash_TableData));
	((intintIdentityPerfectCLHash_TableData *) table->tableData)->hashID =
	    IDENTITY_PERFECT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((intintIdentityPerfectCLHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(intintIdentityPerfectCLHash_TableData) +
			   ((intintIdentityPerfectCLHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(intintIdentityPerfectCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintIdentityPerfectCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(intintIdentityPerfectHash_TableData) +
			   ((intintIdentityPerfectHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(intintIdentityPerfectHash_Bucket), NULL,
			   &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0, sizeof(intintIdentityPerfectHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int intintIdentityPerfectCLHash_CreateFactory(intintHash_Factory * factory,
					      int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintIdentityPerfectCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintIdentityPerfectCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "intintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintIdentityPerfectCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintIdentityPerfectCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int intintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct intintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }intintHash_CompressLCGData; static inline unsigned int intintHash_CompressLCG(intintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct intintIdentityPerfectCLHash_TableData{ int hashID; unsigned int numBuckets; char compressFuncData; }intintIdentityPerfectCLHash_TableData; typedef struct intintIdentityPerfectCLHash_Bucket{ int key; int value; }intintIdentityPerfectCLHash_Bucket; __kernel void intintIdentityPerfectCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global intintIdentityPerfectCLHash_TableData*)tableData)->numBuckets){ return; } __global intintIdentityPerfectCLHash_Bucket *buckets = (__global intintIdentityPerfectCLHash_Bucket*)&tableData[sizeof(intintIdentityPerfectCLHash_TableData)]; buckets[index].key = -1;/*HASH_BUCKET_STATUS_EMPTY*/ }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "intintIdentityPerfectCLHash_Empty", &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "intintIdentityPerfectCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectCLHash_DestroyFactory(intintHash_Factory * factory,
					       int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectCLHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintIdentityPerfectCLHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintIdentityPerfectHash_TableData *) table->
			      tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintIdentityPerfectCLHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintIdentityPerfectHash_TableData *) table->
			      tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintIdentityPerfectCLHash_QuerySingle(intintHash_Table * table, int key,
					    int *valueOutput) {
	return intintIdentityPerfectCLHash_Query(table, 1, &key, valueOutput);
}
int intintIdentityPerfectCLHash_Query(intintHash_Table * table, size_t numKeys,
				      int *keys, int *valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err, "intintIdentityPerfectCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(int) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err, "intintIdentityPerfectCLHash_Query",
				     "clCreateBuffer");
	intintIdentityPerfectCLHash_BufferQuery(table, numKeys, keysBuffer,
						valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(int) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err, "intintIdentityPerfectCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectCLHash_BufferQuery(intintHash_Table * table,
					    size_t numKeys, cl_mem keysBuffer,
					    cl_mem valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectCLHash_InsertSingle(intintHash_Table * table, int key,
					     int value) {
	return intintIdentityPerfectCLHash_Insert(table, 1, &key, &value);
}
int intintIdentityPerfectCLHash_Insert(intintHash_Table * table,
				       size_t numEntries, int *keys,
				       int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err, "intintIdentityPerfectCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err, "intintIdentityPerfectCLHash_Insert",
				     "clCreateBuffer");
	intintIdentityPerfectCLHash_BufferInsert(table, numEntries, keysBuffer,
						 valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectCLHash_BufferInsert(intintHash_Table * table,
					     size_t numEntries,
					     cl_mem keysBuffer,
					     cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int intintIdentityPerfectCLHash_InsertSingleNoOverwrite(intintHash_Table *
							table, int key,
							int value) {
	return intintIdentityPerfectCLHash_InsertNoOverwrite(table, 1, &key,
							     &value);
}
int intintIdentityPerfectCLHash_InsertNoOverwrite(intintHash_Table * table,
						  size_t numEntries, int *keys,
						  int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	intintIdentityPerfectCLHash_BufferInsertNoOverwrite(table, numEntries,
							    keysBuffer,
							    valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectCLHash_BufferInsertNoOverwrite(intintHash_Table *
							table,
							size_t numEntries,
							cl_mem keysBuffer,
							cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct intintIdentityPerfectOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
} intintIdentityPerfectOpenMPHash_TableData;
typedef struct intintIdentityPerfectOpenMPHash_Bucket {
	int key;
	int value;
} intintIdentityPerfectOpenMPHash_Bucket;
intintHash_Table *intintIdentityPerfectOpenMPHash_CreateTable(intintHash_Factory
							      * factory,
							      int hashIndex,
							      size_t keyRange,
							      size_t numEntries,
							      float loadFactor) 
{
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintIdentityPerfectOpenMPHash_DestroyTable;
	table->setupFunc = &intintIdentityPerfectOpenMPHash_SetupTable;
	table->emptyFunc = &intintIdentityPerfectOpenMPHash_EmptyTable;
	table->queryFunc = &intintIdentityPerfectOpenMPHash_Query;
	table->querySingleFunc = &intintIdentityPerfectOpenMPHash_QuerySingle;
	table->insertFunc = &intintIdentityPerfectOpenMPHash_Insert;
	table->insertSingleFunc = &intintIdentityPerfectOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintIdentityPerfectOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintIdentityPerfectOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(intintIdentityPerfectOpenMPHash_TableData));
	((intintIdentityPerfectOpenMPHash_TableData *) table->tableData)->
	    hashID = IDENTITY_PERFECT_OPENMP_HASH_ID;
	((intintIdentityPerfectOpenMPHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(intintIdentityPerfectOpenMPHash_TableData) +
			   ((intintIdentityPerfectOpenMPHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintIdentityPerfectOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintIdentityPerfectOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintIdentityPerfectOpenMPHash_CreateFactory(intintHash_Factory * factory,
						  int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintIdentityPerfectOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintIdentityPerfectOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectOpenMPHash_DestroyFactory(intintHash_Factory * factory,
						   int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentityPerfectOpenMPHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintIdentityPerfectOpenMPHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (int index = 0;
		     index <
		     ((intintIdentityPerfectOpenMPHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentityPerfectOpenMPHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
#pragma omp parallel for
	for (int index = 0;
	     index <
	     ((intintIdentityPerfectOpenMPHash_TableData *) table->tableData)->
	     numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentityPerfectOpenMPHash_InnerQuerySingle(char *tableData, int key,
						     int *valueOutput) {
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintIdentityPerfectOpenMPHash_InnerQuery(char *tableData,
					       unsigned int numKeys, int *keys,
					       int *valuesOutput) {
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    intintHash_CompressIdentity(((intintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentityPerfectOpenMPHash_InnerInsertSingle(char *tableData, int key,
						      int value) {
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (((buckets[index].key ==
	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
					   key,
					   HASH_BUCKET_STATUS_EMPTY) :
	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentityPerfectOpenMPHash_InnerInsert(char *tableData,
						unsigned int numEntries,
						int *keys, int *values) {
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentityPerfectOpenMPHash_InnerInsertSingleNoOverwrite(char
								 *tableData,
								 int key,
								 int value) {
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (((buckets[index].key ==
	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
					   key,
					   HASH_BUCKET_STATUS_EMPTY) :
	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentityPerfectOpenMPHash_InnerInsertNoOverwrite(char *tableData,
							   unsigned int
							   numEntries,
							   int *keys,
							   int *values) {
	intintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentityPerfectOpenMPHash_QuerySingle(intintHash_Table * table,
						int key, int *valueOutput) {
	return intintIdentityPerfectOpenMPHash_InnerQuerySingle(table->
								tableData, key,
								valueOutput);
}
int intintIdentityPerfectOpenMPHash_Query(intintHash_Table * table,
					  size_t numKeys, int *keys,
					  int *valuesOutput) {
	return intintIdentityPerfectOpenMPHash_InnerQuery(table->tableData,
							  numKeys, keys,
							  valuesOutput);
}
int intintIdentityPerfectOpenMPHash_InsertSingle(intintHash_Table * table,
						 int key, int value) {
	return intintIdentityPerfectOpenMPHash_InnerInsertSingle(table->
								 tableData, key,
								 value);
}
int intintIdentityPerfectOpenMPHash_Insert(intintHash_Table * table,
					   size_t numEntries, int *keys,
					   int *values) {
	return intintIdentityPerfectOpenMPHash_InnerInsert(table->tableData,
							   numEntries, keys,
							   values);
}
int intintIdentityPerfectOpenMPHash_InsertSingleNoOverwrite(intintHash_Table *
							    table, int key,
							    int value) {
	return
	    intintIdentityPerfectOpenMPHash_InnerInsertSingleNoOverwrite(table->
									 tableData,
									 key,
									 value);
}
int intintIdentityPerfectOpenMPHash_InsertNoOverwrite(intintHash_Table * table,
						      size_t numEntries,
						      int *keys, int *values) {
	return intintIdentityPerfectOpenMPHash_InnerInsertNoOverwrite(table->
								      tableData,
								      numEntries,
								      keys,
								      values);
}

typedef struct intintIdentitySentinelPerfectHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
	int emptyValue;
} intintIdentitySentinelPerfectHash_TableData;
typedef struct intintIdentitySentinelPerfectHash_Bucket {
	int value;
} intintIdentitySentinelPerfectHash_Bucket;
intintHash_Table
    *intintIdentitySentinelPerfectHash_CreateTable(intintHash_Factory * factory,
						   int hashIndex,
						   size_t keyRange,
						   size_t numEntries,
						   float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintIdentitySentinelPerfectHash_DestroyTable;
	table->setupFunc = &intintIdentitySentinelPerfectHash_SetupTable;
	table->emptyFunc = &intintIdentitySentinelPerfectHash_EmptyTable;
	table->queryFunc = &intintIdentitySentinelPerfectHash_Query;
	table->querySingleFunc = &intintIdentitySentinelPerfectHash_QuerySingle;
	table->insertFunc = &intintIdentitySentinelPerfectHash_Insert;
	table->insertSingleFunc =
	    &intintIdentitySentinelPerfectHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintIdentitySentinelPerfectHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintIdentitySentinelPerfectHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(intintIdentitySentinelPerfectHash_TableData));
	((intintIdentitySentinelPerfectHash_TableData *) table->tableData)->
	    hashID = IDENTITY_SENTINEL_PERFECT_HASH_ID;
	((intintIdentitySentinelPerfectHash_TableData *) table->tableData)->
	    emptyValue = factory->emptyValue;
	((intintIdentitySentinelPerfectHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(intintIdentitySentinelPerfectHash_TableData) +
			   ((intintIdentitySentinelPerfectHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintIdentitySentinelPerfectHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintIdentitySentinelPerfectHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintIdentitySentinelPerfectHash_CreateFactory(intintHash_Factory *
						    factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintIdentitySentinelPerfectHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintIdentitySentinelPerfectHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectHash_DestroyFactory(intintHash_Factory *
						     factory, int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintIdentitySentinelPerfectHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) & table->
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (int index = 0;
		     index <
		     ((intintIdentitySentinelPerfectHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].value =
			    ((intintIdentitySentinelPerfectHash_TableData *)
			     table->tableData)->emptyValue;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentitySentinelPerfectHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) & table->
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	for (int index = 0;
	     index <
	     ((intintIdentitySentinelPerfectHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].value =
		    ((intintIdentitySentinelPerfectHash_TableData *) table->
		     tableData)->emptyValue;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentitySentinelPerfectHash_InnerQuerySingle(char *tableData, int key,
						       int *valueOutput) {
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((intintIdentitySentinelPerfectHash_TableData *) tableData)->
	    emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintIdentitySentinelPerfectHash_InnerQuery(char *tableData,
						 unsigned int numKeys,
						 int *keys, int *valuesOutput) {
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    intintHash_CompressIdentity(((intintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((intintIdentitySentinelPerfectHash_TableData *)
		     tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentitySentinelPerfectHash_InnerInsertSingle(char *tableData,
							int key, int value) {
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((intintIdentitySentinelPerfectHash_TableData *) tableData)->
	    emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentitySentinelPerfectHash_InnerInsert(char *tableData,
						  unsigned int numEntries,
						  int *keys, int *values) {
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((intintIdentitySentinelPerfectHash_TableData *)
		     tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentitySentinelPerfectHash_InnerInsertSingleNoOverwrite(char
								   *tableData,
								   int key,
								   int value) {
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((intintIdentitySentinelPerfectHash_TableData *) tableData)->
	    emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentitySentinelPerfectHash_InnerInsertNoOverwrite(char *tableData,
							     unsigned int
							     numEntries,
							     int *keys,
							     int *values) {
	intintIdentitySentinelPerfectHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((intintIdentitySentinelPerfectHash_TableData *)
		     tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentitySentinelPerfectHash_QuerySingle(intintHash_Table * table,
						  int key, int *valueOutput) {
	return intintIdentitySentinelPerfectHash_InnerQuerySingle(table->
								  tableData,
								  key,
								  valueOutput);
}
int intintIdentitySentinelPerfectHash_Query(intintHash_Table * table,
					    size_t numKeys, int *keys,
					    int *valuesOutput) {
	return intintIdentitySentinelPerfectHash_InnerQuery(table->tableData,
							    numKeys, keys,
							    valuesOutput);
}
int intintIdentitySentinelPerfectHash_InsertSingle(intintHash_Table * table,
						   int key, int value) {
	return intintIdentitySentinelPerfectHash_InnerInsertSingle(table->
								   tableData,
								   key, value);
}
int intintIdentitySentinelPerfectHash_Insert(intintHash_Table * table,
					     size_t numEntries, int *keys,
					     int *values) {
	return intintIdentitySentinelPerfectHash_InnerInsert(table->tableData,
							     numEntries, keys,
							     values);
}
int intintIdentitySentinelPerfectHash_InsertSingleNoOverwrite(intintHash_Table *
							      table, int key,
							      int value) {
	return
	    intintIdentitySentinelPerfectHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int intintIdentitySentinelPerfectHash_InsertNoOverwrite(intintHash_Table *
							table,
							size_t numEntries,
							int *keys,
							int *values) {
	return intintIdentitySentinelPerfectHash_InnerInsertNoOverwrite(table->
									tableData,
									numEntries,
									keys,
									values);
}

typedef struct intintIdentitySentinelPerfectCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
	int emptyValue;
} intintIdentitySentinelPerfectCLHash_TableData;
typedef struct intintIdentitySentinelPerfectCLHash_Bucket {
	int value;
} intintIdentitySentinelPerfectCLHash_Bucket;
intintHash_Table
    *intintIdentitySentinelPerfectCLHash_CreateTable(intintHash_Factory *
						     factory, int hashIndex,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintIdentitySentinelPerfectCLHash_DestroyTable;
	table->setupFunc = &intintIdentitySentinelPerfectCLHash_SetupTable;
	table->emptyFunc = &intintIdentitySentinelPerfectCLHash_EmptyTable;
	table->queryFunc = &intintIdentitySentinelPerfectCLHash_Query;
	table->querySingleFunc =
	    &intintIdentitySentinelPerfectCLHash_QuerySingle;
	table->insertFunc = &intintIdentitySentinelPerfectCLHash_Insert;
	table->insertSingleFunc =
	    &intintIdentitySentinelPerfectCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintIdentitySentinelPerfectCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(intintIdentitySentinelPerfectCLHash_TableData));
	((intintIdentitySentinelPerfectCLHash_TableData *) table->tableData)->
	    hashID = IDENTITY_SENTINEL_PERFECT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((intintIdentitySentinelPerfectCLHash_TableData *) table->tableData)->
	    emptyValue = factory->emptyValue;
	((intintIdentitySentinelPerfectCLHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(intintIdentitySentinelPerfectCLHash_TableData)
			   +
			   ((intintIdentitySentinelPerfectCLHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintIdentitySentinelPerfectCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintIdentitySentinelPerfectCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(intintIdentitySentinelPerfectHash_TableData) +
			   ((intintIdentitySentinelPerfectHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintIdentitySentinelPerfectHash_Bucket),
			   NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0,
				 sizeof
				 (intintIdentitySentinelPerfectHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int intintIdentitySentinelPerfectCLHash_CreateFactory(intintHash_Factory *
						      factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintIdentitySentinelPerfectCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintIdentitySentinelPerfectCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "intintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintIdentitySentinelPerfectCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintIdentitySentinelPerfectCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int intintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct intintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }intintHash_CompressLCGData; static inline unsigned int intintHash_CompressLCG(intintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct intintIdentitySentinelPerfectCLHash_TableData{ int hashID; unsigned int numBuckets; char compressFuncData; int emptyValue; }intintIdentitySentinelPerfectCLHash_TableData; typedef struct intintIdentitySentinelPerfectCLHash_Bucket{ int value; }intintIdentitySentinelPerfectCLHash_Bucket; __kernel void intintIdentitySentinelPerfectCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global intintIdentitySentinelPerfectCLHash_TableData*)tableData)->numBuckets){ return; } __global intintIdentitySentinelPerfectCLHash_Bucket *buckets = (__global intintIdentitySentinelPerfectCLHash_Bucket*)&tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)]; buckets[index].value = ((__global intintIdentitySentinelPerfectCLHash_TableData*)tableData)->emptyValue; }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "intintIdentitySentinelPerfectCLHash_Empty", &error);
	CLHash_Utilities_HandleError(error,
				     "intintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "intintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectCLHash_DestroyFactory(intintHash_Factory *
						       factory,
						       int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectCLHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintIdentitySentinelPerfectCLHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintIdentitySentinelPerfectHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintIdentitySentinelPerfectCLHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintIdentitySentinelPerfectHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintIdentitySentinelPerfectCLHash_QuerySingle(intintHash_Table * table,
						    int key, int *valueOutput) {
	return intintIdentitySentinelPerfectCLHash_Query(table, 1, &key,
							 valueOutput);
}
int intintIdentitySentinelPerfectCLHash_Query(intintHash_Table * table,
					      size_t numKeys, int *keys,
					      int *valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(int) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_Query",
				     "clCreateBuffer");
	intintIdentitySentinelPerfectCLHash_BufferQuery(table, numKeys,
							keysBuffer,
							valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(int) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectCLHash_BufferQuery(intintHash_Table * table,
						    size_t numKeys,
						    cl_mem keysBuffer,
						    cl_mem valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectCLHash_InsertSingle(intintHash_Table * table,
						     int key, int value) {
	return intintIdentitySentinelPerfectCLHash_Insert(table, 1, &key,
							  &value);
}
int intintIdentitySentinelPerfectCLHash_Insert(intintHash_Table * table,
					       size_t numEntries, int *keys,
					       int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_Insert",
				     "clCreateBuffer");
	intintIdentitySentinelPerfectCLHash_BufferInsert(table, numEntries,
							 keysBuffer,
							 valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectCLHash_BufferInsert(intintHash_Table * table,
						     size_t numEntries,
						     cl_mem keysBuffer,
						     cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(intintHash_Table
								* table,
								int key,
								int value) {
	return intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(table, 1,
								     &key,
								     &value);
}
int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(intintHash_Table *
							  table,
							  size_t numEntries,
							  int *keys,
							  int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite(table,
								    numEntries,
								    keysBuffer,
								    valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite(intintHash_Table
								* table,
								size_t
								numEntries,
								cl_mem
								keysBuffer,
								cl_mem
								valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct intintIdentitySentinelPerfectOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
	int emptyValue;
} intintIdentitySentinelPerfectOpenMPHash_TableData;
typedef struct intintIdentitySentinelPerfectOpenMPHash_Bucket {
	int value;
} intintIdentitySentinelPerfectOpenMPHash_Bucket;
intintHash_Table
    *intintIdentitySentinelPerfectOpenMPHash_CreateTable(intintHash_Factory *
							 factory, int hashIndex,
							 size_t keyRange,
							 size_t numEntries,
							 float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc =
	    &intintIdentitySentinelPerfectOpenMPHash_DestroyTable;
	table->setupFunc = &intintIdentitySentinelPerfectOpenMPHash_SetupTable;
	table->emptyFunc = &intintIdentitySentinelPerfectOpenMPHash_EmptyTable;
	table->queryFunc = &intintIdentitySentinelPerfectOpenMPHash_Query;
	table->querySingleFunc =
	    &intintIdentitySentinelPerfectOpenMPHash_QuerySingle;
	table->insertFunc = &intintIdentitySentinelPerfectOpenMPHash_Insert;
	table->insertSingleFunc =
	    &intintIdentitySentinelPerfectOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintIdentitySentinelPerfectOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintIdentitySentinelPerfectOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(intintIdentitySentinelPerfectOpenMPHash_TableData));
	((intintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	 tableData)->hashID = IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID;
	((intintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	 tableData)->emptyValue = factory->emptyValue;
	((intintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	 tableData)->numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)
	    malloc(sizeof(intintIdentitySentinelPerfectOpenMPHash_TableData) +
		   ((intintIdentitySentinelPerfectOpenMPHash_TableData *)
		    table->tableData)->numBuckets *
		   sizeof(intintIdentitySentinelPerfectOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintIdentitySentinelPerfectOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintIdentitySentinelPerfectOpenMPHash_CreateFactory(intintHash_Factory *
							  factory,
							  int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintIdentitySentinelPerfectOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintIdentitySentinelPerfectOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectOpenMPHash_DestroyFactory(intintHash_Factory *
							   factory,
							   int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintIdentitySentinelPerfectOpenMPHash_DestroyTable(intintHash_Table *
							 table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintIdentitySentinelPerfectOpenMPHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (int index = 0;
		     index <
		     ((intintIdentitySentinelPerfectOpenMPHash_TableData *)
		      table->tableData)->numBuckets; index++) {
			buckets[index].value =
			    ((intintIdentitySentinelPerfectOpenMPHash_TableData
			      *) table->tableData)->emptyValue;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentitySentinelPerfectOpenMPHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
#pragma omp parallel for
	for (int index = 0;
	     index <
	     ((intintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].value =
		    ((intintIdentitySentinelPerfectOpenMPHash_TableData *)
		     table->tableData)->emptyValue;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintIdentitySentinelPerfectOpenMPHash_InnerQuerySingle(char *tableData,
							     int key,
							     int *valueOutput) {
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->
	    emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintIdentitySentinelPerfectOpenMPHash_InnerQuery(char *tableData,
						       unsigned int numKeys,
						       int *keys,
						       int *valuesOutput) {
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    intintHash_CompressIdentity(((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((intintIdentitySentinelPerfectOpenMPHash_TableData *)
		     tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentitySentinelPerfectOpenMPHash_InnerInsertSingle(char *tableData,
							      int key,
							      int value) {
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->
	    emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentitySentinelPerfectOpenMPHash_InnerInsert(char *tableData,
							unsigned int numEntries,
							int *keys,
							int *values) {
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((intintIdentitySentinelPerfectOpenMPHash_TableData *)
		     tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentitySentinelPerfectOpenMPHash_InnerInsertSingleNoOverwrite(char
									 *tableData,
									 int
									 key,
									 int
									 value) 
{
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->
	    emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintIdentitySentinelPerfectOpenMPHash_InnerInsertNoOverwrite(char
								   *tableData,
								   unsigned int
								   numEntries,
								   int *keys,
								   int *values) 
{
	intintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (intintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintIdentitySentinelPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((intintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((intintIdentitySentinelPerfectOpenMPHash_TableData *)
		     tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintIdentitySentinelPerfectOpenMPHash_QuerySingle(intintHash_Table *
							table, int key,
							int *valueOutput) {
	return intintIdentitySentinelPerfectOpenMPHash_InnerQuerySingle(table->
									tableData,
									key,
									valueOutput);
}
int intintIdentitySentinelPerfectOpenMPHash_Query(intintHash_Table * table,
						  size_t numKeys, int *keys,
						  int *valuesOutput) {
	return intintIdentitySentinelPerfectOpenMPHash_InnerQuery(table->
								  tableData,
								  numKeys, keys,
								  valuesOutput);
}
int intintIdentitySentinelPerfectOpenMPHash_InsertSingle(intintHash_Table *
							 table, int key,
							 int value) {
	return intintIdentitySentinelPerfectOpenMPHash_InnerInsertSingle(table->
									 tableData,
									 key,
									 value);
}
int intintIdentitySentinelPerfectOpenMPHash_Insert(intintHash_Table * table,
						   size_t numEntries, int *keys,
						   int *values) {
	return intintIdentitySentinelPerfectOpenMPHash_InnerInsert(table->
								   tableData,
								   numEntries,
								   keys,
								   values);
}
int
intintIdentitySentinelPerfectOpenMPHash_InsertSingleNoOverwrite(intintHash_Table
								* table,
								int key,
								int value) {
	return
	    intintIdentitySentinelPerfectOpenMPHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int intintIdentitySentinelPerfectOpenMPHash_InsertNoOverwrite(intintHash_Table *
							      table,
							      size_t numEntries,
							      int *keys,
							      int *values) {
	return
	    intintIdentitySentinelPerfectOpenMPHash_InnerInsertNoOverwrite
	    (table->tableData, numEntries, keys, values);
}

typedef struct intintLCGLinearOpenCompactHash_TableData {
	int hashID;
	unsigned int numBuckets;
	intintHash_CompressLCGData compressFuncData;
} intintLCGLinearOpenCompactHash_TableData;
typedef struct intintLCGLinearOpenCompactHash_Bucket {
	int key;
	int value;
} intintLCGLinearOpenCompactHash_Bucket;
intintHash_Table *intintLCGLinearOpenCompactHash_CreateTable(intintHash_Factory
							     * factory,
							     int hashIndex,
							     size_t keyRange,
							     size_t numEntries,
							     float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintLCGLinearOpenCompactHash_DestroyTable;
	table->setupFunc = &intintLCGLinearOpenCompactHash_SetupTable;
	table->emptyFunc = &intintLCGLinearOpenCompactHash_EmptyTable;
	table->queryFunc = &intintLCGLinearOpenCompactHash_Query;
	table->querySingleFunc = &intintLCGLinearOpenCompactHash_QuerySingle;
	table->insertFunc = &intintLCGLinearOpenCompactHash_Insert;
	table->insertSingleFunc = &intintLCGLinearOpenCompactHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintLCGLinearOpenCompactHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintLCGLinearOpenCompactHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(intintLCGLinearOpenCompactHash_TableData));
	((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    hashID = LCG_LINEAR_OPEN_COMPACT_HASH_ID;
	((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    numBuckets;
	char *tempHashData =
	    (char *)malloc(sizeof(intintLCGLinearOpenCompactHash_TableData) +
			   ((intintLCGLinearOpenCompactHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(intintLCGLinearOpenCompactHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintLCGLinearOpenCompactHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintLCGLinearOpenCompactHash_CreateFactory(intintHash_Factory * factory,
						 int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintLCGLinearOpenCompactHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintLCGLinearOpenCompactHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactHash_DestroyFactory(intintHash_Factory * factory,
						  int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintLCGLinearOpenCompactHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (int index = 0;
		     index <
		     ((intintLCGLinearOpenCompactHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGLinearOpenCompactHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	for (int index = 0;
	     index <
	     ((intintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	     numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGLinearOpenCompactHash_InnerQuerySingle(char *tableData, int key,
						    int *valueOutput) {
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	int index;
	int exitCode;
	intintLCGLinearOpenCompactHash_TableData *mytableData =
	    (intintLCGLinearOpenCompactHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((intintLCGLinearOpenCompactHash_TableData *) tableData)->
		     numBuckets);
		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintLCGLinearOpenCompactHash_InnerQuery(char *tableData,
					      unsigned int numKeys, int *keys,
					      int *valuesOutput) {
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		intintLCGLinearOpenCompactHash_TableData *mytableData =
		    (intintLCGLinearOpenCompactHash_TableData *) tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((intintLCGLinearOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGLinearOpenCompactHash_InnerInsertSingle(char *tableData, int key,
						     int value) {
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	int index;
	int exitCode;
	intintLCGLinearOpenCompactHash_TableData *mytableData =
	    (intintLCGLinearOpenCompactHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((intintLCGLinearOpenCompactHash_TableData *) tableData)->
		     numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGLinearOpenCompactHash_InnerInsert(char *tableData,
					       unsigned int numEntries,
					       int *keys, int *values) {
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGLinearOpenCompactHash_TableData *mytableData =
		    (intintLCGLinearOpenCompactHash_TableData *) tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((intintLCGLinearOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGLinearOpenCompactHash_InnerInsertSingleNoOverwrite(char *tableData,
								int key,
								int value) {
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	int index;
	int exitCode;
	intintLCGLinearOpenCompactHash_TableData *mytableData =
	    (intintLCGLinearOpenCompactHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((intintLCGLinearOpenCompactHash_TableData *) tableData)->
		     numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGLinearOpenCompactHash_InnerInsertNoOverwrite(char *tableData,
							  unsigned int
							  numEntries, int *keys,
							  int *values) {
	intintLCGLinearOpenCompactHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGLinearOpenCompactHash_TableData *mytableData =
		    (intintLCGLinearOpenCompactHash_TableData *) tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((intintLCGLinearOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGLinearOpenCompactHash_QuerySingle(intintHash_Table * table,
					       int key, int *valueOutput) {
	return intintLCGLinearOpenCompactHash_InnerQuerySingle(table->tableData,
							       key,
							       valueOutput);
}
int intintLCGLinearOpenCompactHash_Query(intintHash_Table * table,
					 size_t numKeys, int *keys,
					 int *valuesOutput) {
	return intintLCGLinearOpenCompactHash_InnerQuery(table->tableData,
							 numKeys, keys,
							 valuesOutput);
}
int intintLCGLinearOpenCompactHash_InsertSingle(intintHash_Table * table,
						int key, int value) {
	return intintLCGLinearOpenCompactHash_InnerInsertSingle(table->
								tableData, key,
								value);
}
int intintLCGLinearOpenCompactHash_Insert(intintHash_Table * table,
					  size_t numEntries, int *keys,
					  int *values) {
	return intintLCGLinearOpenCompactHash_InnerInsert(table->tableData,
							  numEntries, keys,
							  values);
}
int intintLCGLinearOpenCompactHash_InsertSingleNoOverwrite(intintHash_Table *
							   table, int key,
							   int value) {
	return
	    intintLCGLinearOpenCompactHash_InnerInsertSingleNoOverwrite(table->
									tableData,
									key,
									value);
}
int intintLCGLinearOpenCompactHash_InsertNoOverwrite(intintHash_Table * table,
						     size_t numEntries,
						     int *keys, int *values) {
	return intintLCGLinearOpenCompactHash_InnerInsertNoOverwrite(table->
								     tableData,
								     numEntries,
								     keys,
								     values);
}

typedef struct intintLCGLinearOpenCompactCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	intintHash_CompressLCGData compressFuncData;
} intintLCGLinearOpenCompactCLHash_TableData;
typedef struct intintLCGLinearOpenCompactCLHash_Bucket {
	int key;
	int value;
} intintLCGLinearOpenCompactCLHash_Bucket;
intintHash_Table
    *intintLCGLinearOpenCompactCLHash_CreateTable(intintHash_Factory * factory,
						  int hashIndex,
						  size_t keyRange,
						  size_t numEntries,
						  float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintLCGLinearOpenCompactCLHash_DestroyTable;
	table->setupFunc = &intintLCGLinearOpenCompactCLHash_SetupTable;
	table->emptyFunc = &intintLCGLinearOpenCompactCLHash_EmptyTable;
	table->queryFunc = &intintLCGLinearOpenCompactCLHash_Query;
	table->querySingleFunc = &intintLCGLinearOpenCompactCLHash_QuerySingle;
	table->insertFunc = &intintLCGLinearOpenCompactCLHash_Insert;
	table->insertSingleFunc =
	    &intintLCGLinearOpenCompactCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintLCGLinearOpenCompactCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(intintLCGLinearOpenCompactCLHash_TableData));
	((intintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    hashID = LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((intintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((intintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((intintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((intintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((intintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((intintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    numBuckets;
	char *tempHashData =
	    (char *)malloc(sizeof(intintLCGLinearOpenCompactCLHash_TableData) +
			   ((intintLCGLinearOpenCompactCLHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintLCGLinearOpenCompactCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintLCGLinearOpenCompactCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(intintLCGLinearOpenCompactHash_TableData) +
			   ((intintLCGLinearOpenCompactHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(intintLCGLinearOpenCompactHash_Bucket), NULL,
			   &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0,
				 sizeof
				 (intintLCGLinearOpenCompactHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int intintLCGLinearOpenCompactCLHash_CreateFactory(intintHash_Factory * factory,
						   int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintLCGLinearOpenCompactCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintLCGLinearOpenCompactCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "intintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintLCGLinearOpenCompactCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintLCGLinearOpenCompactCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int intintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct intintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }intintHash_CompressLCGData; static inline unsigned int intintHash_CompressLCG(intintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct intintLCGLinearOpenCompactCLHash_TableData{ int hashID; unsigned int numBuckets; intintHash_CompressLCGData compressFuncData; }intintLCGLinearOpenCompactCLHash_TableData; typedef struct intintLCGLinearOpenCompactCLHash_Bucket{ int key; int value; }intintLCGLinearOpenCompactCLHash_Bucket; __kernel void intintLCGLinearOpenCompactCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global intintLCGLinearOpenCompactCLHash_TableData*)tableData)->numBuckets){ return; } __global intintLCGLinearOpenCompactCLHash_Bucket *buckets = (__global intintLCGLinearOpenCompactCLHash_Bucket*)&tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)]; buckets[index].key = -1;/*HASH_BUCKET_STATUS_EMPTY*/ }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "intintLCGLinearOpenCompactCLHash_Empty", &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "intintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactCLHash_DestroyFactory(intintHash_Factory *
						    factory, int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactCLHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintLCGLinearOpenCompactCLHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintLCGLinearOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintLCGLinearOpenCompactCLHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintLCGLinearOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintLCGLinearOpenCompactCLHash_QuerySingle(intintHash_Table * table,
						 int key, int *valueOutput) {
	return intintLCGLinearOpenCompactCLHash_Query(table, 1, &key,
						      valueOutput);
}
int intintLCGLinearOpenCompactCLHash_Query(intintHash_Table * table,
					   size_t numKeys, int *keys,
					   int *valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(int) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_Query",
				     "clCreateBuffer");
	intintLCGLinearOpenCompactCLHash_BufferQuery(table, numKeys, keysBuffer,
						     valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(int) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactCLHash_BufferQuery(intintHash_Table * table,
						 size_t numKeys,
						 cl_mem keysBuffer,
						 cl_mem valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactCLHash_InsertSingle(intintHash_Table * table,
						  int key, int value) {
	return intintLCGLinearOpenCompactCLHash_Insert(table, 1, &key, &value);
}
int intintLCGLinearOpenCompactCLHash_Insert(intintHash_Table * table,
					    size_t numEntries, int *keys,
					    int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	intintLCGLinearOpenCompactCLHash_BufferInsert(table, numEntries,
						      keysBuffer, valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactCLHash_BufferInsert(intintHash_Table * table,
						  size_t numEntries,
						  cl_mem keysBuffer,
						  cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(intintHash_Table *
							     table, int key,
							     int value) {
	return intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(table, 1,
								  &key, &value);
}
int intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(intintHash_Table * table,
						       size_t numEntries,
						       int *keys, int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite(table,
								 numEntries,
								 keysBuffer,
								 valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite(intintHash_Table *
							     table,
							     size_t numEntries,
							     cl_mem keysBuffer,
							     cl_mem
							     valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct intintLCGLinearOpenCompactOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	intintHash_CompressLCGData compressFuncData;
} intintLCGLinearOpenCompactOpenMPHash_TableData;
typedef struct intintLCGLinearOpenCompactOpenMPHash_Bucket {
	int key;
	int value;
} intintLCGLinearOpenCompactOpenMPHash_Bucket;
intintHash_Table
    *intintLCGLinearOpenCompactOpenMPHash_CreateTable(intintHash_Factory *
						      factory, int hashIndex,
						      size_t keyRange,
						      size_t numEntries,
						      float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintLCGLinearOpenCompactOpenMPHash_DestroyTable;
	table->setupFunc = &intintLCGLinearOpenCompactOpenMPHash_SetupTable;
	table->emptyFunc = &intintLCGLinearOpenCompactOpenMPHash_EmptyTable;
	table->queryFunc = &intintLCGLinearOpenCompactOpenMPHash_Query;
	table->querySingleFunc =
	    &intintLCGLinearOpenCompactOpenMPHash_QuerySingle;
	table->insertFunc = &intintLCGLinearOpenCompactOpenMPHash_Insert;
	table->insertSingleFunc =
	    &intintLCGLinearOpenCompactOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintLCGLinearOpenCompactOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintLCGLinearOpenCompactOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData));
	((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->tableData)->
	    hashID = LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID;
	((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	     tableData)->numBuckets;
	char *tempHashData =
	    (char *)
	    malloc(sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData) +
		   ((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->
		    tableData)->numBuckets *
		   sizeof(intintLCGLinearOpenCompactOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintLCGLinearOpenCompactOpenMPHash_CreateFactory(intintHash_Factory *
						       factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintLCGLinearOpenCompactOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintLCGLinearOpenCompactOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactOpenMPHash_DestroyFactory(intintHash_Factory *
							factory,
							int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGLinearOpenCompactOpenMPHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintLCGLinearOpenCompactOpenMPHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (int index = 0;
		     index <
		     ((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGLinearOpenCompactOpenMPHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
#pragma omp parallel for
	for (int index = 0;
	     index <
	     ((intintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGLinearOpenCompactOpenMPHash_InnerQuerySingle(char *tableData,
							  int key,
							  int *valueOutput) {
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
	int index;
	int exitCode;
	intintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
	    (intintLCGLinearOpenCompactOpenMPHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((intintLCGLinearOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		int old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintLCGLinearOpenCompactOpenMPHash_InnerQuery(char *tableData,
						    unsigned int numKeys,
						    int *keys,
						    int *valuesOutput) {
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		intintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
		    (intintLCGLinearOpenCompactOpenMPHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((intintLCGLinearOpenCompactOpenMPHash_TableData *)
			      tableData)->numBuckets);
			int old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGLinearOpenCompactOpenMPHash_InnerInsertSingle(char *tableData,
							   int key, int value) {
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
	int index;
	int exitCode;
	intintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
	    (intintLCGLinearOpenCompactOpenMPHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((intintLCGLinearOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		int old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGLinearOpenCompactOpenMPHash_InnerInsert(char *tableData,
						     unsigned int numEntries,
						     int *keys, int *values) {
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
		    (intintLCGLinearOpenCompactOpenMPHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((intintLCGLinearOpenCompactOpenMPHash_TableData *)
			      tableData)->numBuckets);
			int old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGLinearOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite(char
								      *tableData,
								      int key,
								      int
								      value) {
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
	int index;
	int exitCode;
	intintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
	    (intintLCGLinearOpenCompactOpenMPHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((intintLCGLinearOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		int old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGLinearOpenCompactOpenMPHash_InnerInsertNoOverwrite(char *tableData,
								unsigned int
								numEntries,
								int *keys,
								int *values) {
	intintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
		    (intintLCGLinearOpenCompactOpenMPHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((intintLCGLinearOpenCompactOpenMPHash_TableData *)
			      tableData)->numBuckets);
			int old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGLinearOpenCompactOpenMPHash_QuerySingle(intintHash_Table * table,
						     int key,
						     int *valueOutput) {
	return intintLCGLinearOpenCompactOpenMPHash_InnerQuerySingle(table->
								     tableData,
								     key,
								     valueOutput);
}
int intintLCGLinearOpenCompactOpenMPHash_Query(intintHash_Table * table,
					       size_t numKeys, int *keys,
					       int *valuesOutput) {
	return intintLCGLinearOpenCompactOpenMPHash_InnerQuery(table->tableData,
							       numKeys, keys,
							       valuesOutput);
}
int intintLCGLinearOpenCompactOpenMPHash_InsertSingle(intintHash_Table * table,
						      int key, int value) {
	return intintLCGLinearOpenCompactOpenMPHash_InnerInsertSingle(table->
								      tableData,
								      key,
								      value);
}
int intintLCGLinearOpenCompactOpenMPHash_Insert(intintHash_Table * table,
						size_t numEntries, int *keys,
						int *values) {
	return intintLCGLinearOpenCompactOpenMPHash_InnerInsert(table->
								tableData,
								numEntries,
								keys, values);
}
int
intintLCGLinearOpenCompactOpenMPHash_InsertSingleNoOverwrite(intintHash_Table *
							     table, int key,
							     int value) {
	return
	    intintLCGLinearOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int intintLCGLinearOpenCompactOpenMPHash_InsertNoOverwrite(intintHash_Table *
							   table,
							   size_t numEntries,
							   int *keys,
							   int *values) {
	return
	    intintLCGLinearOpenCompactOpenMPHash_InnerInsertNoOverwrite(table->
									tableData,
									numEntries,
									keys,
									values);
}

typedef struct intintLCGQuadraticOpenCompactHash_TableData {
	int hashID;
	unsigned int numBuckets;
	intintHash_CompressLCGData compressFuncData;
} intintLCGQuadraticOpenCompactHash_TableData;
typedef struct intintLCGQuadraticOpenCompactHash_Bucket {
	int key;
	int value;
} intintLCGQuadraticOpenCompactHash_Bucket;
intintHash_Table
    *intintLCGQuadraticOpenCompactHash_CreateTable(intintHash_Factory * factory,
						   int hashIndex,
						   size_t keyRange,
						   size_t numEntries,
						   float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintLCGQuadraticOpenCompactHash_DestroyTable;
	table->setupFunc = &intintLCGQuadraticOpenCompactHash_SetupTable;
	table->emptyFunc = &intintLCGQuadraticOpenCompactHash_EmptyTable;
	table->queryFunc = &intintLCGQuadraticOpenCompactHash_Query;
	table->querySingleFunc = &intintLCGQuadraticOpenCompactHash_QuerySingle;
	table->insertFunc = &intintLCGQuadraticOpenCompactHash_Insert;
	table->insertSingleFunc =
	    &intintLCGQuadraticOpenCompactHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintLCGQuadraticOpenCompactHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintLCGQuadraticOpenCompactHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(intintLCGQuadraticOpenCompactHash_TableData));
	((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    hashID = LCG_QUADRATIC_OPEN_COMPACT_HASH_ID;
	((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    numBuckets;
	((intintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    numBuckets =
	    largestProthPrimeUnder(((intintLCGQuadraticOpenCompactHash_TableData
				     *) table->tableData)->numBuckets);
	char *tempHashData =
	    (char *)malloc(sizeof(intintLCGQuadraticOpenCompactHash_TableData) +
			   ((intintLCGQuadraticOpenCompactHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintLCGQuadraticOpenCompactHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintLCGQuadraticOpenCompactHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintLCGQuadraticOpenCompactHash_CreateFactory(intintHash_Factory *
						    factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintLCGQuadraticOpenCompactHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintLCGQuadraticOpenCompactHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactHash_DestroyFactory(intintHash_Factory *
						     factory, int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintLCGQuadraticOpenCompactHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (int index = 0;
		     index <
		     ((intintLCGQuadraticOpenCompactHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGQuadraticOpenCompactHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	for (int index = 0;
	     index <
	     ((intintLCGQuadraticOpenCompactHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGQuadraticOpenCompactHash_InnerQuerySingle(char *tableData, int key,
						       int *valueOutput) {
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	int index;
	int exitCode;
	intintLCGQuadraticOpenCompactHash_TableData *mytableData =
	    (intintLCGQuadraticOpenCompactHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((intintLCGQuadraticOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((intintLCGQuadraticOpenCompactHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintLCGQuadraticOpenCompactHash_InnerQuery(char *tableData,
						 unsigned int numKeys,
						 int *keys, int *valuesOutput) {
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		intintLCGQuadraticOpenCompactHash_TableData *mytableData =
		    (intintLCGQuadraticOpenCompactHash_TableData *) tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((intintLCGQuadraticOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((intintLCGQuadraticOpenCompactHash_TableData
				   *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGQuadraticOpenCompactHash_InnerInsertSingle(char *tableData,
							int key, int value) {
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	int index;
	int exitCode;
	intintLCGQuadraticOpenCompactHash_TableData *mytableData =
	    (intintLCGQuadraticOpenCompactHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((intintLCGQuadraticOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((intintLCGQuadraticOpenCompactHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGQuadraticOpenCompactHash_InnerInsert(char *tableData,
						  unsigned int numEntries,
						  int *keys, int *values) {
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGQuadraticOpenCompactHash_TableData *mytableData =
		    (intintLCGQuadraticOpenCompactHash_TableData *) tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((intintLCGQuadraticOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((intintLCGQuadraticOpenCompactHash_TableData
				   *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGQuadraticOpenCompactHash_InnerInsertSingleNoOverwrite(char
								   *tableData,
								   int key,
								   int value) {
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	int index;
	int exitCode;
	intintLCGQuadraticOpenCompactHash_TableData *mytableData =
	    (intintLCGQuadraticOpenCompactHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((intintLCGQuadraticOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((intintLCGQuadraticOpenCompactHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGQuadraticOpenCompactHash_InnerInsertNoOverwrite(char *tableData,
							     unsigned int
							     numEntries,
							     int *keys,
							     int *values) {
	intintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGQuadraticOpenCompactHash_TableData *mytableData =
		    (intintLCGQuadraticOpenCompactHash_TableData *) tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((intintLCGQuadraticOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((intintLCGQuadraticOpenCompactHash_TableData
				   *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGQuadraticOpenCompactHash_QuerySingle(intintHash_Table * table,
						  int key, int *valueOutput) {
	return intintLCGQuadraticOpenCompactHash_InnerQuerySingle(table->
								  tableData,
								  key,
								  valueOutput);
}
int intintLCGQuadraticOpenCompactHash_Query(intintHash_Table * table,
					    size_t numKeys, int *keys,
					    int *valuesOutput) {
	return intintLCGQuadraticOpenCompactHash_InnerQuery(table->tableData,
							    numKeys, keys,
							    valuesOutput);
}
int intintLCGQuadraticOpenCompactHash_InsertSingle(intintHash_Table * table,
						   int key, int value) {
	return intintLCGQuadraticOpenCompactHash_InnerInsertSingle(table->
								   tableData,
								   key, value);
}
int intintLCGQuadraticOpenCompactHash_Insert(intintHash_Table * table,
					     size_t numEntries, int *keys,
					     int *values) {
	return intintLCGQuadraticOpenCompactHash_InnerInsert(table->tableData,
							     numEntries, keys,
							     values);
}
int intintLCGQuadraticOpenCompactHash_InsertSingleNoOverwrite(intintHash_Table *
							      table, int key,
							      int value) {
	return
	    intintLCGQuadraticOpenCompactHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int intintLCGQuadraticOpenCompactHash_InsertNoOverwrite(intintHash_Table *
							table,
							size_t numEntries,
							int *keys,
							int *values) {
	return intintLCGQuadraticOpenCompactHash_InnerInsertNoOverwrite(table->
									tableData,
									numEntries,
									keys,
									values);
}

typedef struct intintLCGQuadraticOpenCompactCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	intintHash_CompressLCGData compressFuncData;
} intintLCGQuadraticOpenCompactCLHash_TableData;
typedef struct intintLCGQuadraticOpenCompactCLHash_Bucket {
	int key;
	int value;
} intintLCGQuadraticOpenCompactCLHash_Bucket;
intintHash_Table
    *intintLCGQuadraticOpenCompactCLHash_CreateTable(intintHash_Factory *
						     factory, int hashIndex,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc = &intintLCGQuadraticOpenCompactCLHash_DestroyTable;
	table->setupFunc = &intintLCGQuadraticOpenCompactCLHash_SetupTable;
	table->emptyFunc = &intintLCGQuadraticOpenCompactCLHash_EmptyTable;
	table->queryFunc = &intintLCGQuadraticOpenCompactCLHash_Query;
	table->querySingleFunc =
	    &intintLCGQuadraticOpenCompactCLHash_QuerySingle;
	table->insertFunc = &intintLCGQuadraticOpenCompactCLHash_Insert;
	table->insertSingleFunc =
	    &intintLCGQuadraticOpenCompactCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(intintLCGQuadraticOpenCompactCLHash_TableData));
	((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    hashID = LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((intintLCGQuadraticOpenCompactCLHash_TableData *) table->
	     tableData)->numBuckets;
	((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    numBuckets =
	    largestProthPrimeUnder(((intintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->numBuckets);
	char *tempHashData =
	    (char *)malloc(sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)
			   +
			   ((intintLCGQuadraticOpenCompactCLHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintLCGQuadraticOpenCompactCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintLCGQuadraticOpenCompactCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(intintLCGQuadraticOpenCompactHash_TableData) +
			   ((intintLCGQuadraticOpenCompactHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(intintLCGQuadraticOpenCompactHash_Bucket),
			   NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0,
				 sizeof
				 (intintLCGQuadraticOpenCompactHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int intintLCGQuadraticOpenCompactCLHash_CreateFactory(intintHash_Factory *
						      factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintLCGQuadraticOpenCompactCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintLCGQuadraticOpenCompactCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "intintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintLCGQuadraticOpenCompactCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintLCGQuadraticOpenCompactCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "intintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int intintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct intintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }intintHash_CompressLCGData; static inline unsigned int intintHash_CompressLCG(intintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct intintLCGQuadraticOpenCompactCLHash_TableData{ int hashID; unsigned int numBuckets; intintHash_CompressLCGData compressFuncData; }intintLCGQuadraticOpenCompactCLHash_TableData; typedef struct intintLCGQuadraticOpenCompactCLHash_Bucket{ int key; int value; }intintLCGQuadraticOpenCompactCLHash_Bucket; __kernel void intintLCGQuadraticOpenCompactCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global intintLCGQuadraticOpenCompactCLHash_TableData*)tableData)->numBuckets){ return; } __global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets = (__global intintLCGQuadraticOpenCompactCLHash_Bucket*)&tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)]; buckets[index].key = -1;/*HASH_BUCKET_STATUS_EMPTY*/ }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "intintLCGQuadraticOpenCompactCLHash_Empty", &error);
	CLHash_Utilities_HandleError(error,
				     "intintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "intintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactCLHash_DestroyFactory(intintHash_Factory *
						       factory,
						       int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactCLHash_DestroyTable(intintHash_Table * table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintLCGQuadraticOpenCompactCLHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintLCGQuadraticOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintLCGQuadraticOpenCompactCLHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((intintLCGQuadraticOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int intintLCGQuadraticOpenCompactCLHash_QuerySingle(intintHash_Table * table,
						    int key, int *valueOutput) {
	return intintLCGQuadraticOpenCompactCLHash_Query(table, 1, &key,
							 valueOutput);
}
int intintLCGQuadraticOpenCompactCLHash_Query(intintHash_Table * table,
					      size_t numKeys, int *keys,
					      int *valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(int) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_Query",
				     "clCreateBuffer");
	intintLCGQuadraticOpenCompactCLHash_BufferQuery(table, numKeys,
							keysBuffer,
							valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(int) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactCLHash_BufferQuery(intintHash_Table * table,
						    size_t numKeys,
						    cl_mem keysBuffer,
						    cl_mem valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactCLHash_InsertSingle(intintHash_Table * table,
						     int key, int value) {
	return intintLCGQuadraticOpenCompactCLHash_Insert(table, 1, &key,
							  &value);
}
int intintLCGQuadraticOpenCompactCLHash_Insert(intintHash_Table * table,
					       size_t numEntries, int *keys,
					       int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	intintLCGQuadraticOpenCompactCLHash_BufferInsert(table, numEntries,
							 keysBuffer,
							 valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactCLHash_BufferInsert(intintHash_Table * table,
						     size_t numEntries,
						     cl_mem keysBuffer,
						     cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(intintHash_Table
								* table,
								int key,
								int value) {
	return intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(table, 1,
								     &key,
								     &value);
}
int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(intintHash_Table *
							  table,
							  size_t numEntries,
							  int *keys,
							  int *values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(int) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite(table,
								    numEntries,
								    keysBuffer,
								    valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite(intintHash_Table
								* table,
								size_t
								numEntries,
								cl_mem
								keysBuffer,
								cl_mem
								valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct intintLCGQuadraticOpenCompactOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	intintHash_CompressLCGData compressFuncData;
} intintLCGQuadraticOpenCompactOpenMPHash_TableData;
typedef struct intintLCGQuadraticOpenCompactOpenMPHash_Bucket {
	int key;
	int value;
} intintLCGQuadraticOpenCompactOpenMPHash_Bucket;
intintHash_Table
    *intintLCGQuadraticOpenCompactOpenMPHash_CreateTable(intintHash_Factory *
							 factory, int hashIndex,
							 size_t keyRange,
							 size_t numEntries,
							 float loadFactor) {
	intintHash_Table *table =
	    (intintHash_Table *) malloc(sizeof(intintHash_Table));
	table->destroyFunc =
	    &intintLCGQuadraticOpenCompactOpenMPHash_DestroyTable;
	table->setupFunc = &intintLCGQuadraticOpenCompactOpenMPHash_SetupTable;
	table->emptyFunc = &intintLCGQuadraticOpenCompactOpenMPHash_EmptyTable;
	table->queryFunc = &intintLCGQuadraticOpenCompactOpenMPHash_Query;
	table->querySingleFunc =
	    &intintLCGQuadraticOpenCompactOpenMPHash_QuerySingle;
	table->insertFunc = &intintLCGQuadraticOpenCompactOpenMPHash_Insert;
	table->insertSingleFunc =
	    &intintLCGQuadraticOpenCompactOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &intintLCGQuadraticOpenCompactOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &intintLCGQuadraticOpenCompactOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(intintLCGQuadraticOpenCompactOpenMPHash_TableData));
	((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->hashID = LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID;
	((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->numBuckets =
(unsigned int)((double)numEntries / loadFactor);
	((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.a = HASH_LCG_A;
	((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.c = HASH_LCG_C;
	((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.m = HASH_LCG_M;
	((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.n =
((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->tableData)->numBuckets;
	((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->numBuckets =
largestProthPrimeUnder(((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->tableData)->numBuckets);
	char *tempHashData =
	    (char *)
	    malloc(sizeof(intintLCGQuadraticOpenCompactOpenMPHash_TableData) +
		   ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		    table->tableData)->numBuckets *
		   sizeof(intintLCGQuadraticOpenCompactOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(intintLCGQuadraticOpenCompactOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int intintLCGQuadraticOpenCompactOpenMPHash_CreateFactory(intintHash_Factory *
							  factory,
							  int hashIndex) {
	factory->createFunc[hashIndex] =
	    &intintLCGQuadraticOpenCompactOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &intintLCGQuadraticOpenCompactOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactOpenMPHash_DestroyFactory(intintHash_Factory *
							   factory,
							   int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int intintLCGQuadraticOpenCompactOpenMPHash_DestroyTable(intintHash_Table *
							 table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int intintLCGQuadraticOpenCompactOpenMPHash_SetupTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	if (intintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (int index = 0;
		     index <
		     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      table->tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGQuadraticOpenCompactOpenMPHash_EmptyTable(intintHash_Table * table) {
	int exitCode = 0;
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
#pragma omp parallel for
	for (int index = 0;
	     index <
	     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	} exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int intintLCGQuadraticOpenCompactOpenMPHash_InnerQuerySingle(char *tableData,
							     int key,
							     int *valueOutput) {
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int index;
	int exitCode;
	intintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
	    (intintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		int old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintLCGQuadraticOpenCompactOpenMPHash_InnerQuery(char *tableData,
						       unsigned int numKeys,
						       int *keys,
						       int *valuesOutput) {
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int key;
	int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		intintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
		    (intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData
			       *) tableData)->numBuckets);
			int old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingle(char *tableData,
							      int key,
							      int value) {
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int index;
	int exitCode;
	intintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
	    (intintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		int old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGQuadraticOpenCompactOpenMPHash_InnerInsert(char *tableData,
							unsigned int numEntries,
							int *keys,
							int *values) {
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
		    (intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData
			       *) tableData)->numBuckets);
			int old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite(char
									 *tableData,
									 int
									 key,
									 int
									 value) 
{
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int index;
	int exitCode;
	intintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
	    (intintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		int old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int intintLCGQuadraticOpenCompactOpenMPHash_InnerInsertNoOverwrite(char
								   *tableData,
								   unsigned int
								   numEntries,
								   int *keys,
								   int *values) 
{
	intintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (intintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (intintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		intintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
		    (intintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((intintLCGQuadraticOpenCompactOpenMPHash_TableData
			       *) tableData)->numBuckets);
			int old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((intintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int intintLCGQuadraticOpenCompactOpenMPHash_QuerySingle(intintHash_Table *
							table, int key,
							int *valueOutput) {
	return intintLCGQuadraticOpenCompactOpenMPHash_InnerQuerySingle(table->
									tableData,
									key,
									valueOutput);
}
int intintLCGQuadraticOpenCompactOpenMPHash_Query(intintHash_Table * table,
						  size_t numKeys, int *keys,
						  int *valuesOutput) {
	return intintLCGQuadraticOpenCompactOpenMPHash_InnerQuery(table->
								  tableData,
								  numKeys, keys,
								  valuesOutput);
}
int intintLCGQuadraticOpenCompactOpenMPHash_InsertSingle(intintHash_Table *
							 table, int key,
							 int value) {
	return intintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingle(table->
									 tableData,
									 key,
									 value);
}
int intintLCGQuadraticOpenCompactOpenMPHash_Insert(intintHash_Table * table,
						   size_t numEntries, int *keys,
						   int *values) {
	return intintLCGQuadraticOpenCompactOpenMPHash_InnerInsert(table->
								   tableData,
								   numEntries,
								   keys,
								   values);
}
int
intintLCGQuadraticOpenCompactOpenMPHash_InsertSingleNoOverwrite(intintHash_Table
								* table,
								int key,
								int value) {
	return
	    intintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int intintLCGQuadraticOpenCompactOpenMPHash_InsertNoOverwrite(intintHash_Table *
							      table,
							      size_t numEntries,
							      int *keys,
							      int *values) {
	return
	    intintLCGQuadraticOpenCompactOpenMPHash_InnerInsertNoOverwrite
	    (table->tableData, numEntries, keys, values);
}
const char *HashFactory_source =
"\n"
"/* Copyright (C) 1991-2012 Free Software Foundation, Inc.\n"
"   This file is part of the GNU C Library.\n"
"\n"
"   The GNU C Library is free software; you can redistribute it and/or\n"
"   modify it under the terms of the GNU Lesser General Public\n"
"   License as published by the Free Software Foundation; either\n"
"   version 2.1 of the License, or (at your option) any later version.\n"
"\n"
"   The GNU C Library is distributed in the hope that it will be useful,\n"
"   but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n"
"   Lesser General Public License for more details.\n"
"\n"
"   You should have received a copy of the GNU Lesser General Public\n"
"   License along with the GNU C Library; if not, see\n"
"   <http://www.gnu.org/licenses/>.  */\n"
"/* This header is separate from features.h so that the compiler can\n"
"   include it implicitly at the start of every compilation.  It must\n"
"   not itself include <features.h> or any other header that includes\n"
"   <features.h> because the implicit include comes before any feature\n"
"   test macros that may be defined in a source file before it first\n"
"   explicitly includes a system header.  GCC knows the name of this\n"
"   header in order to preinclude it.  */\n"
"/* We do support the IEC 559 math functionality, real and complex.  */\n"
"/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /\n"
"   Unicode 6.0.  */\n"
"/* We do not support C11 <threads.h>.  */\n"
"/* Copyright 2013-14.  Los Alamos National Security, LLC. This material was produced\n"
" * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National \n"
" * Laboratory (LANL), which is operated by Los Alamos National Security, LLC\n"
" * for the U.S. Department of Energy. The U.S. Government has rights to use,\n"
" * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS\n"
" * ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR\n"
" * ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified\n"
" * to produce derivative works, such modified software should be clearly marked,\n"
" * so as not to confuse it with the version available from LANL.   \n"
" *\n"
" * Licensed under the Apache License, Version 2.0 (the ""License""); you may not\n"
" * use this file except in compliance with the License. You may obtain a copy\n"
" * of the License at \n"
" *\n"
" * http://www.apache.org/licenses/LICENSE-2.0\n"
" *\n"
" * Unless required by applicable law or agreed to in writing, software distributed\n"
" * under the License is distributed on an ""AS IS"" BASIS, WITHOUT WARRANTIES OR\n"
" * CONDITIONS OF ANY KIND, either express or implied. See the License for the\n"
" * specific language governing permissions and limitations under the License.\n"
" *\n"
" * Under this license, it is required to include a reference to this work. We\n"
" * request that each derivative work contain a reference to LANL Copyright \n"
" * Disclosure C14043/LA-CC-14-003 so that this work's impact can be roughly\n"
" * measured. In addition, it is requested that a modifier is included as in\n"
" * the following example:\n"
" *\n"
" * //<Uses | improves on | modified from> LANL Copyright Disclosure C14043/LA-CC-14-003\n"
" *\n"
" * This is LANL Copyright Disclosure C14043/LA-CC-14-003\n"
" */\n"
"int intintIdentityPerfectCLHash_InsertSingle(__global char *tableData,\n"
"					     int key, int value);\n"
"int intintIdentityPerfectCLHash_InnerInsertSingle(__global char *tableData,\n"
"						  int key, int value);\n"
"int intintHash_InsertSingle(__global char *tableData, int key, int value);\n"
"int intintIdentityPerfectCLHash_InnerQuery(__global char *tableData,\n"
"					   unsigned int numKeys,\n"
"					   __global int *keys,\n"
"					   __global int *valuesOutput);\n"
"int intintIdentityPerfectCLHash_InnerQuerySingle(__global char *tableData,\n"
"						 int key,\n"
"						 __global int *valueOutput);\n"
"int intintIdentityPerfectCLHash_InnerInsert(__global char *tableData,\n"
"					    unsigned int numEntries,\n"
"					    __global int *keys,\n"
"					    __global int *values);\n"
"int intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(__global char\n"
"							     *tableData,\n"
"							     int key,\n"
"							     int value);\n"
"int intintIdentityPerfectCLHash_InnerInsertNoOverwrite(__global char *tableData,\n"
"						       unsigned int numEntries,\n"
"						       __global int *keys,\n"
"						       __global int *values);\n"
"int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,\n"
"					    __global int *valueOutput);\n"
"int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,\n"
"					    __global int *valueOutput);\n"
"int intintIdentityPerfectCLHash_Query(__global char *tableData, size_t numKeys,\n"
"				      __global int *keys,\n"
"				      __global int *valuesOutput);\n"
"int intintIdentityPerfectCLHash_Insert(__global char *tableData,\n"
"				       size_t numEntries, __global int *keys,\n"
"				       __global int *values);\n"
"int intintIdentityPerfectCLHash_InsertSingleNoOverwrite(__global char\n"
"							*tableData, int key,\n"
"							int value);\n"
"int intintIdentityPerfectCLHash_InsertNoOverwrite(__global char *tableData,\n"
"						  size_t numEntries,\n"
"						  __global int *keys,\n"
"						  __global int *values);\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(__global char\n"
"							       *tableData,\n"
"							       unsigned int\n"
"							       numEntries,\n"
"							       __global int\n"
"							       *keys,\n"
"							       __global int\n"
"							       *values);\n"
"int intintIdentitySentinelPerfectCLHash_InnerQuerySingle(__global char\n"
"							 *tableData, int key,\n"
"							 __global int\n"
"							 *valueOutput);\n"
"int intintIdentitySentinelPerfectCLHash_InnerQuery(__global char *tableData,\n"
"						   unsigned int numKeys,\n"
"						   __global int *keys,\n"
"						   __global int *valuesOutput);\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsertSingle(__global char\n"
"							  *tableData, int key,\n"
"							  int value);\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsert(__global char *tableData,\n"
"						    unsigned int numEntries,\n"
"						    __global int *keys,\n"
"						    __global int *values);\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite(__global\n"
"								     char\n"
"								     *tableData,\n"
"								     int key,\n"
"								     int value);\n"
"int intintIdentitySentinelPerfectCLHash_QuerySingle(__global char *tableData,\n"
"						    int key,\n"
"						    __global int *valueOutput);\n"
"int intintIdentitySentinelPerfectCLHash_Query(__global char *tableData,\n"
"					      size_t numKeys,\n"
"					      __global int *keys,\n"
"					      __global int *valuesOutput);\n"
"int intintIdentitySentinelPerfectCLHash_InsertSingle(__global char *tableData,\n"
"						     int key, int value);\n"
"int intintIdentitySentinelPerfectCLHash_Insert(__global char *tableData,\n"
"					       size_t numEntries,\n"
"					       __global int *keys,\n"
"					       __global int *values);\n"
"int intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(__global char\n"
"								*tableData,\n"
"								int key,\n"
"								int value);\n"
"int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(__global char\n"
"							  *tableData,\n"
"							  size_t numEntries,\n"
"							  __global int *keys,\n"
"							  __global int *values);\n"
"int intintLCGLinearOpenCompactCLHash_InnerQuerySingle(__global char *tableData,\n"
"						      int key,\n"
"						      __global int\n"
"						      *valueOutput);\n"
"int intintLCGLinearOpenCompactCLHash_QuerySingle(__global char *tableData,\n"
"						 int key,\n"
"						 __global int *valueOutput);\n"
"int intintLCGLinearOpenCompactCLHash_Query(__global char *tableData,\n"
"					   size_t numKeys, __global int *keys,\n"
"					   __global int *valuesOutput);\n"
"int intintLCGLinearOpenCompactCLHash_InsertSingle(__global char *tableData,\n"
"						  int key, int value);\n"
"int intintLCGLinearOpenCompactCLHash_Insert(__global char *tableData,\n"
"					    size_t numEntries,\n"
"					    __global int *keys,\n"
"					    __global int *values);\n"
"int intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(__global char\n"
"							     *tableData,\n"
"							     int key,\n"
"							     int value);\n"
"int intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(__global char *tableData,\n"
"						       size_t numEntries,\n"
"						       __global int *keys,\n"
"						       __global int *values);\n"
"int intintLCGLinearOpenCompactCLHash_InnerQuery(__global char *tableData,\n"
"						unsigned int numKeys,\n"
"						__global int *keys,\n"
"						__global int *valuesOutput);\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(__global char\n"
"							    *tableData,\n"
"							    unsigned int\n"
"							    numEntries,\n"
"							    __global int *keys,\n"
"							    __global int\n"
"							    *values);\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsertSingle(__global char *tableData,\n"
"						       int key, int value);\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global char\n"
"								  *tableData,\n"
"								  int key,\n"
"								  int value);\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsert(__global char *tableData,\n"
"						 unsigned int numEntries,\n"
"						 __global int *keys,\n"
"						 __global int *values);\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(__global char\n"
"							 *tableData, int key,\n"
"							 __global int\n"
"							 *valueOutput);\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerQuery(__global char *tableData,\n"
"						   unsigned int numKeys,\n"
"						   __global int *keys,\n"
"						   __global int *valuesOutput);\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(__global char\n"
"							  *tableData, int key,\n"
"							  int value);\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsert(__global char *tableData,\n"
"						    unsigned int numEntries,\n"
"						    __global int *keys,\n"
"						    __global int *values);\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global\n"
"								     char\n"
"								     *tableData,\n"
"								     int key,\n"
"								     int value);\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(__global char\n"
"							       *tableData,\n"
"							       unsigned int\n"
"							       numEntries,\n"
"							       __global int\n"
"							       *keys,\n"
"							       __global int\n"
"							       *values);\n"
"int intintLCGQuadraticOpenCompactCLHash_QuerySingle(__global char *tableData,\n"
"						    int key,\n"
"						    __global int *valueOutput);\n"
"int intintLCGQuadraticOpenCompactCLHash_Query(__global char *tableData,\n"
"					      size_t numKeys,\n"
"					      __global int *keys,\n"
"					      __global int *valuesOutput);\n"
"int intintLCGQuadraticOpenCompactCLHash_InsertSingle(__global char *tableData,\n"
"						     int key, int value);\n"
"int intintLCGQuadraticOpenCompactCLHash_Insert(__global char *tableData,\n"
"					       size_t numEntries,\n"
"					       __global int *keys,\n"
"					       __global int *values);\n"
"int intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(__global char\n"
"								*tableData,\n"
"								int key,\n"
"								int value);\n"
"int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(__global char\n"
"							  *tableData,\n"
"							  size_t numEntries,\n"
"							  __global int *keys,\n"
"							  __global int *values);\n"
"int intintHash_Query(__global char *tableData, unsigned int numKeys,\n"
"		     __global int *keys, __global int *valuesOutput);\n"
"int intintHash_QuerySingle(__global char *tableData, int key,\n"
"			   __global int *valueOutput);\n"
"int intintHash_Insert(__global char *tableData, unsigned int numEntries,\n"
"		      __global int *keys, __global int *values);\n"
"int intintHash_InsertNoOverwrite(__global char *tableData,\n"
"				 unsigned int numEntries, __global int *keys,\n"
"				 __global int *values);\n"
"int intintHash_InsertSingleNoOverwrite(__global char *tableData, int key,\n"
"				       int value);\n"
"#define HASH_REPORT_NEVER /**/ 0\n"
"#define HASH_REPORT_CYCLE /**/ 1\n"
"#define HASH_REPORT_END /****/ 2\n"
"//\n"
"#define HASH_EXIT_CODE_NORMAL /****************/ -1\n"
"#define HASH_EXIT_CODE_ERROR /*****************/ -2\n"
"#define HASH_EXIT_CODE_OVERWRITE /*************/ -3\n"
"#define HASH_EXIT_CODE_KEY_DNE /***************/ -4\n"
"#define HASH_EXIT_CODE_CYCLE /*****************/ -5\n"
"#define HASH_EXIT_CODE_MAX_ENTRIES_EXCEEDED /**/ -6\n"
"#define HASH_EXIT_CODE_BUCKET_INDEX_OOB /******/ -7\n"
"//\n"
"#define HASH_SEARCH_CODE_MATCH /*****/ 0\n"
"#define HASH_SEARCH_CODE_MISMATCH /**/ 1\n"
"#define HASH_SEARCH_CODE_EMPTY /*****/ 2\n"
"//\n"
"#define IDENTITY_PERFECT_CL_HASH_ID /****************/ 16\n"
"#define IDENTITY_SENTINEL_PERFECT_CL_HASH_ID /*******/ 32\n"
"#define LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID /*********/ 64\n"
"#define LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID /******/ 128\n"
"//\n"
"#define HASH_BUCKET_STATUS_EMPTY /**/ -1\n"
"#define HASH_BUCKET_STATUS_FULL /***/ -2\n"
"#define HASH_BUCKET_STATUS_LOCK /***/ -3\n"
"static inline unsigned int intintHash_CompressIdentity(char data, int hashCode) {\n"
"	return hashCode;\n"
"}\n"
"\n"
"typedef struct intintHash_CompressLCGData {\n"
"	long unsigned int a;\n"
"	long unsigned int c;\n"
"	unsigned int m;\n"
"	unsigned int n;\n"
"} intintHash_CompressLCGData;\n"
"static inline unsigned int intintHash_CompressLCG(intintHash_CompressLCGData\n"
"						  compressLCGData,\n"
"						  int hashCode) {\n"
"	return ((compressLCGData.a * hashCode +\n"
"		 compressLCGData.c) % compressLCGData.m) % compressLCGData.n;\n"
"}\n"
"\n"
"typedef struct intintIdentityPerfectCLHash_TableData {\n"
"	int hashID;\n"
"	unsigned int numBuckets;\n"
"	char compressFuncData;\n"
"} intintIdentityPerfectCLHash_TableData;\n"
"typedef struct intintIdentityPerfectCLHash_Bucket {\n"
"	int key;\n"
"	int value;\n"
"} intintIdentityPerfectCLHash_Bucket;\n"
"int intintIdentityPerfectCLHash_InnerQuerySingle(__global char *tableData,\n"
"						 int key,\n"
"						 __global int *valueOutput) {\n"
"	__global intintIdentityPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentityPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	index =\n"
"	    intintHash_CompressIdentity(((__global\n"
"					  intintIdentityPerfectCLHash_TableData\n"
"					  *) tableData)->compressFuncData, key);\n"
"	if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {\n"
"		if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_MISMATCH;\n"
"		}\n"
"	} else {\n"
"		exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"		*valueOutput = buckets[index].value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		return HASH_EXIT_CODE_KEY_DNE;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintIdentityPerfectCLHash_InnerQuery(__global char *tableData,\n"
"					   unsigned int numKeys,\n"
"					   __global int *keys,\n"
"					   __global int *valuesOutput) {\n"
"	__global intintIdentityPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentityPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];\n"
"	int key;\n"
"	__global int *valueOutput;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	for (i = 0; i < numKeys; i++) {\n"
"		key = keys[i];\n"
"		valueOutput = &valuesOutput[i];\n"
"		index =\n"
"		    intintHash_CompressIdentity(((__global\n"
"						  intintIdentityPerfectCLHash_TableData\n"
"						  *) tableData)->\n"
"						compressFuncData, key);\n"
"		if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {\n"
"			if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			} else {\n"
"				exitCode = HASH_SEARCH_CODE_MISMATCH;\n"
"			}\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"			*valueOutput = buckets[index].value;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			resultExitCode = HASH_EXIT_CODE_KEY_DNE;\n"
"			break;\n"
"		default:\n"
"			return exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintIdentityPerfectCLHash_InnerInsertSingle(__global char *tableData,\n"
"						  int key, int value) {\n"
"	__global intintIdentityPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentityPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	index =\n"
"	    intintHash_CompressIdentity(((__global\n"
"					  intintIdentityPerfectCLHash_TableData\n"
"					  *) tableData)->compressFuncData, key);\n"
"	if (((buckets[index].key ==\n"
"	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =\n"
"					   key,\n"
"					   HASH_BUCKET_STATUS_EMPTY) :\n"
"	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {\n"
"		if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_MISMATCH;\n"
"		}\n"
"	} else {\n"
"		exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintIdentityPerfectCLHash_InnerInsert(__global char *tableData,\n"
"					    unsigned int numEntries,\n"
"					    __global int *keys,\n"
"					    __global int *values) {\n"
"	__global intintIdentityPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentityPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		index =\n"
"		    intintHash_CompressIdentity(((__global\n"
"						  intintIdentityPerfectCLHash_TableData\n"
"						  *) tableData)->\n"
"						compressFuncData, key);\n"
"		if (((buckets[index].key ==\n"
"		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =\n"
"						   key,\n"
"						   HASH_BUCKET_STATUS_EMPTY) :\n"
"		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {\n"
"			if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			} else {\n"
"				exitCode = HASH_SEARCH_CODE_MISMATCH;\n"
"			}\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(__global char\n"
"							     *tableData,\n"
"							     int key,\n"
"							     int value) {\n"
"	__global intintIdentityPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentityPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	index =\n"
"	    intintHash_CompressIdentity(((__global\n"
"					  intintIdentityPerfectCLHash_TableData\n"
"					  *) tableData)->compressFuncData, key);\n"
"	if (((buckets[index].key ==\n"
"	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =\n"
"					   key,\n"
"					   HASH_BUCKET_STATUS_EMPTY) :\n"
"	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {\n"
"		if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_MISMATCH;\n"
"		}\n"
"	} else {\n"
"		exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintIdentityPerfectCLHash_InnerInsertNoOverwrite(__global char *tableData,\n"
"						       unsigned int numEntries,\n"
"						       __global int *keys,\n"
"						       __global int *values) {\n"
"	__global intintIdentityPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentityPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		index =\n"
"		    intintHash_CompressIdentity(((__global\n"
"						  intintIdentityPerfectCLHash_TableData\n"
"						  *) tableData)->\n"
"						compressFuncData, key);\n"
"		if (((buckets[index].key ==\n"
"		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =\n"
"						   key,\n"
"						   HASH_BUCKET_STATUS_EMPTY) :\n"
"		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {\n"
"			if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			} else {\n"
"				exitCode = HASH_SEARCH_CODE_MISMATCH;\n"
"			}\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,\n"
"					    __global int *valueOutput) {\n"
"	return intintIdentityPerfectCLHash_InnerQuerySingle(tableData, key,\n"
"							    valueOutput);\n"
"}\n"
"int intintIdentityPerfectCLHash_Query(__global char *tableData, size_t numKeys,\n"
"				      __global int *keys,\n"
"				      __global int *valuesOutput) {\n"
"	return intintIdentityPerfectCLHash_InnerQuery(tableData, numKeys, keys,\n"
"						      valuesOutput);\n"
"}\n"
"int intintIdentityPerfectCLHash_InsertSingle(__global char *tableData, int key,\n"
"					     int value) {\n"
"	return intintIdentityPerfectCLHash_InnerInsertSingle(tableData, key,\n"
"							     value);\n"
"}\n"
"int intintIdentityPerfectCLHash_Insert(__global char *tableData,\n"
"				       size_t numEntries, __global int *keys,\n"
"				       __global int *values) {\n"
"	return intintIdentityPerfectCLHash_InnerInsert(tableData, numEntries,\n"
"						       keys, values);\n"
"}\n"
"int intintIdentityPerfectCLHash_InsertSingleNoOverwrite(__global char\n"
"							*tableData, int key,\n"
"							int value) {\n"
"	return\n"
"	    intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(tableData,\n"
"								     key,\n"
"								     value);\n"
"}\n"
"int intintIdentityPerfectCLHash_InsertNoOverwrite(__global char *tableData,\n"
"						  size_t numEntries,\n"
"						  __global int *keys,\n"
"						  __global int *values) {\n"
"	return intintIdentityPerfectCLHash_InnerInsertNoOverwrite(tableData,\n"
"								  numEntries,\n"
"								  keys, values);\n"
"}\n"
"__kernel void intintIdentityPerfectCLHash_RangeQuerySingle(__global char\n"
"							   *tableData,\n"
"							   unsigned int\n"
"							   numQueries,\n"
"							   __global int *keys,\n"
"							   __global int\n"
"							   *valuesOutput) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintIdentityPerfectCLHash_InnerQuerySingle(tableData, keys[i],\n"
"						     valuesOutput + i);\n"
"}\n"
"__kernel void intintIdentityPerfectCLHash_RangeQuery(__global char *tableData,\n"
"						     unsigned int numQueries,\n"
"						     unsigned int numKeys,\n"
"						     __global int *keys,\n"
"						     __global int\n"
"						     *valuesOutput) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintIdentityPerfectCLHash_InnerQuery(tableData, numKeys,\n"
"					       keys + (i * numKeys),\n"
"					       valuesOutput + (i * numKeys));\n"
"}\n"
"__kernel void intintIdentityPerfectCLHash_RangeInsertSingle(__global char\n"
"							    *tableData,\n"
"							    unsigned int\n"
"							    numInsertions,\n"
"							    __global int *keys,\n"
"							    __global int\n"
"							    *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentityPerfectCLHash_InnerInsertSingle(tableData, keys[i],\n"
"						      values[i]);\n"
"}\n"
"__kernel void intintIdentityPerfectCLHash_RangeInsert(__global char *tableData,\n"
"						      unsigned int\n"
"						      numInsertions,\n"
"						      unsigned int numEntries,\n"
"						      __global int *keys,\n"
"						      __global int *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentityPerfectCLHash_InnerInsert(tableData, numEntries,\n"
"						keys + (i * numEntries),\n"
"						values + (i * numEntries));\n"
"}\n"
"__kernel void intintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite(__global\n"
"								       char\n"
"								       *tableData,\n"
"								       unsigned\n"
"								       int\n"
"								       numInsertions,\n"
"								       __global\n"
"								       int\n"
"								       *keys,\n"
"								       __global\n"
"								       int\n"
"								       *values) \n"
"{\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(tableData,\n"
"								 keys[i],\n"
"								 values[i]);\n"
"}\n"
"__kernel void intintIdentityPerfectCLHash_RangeInsertNoOverwrite(__global char\n"
"								 *tableData,\n"
"								 unsigned int\n"
"								 numInsertions,\n"
"								 unsigned int\n"
"								 numEntries,\n"
"								 __global int\n"
"								 *keys,\n"
"								 __global int\n"
"								 *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentityPerfectCLHash_InnerInsertNoOverwrite(tableData,\n"
"							   numEntries,\n"
"							   keys +\n"
"							   (i * numEntries),\n"
"							   values +\n"
"							   (i * numEntries));\n"
"}\n"
"\n"
"typedef struct intintIdentitySentinelPerfectCLHash_TableData {\n"
"	int hashID;\n"
"	unsigned int numBuckets;\n"
"	char compressFuncData;\n"
"	int emptyValue;\n"
"} intintIdentitySentinelPerfectCLHash_TableData;\n"
"typedef struct intintIdentitySentinelPerfectCLHash_Bucket {\n"
"	int value;\n"
"} intintIdentitySentinelPerfectCLHash_Bucket;\n"
"int intintIdentitySentinelPerfectCLHash_InnerQuerySingle(__global char\n"
"							 *tableData, int key,\n"
"							 __global int\n"
"							 *valueOutput) {\n"
"	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	index =\n"
"	    intintHash_CompressIdentity(((__global\n"
"					  intintIdentitySentinelPerfectCLHash_TableData\n"
"					  *) tableData)->compressFuncData, key);\n"
"	if (buckets[index].value !=\n"
"	    ((__global intintIdentitySentinelPerfectCLHash_TableData *)\n"
"	     tableData)->emptyValue) {\n"
"		exitCode = HASH_SEARCH_CODE_MATCH;\n"
"	} else {\n"
"		exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"		*valueOutput = buckets[index].value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		return HASH_EXIT_CODE_KEY_DNE;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InnerQuery(__global char *tableData,\n"
"						   unsigned int numKeys,\n"
"						   __global int *keys,\n"
"						   __global int *valuesOutput) {\n"
"	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];\n"
"	int key;\n"
"	__global int *valueOutput;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	for (i = 0; i < numKeys; i++) {\n"
"		key = keys[i];\n"
"		valueOutput = &valuesOutput[i];\n"
"		index =\n"
"		    intintHash_CompressIdentity(((__global\n"
"						  intintIdentitySentinelPerfectCLHash_TableData\n"
"						  *) tableData)->\n"
"						compressFuncData, key);\n"
"		if (buckets[index].value !=\n"
"		    ((__global intintIdentitySentinelPerfectCLHash_TableData *)\n"
"		     tableData)->emptyValue) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"			*valueOutput = buckets[index].value;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			resultExitCode = HASH_EXIT_CODE_KEY_DNE;\n"
"			break;\n"
"		default:\n"
"			return exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsertSingle(__global char\n"
"							  *tableData, int key,\n"
"							  int value) {\n"
"	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	index =\n"
"	    intintHash_CompressIdentity(((__global\n"
"					  intintIdentitySentinelPerfectCLHash_TableData\n"
"					  *) tableData)->compressFuncData, key);\n"
"	if (buckets[index].value !=\n"
"	    ((__global intintIdentitySentinelPerfectCLHash_TableData *)\n"
"	     tableData)->emptyValue) {\n"
"		exitCode = HASH_SEARCH_CODE_MATCH;\n"
"	} else {\n"
"		exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsert(__global char *tableData,\n"
"						    unsigned int numEntries,\n"
"						    __global int *keys,\n"
"						    __global int *values) {\n"
"	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		index =\n"
"		    intintHash_CompressIdentity(((__global\n"
"						  intintIdentitySentinelPerfectCLHash_TableData\n"
"						  *) tableData)->\n"
"						compressFuncData, key);\n"
"		if (buckets[index].value !=\n"
"		    ((__global intintIdentitySentinelPerfectCLHash_TableData *)\n"
"		     tableData)->emptyValue) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite(__global\n"
"								     char\n"
"								     *tableData,\n"
"								     int key,\n"
"								     int value) \n"
"{\n"
"	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	index =\n"
"	    intintHash_CompressIdentity(((__global\n"
"					  intintIdentitySentinelPerfectCLHash_TableData\n"
"					  *) tableData)->compressFuncData, key);\n"
"	if (buckets[index].value !=\n"
"	    ((__global intintIdentitySentinelPerfectCLHash_TableData *)\n"
"	     tableData)->emptyValue) {\n"
"		exitCode = HASH_SEARCH_CODE_MATCH;\n"
"	} else {\n"
"		exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(__global char\n"
"							       *tableData,\n"
"							       unsigned int\n"
"							       numEntries,\n"
"							       __global int\n"
"							       *keys,\n"
"							       __global int\n"
"							       *values) {\n"
"	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =\n"
"	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		index =\n"
"		    intintHash_CompressIdentity(((__global\n"
"						  intintIdentitySentinelPerfectCLHash_TableData\n"
"						  *) tableData)->\n"
"						compressFuncData, key);\n"
"		if (buckets[index].value !=\n"
"		    ((__global intintIdentitySentinelPerfectCLHash_TableData *)\n"
"		     tableData)->emptyValue) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"		} else {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_QuerySingle(__global char *tableData,\n"
"						    int key,\n"
"						    __global int *valueOutput) {\n"
"	return intintIdentitySentinelPerfectCLHash_InnerQuerySingle(tableData,\n"
"								    key,\n"
"								    valueOutput);\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_Query(__global char *tableData,\n"
"					      size_t numKeys,\n"
"					      __global int *keys,\n"
"					      __global int *valuesOutput) {\n"
"	return intintIdentitySentinelPerfectCLHash_InnerQuery(tableData,\n"
"							      numKeys, keys,\n"
"							      valuesOutput);\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InsertSingle(__global char *tableData,\n"
"						     int key, int value) {\n"
"	return intintIdentitySentinelPerfectCLHash_InnerInsertSingle(tableData,\n"
"								     key,\n"
"								     value);\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_Insert(__global char *tableData,\n"
"					       size_t numEntries,\n"
"					       __global int *keys,\n"
"					       __global int *values) {\n"
"	return intintIdentitySentinelPerfectCLHash_InnerInsert(tableData,\n"
"							       numEntries, keys,\n"
"							       values);\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(__global char\n"
"								*tableData,\n"
"								int key,\n"
"								int value) {\n"
"	return\n"
"	    intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite\n"
"	    (tableData, key, value);\n"
"}\n"
"int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(__global char\n"
"							  *tableData,\n"
"							  size_t numEntries,\n"
"							  __global int *keys,\n"
"							  __global int *values) \n"
"{\n"
"	return\n"
"	    intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite\n"
"	    (tableData, numEntries, keys, values);\n"
"}\n"
"__kernel void intintIdentitySentinelPerfectCLHash_RangeQuerySingle(__global char\n"
"								   *tableData,\n"
"								   unsigned int\n"
"								   numQueries,\n"
"								   __global int\n"
"								   *keys,\n"
"								   __global int\n"
"								   *valuesOutput) \n"
"{\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintIdentitySentinelPerfectCLHash_InnerQuerySingle(tableData, keys[i],\n"
"							     valuesOutput + i);\n"
"}\n"
"__kernel void intintIdentitySentinelPerfectCLHash_RangeQuery(__global char\n"
"							     *tableData,\n"
"							     unsigned int\n"
"							     numQueries,\n"
"							     unsigned int\n"
"							     numKeys,\n"
"							     __global int *keys,\n"
"							     __global int\n"
"							     *valuesOutput) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintIdentitySentinelPerfectCLHash_InnerQuery(tableData, numKeys,\n"
"						       keys + (i * numKeys),\n"
"						       valuesOutput +\n"
"						       (i * numKeys));\n"
"}\n"
"__kernel void intintIdentitySentinelPerfectCLHash_RangeInsertSingle(__global\n"
"								    char\n"
"								    *tableData,\n"
"								    unsigned int\n"
"								    numInsertions,\n"
"								    __global int\n"
"								    *keys,\n"
"								    __global int\n"
"								    *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentitySentinelPerfectCLHash_InnerInsertSingle(tableData,\n"
"							      keys[i],\n"
"							      values[i]);\n"
"}\n"
"__kernel void intintIdentitySentinelPerfectCLHash_RangeInsert(__global char\n"
"							      *tableData,\n"
"							      unsigned int\n"
"							      numInsertions,\n"
"							      unsigned int\n"
"							      numEntries,\n"
"							      __global int\n"
"							      *keys,\n"
"							      __global int\n"
"							      *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentitySentinelPerfectCLHash_InnerInsert(tableData, numEntries,\n"
"							keys + (i * numEntries),\n"
"							values +\n"
"							(i * numEntries));\n"
"}\n"
"__kernel void\n"
"intintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite(__global char\n"
"								 *tableData,\n"
"								 unsigned int\n"
"								 numInsertions,\n"
"								 __global int\n"
"								 *keys,\n"
"								 __global int\n"
"								 *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite\n"
"	    (tableData, keys[i], values[i]);\n"
"}\n"
"__kernel void\n"
"intintIdentitySentinelPerfectCLHash_RangeInsertNoOverwrite(__global char\n"
"							   *tableData,\n"
"							   unsigned int\n"
"							   numInsertions,\n"
"							   unsigned int\n"
"							   numEntries,\n"
"							   __global int *keys,\n"
"							   __global int\n"
"							   *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(tableData,\n"
"								   numEntries,\n"
"								   keys +\n"
"								   (i *\n"
"								    numEntries),\n"
"								   values +\n"
"								   (i *\n"
"								    numEntries));\n"
"}\n"
"\n"
"typedef struct intintLCGLinearOpenCompactCLHash_TableData {\n"
"	int hashID;\n"
"	unsigned int numBuckets;\n"
"	intintHash_CompressLCGData compressFuncData;\n"
"} intintLCGLinearOpenCompactCLHash_TableData;\n"
"typedef struct intintLCGLinearOpenCompactCLHash_Bucket {\n"
"	int key;\n"
"	int value;\n"
"} intintLCGLinearOpenCompactCLHash_Bucket;\n"
"int intintLCGLinearOpenCompactCLHash_InnerQuerySingle(__global char *tableData,\n"
"						      int key,\n"
"						      __global int\n"
"						      *valueOutput) {\n"
"	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	__global intintLCGLinearOpenCompactCLHash_TableData *mytableData =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_TableData *) tableData;\n"
"	intintHash_CompressLCGData compressFuncData =\n"
"	    mytableData->compressFuncData;\n"
"	unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"	unsigned long int iteration = 0;\n"
"	for (;;) {\n"
"		index =\n"
"		    ((1 * iteration +\n"
"		      c) %\n"
"		     ((__global intintLCGLinearOpenCompactCLHash_TableData *)\n"
"		      tableData)->numBuckets);\n"
"		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"			break;\n"
"		} else if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			break;\n"
"		} else if ((index == c && iteration > 0)) {\n"
"			exitCode = HASH_EXIT_CODE_CYCLE;\n"
"			break;\n"
"		}\n"
"		iteration++;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"		*valueOutput = buckets[index].value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		return HASH_EXIT_CODE_KEY_DNE;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InnerQuery(__global char *tableData,\n"
"						unsigned int numKeys,\n"
"						__global int *keys,\n"
"						__global int *valuesOutput) {\n"
"	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];\n"
"	int key;\n"
"	__global int *valueOutput;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	for (i = 0; i < numKeys; i++) {\n"
"		key = keys[i];\n"
"		valueOutput = &valuesOutput[i];\n"
"		__global intintLCGLinearOpenCompactCLHash_TableData *mytableData\n"
"		    =\n"
"		    (__global intintLCGLinearOpenCompactCLHash_TableData *)\n"
"		    tableData;\n"
"		intintHash_CompressLCGData compressFuncData =\n"
"		    mytableData->compressFuncData;\n"
"		unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"		unsigned long int iteration = 0;\n"
"		for (;;) {\n"
"			index =\n"
"			    ((1 * iteration +\n"
"			      c) %\n"
"			     ((__global\n"
"			       intintLCGLinearOpenCompactCLHash_TableData *)\n"
"			      tableData)->numBuckets);\n"
"			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {\n"
"				exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"				break;\n"
"			} else if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"				break;\n"
"			} else if ((index == c && iteration > 0)) {\n"
"				exitCode = HASH_EXIT_CODE_CYCLE;\n"
"				break;\n"
"			}\n"
"			iteration++;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"			*valueOutput = buckets[index].value;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			resultExitCode = HASH_EXIT_CODE_KEY_DNE;\n"
"			break;\n"
"		default:\n"
"			return exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsertSingle(__global char *tableData,\n"
"						       int key, int value) {\n"
"	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	__global intintLCGLinearOpenCompactCLHash_TableData *mytableData =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_TableData *) tableData;\n"
"	intintHash_CompressLCGData compressFuncData =\n"
"	    mytableData->compressFuncData;\n"
"	unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"	unsigned long int iteration = 0;\n"
"	for (;;) {\n"
"		index =\n"
"		    ((1 * iteration +\n"
"		      c) %\n"
"		     ((__global intintLCGLinearOpenCompactCLHash_TableData *)\n"
"		      tableData)->numBuckets);\n"
"		if ((atomic_cmpxchg\n"
"		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"		      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"			break;\n"
"		} else if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			break;\n"
"		} else if ((index == c && iteration > 0)) {\n"
"			exitCode = HASH_EXIT_CODE_CYCLE;\n"
"			break;\n"
"		}\n"
"		iteration++;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsert(__global char *tableData,\n"
"						 unsigned int numEntries,\n"
"						 __global int *keys,\n"
"						 __global int *values) {\n"
"	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		__global intintLCGLinearOpenCompactCLHash_TableData *mytableData\n"
"		    =\n"
"		    (__global intintLCGLinearOpenCompactCLHash_TableData *)\n"
"		    tableData;\n"
"		intintHash_CompressLCGData compressFuncData =\n"
"		    mytableData->compressFuncData;\n"
"		unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"		unsigned long int iteration = 0;\n"
"		for (;;) {\n"
"			index =\n"
"			    ((1 * iteration +\n"
"			      c) %\n"
"			     ((__global\n"
"			       intintLCGLinearOpenCompactCLHash_TableData *)\n"
"			      tableData)->numBuckets);\n"
"			if ((atomic_cmpxchg\n"
"			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"			      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"				exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"				break;\n"
"			} else if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"				break;\n"
"			} else if ((index == c && iteration > 0)) {\n"
"				exitCode = HASH_EXIT_CODE_CYCLE;\n"
"				break;\n"
"			}\n"
"			iteration++;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global char\n"
"								  *tableData,\n"
"								  int key,\n"
"								  int value) {\n"
"	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	__global intintLCGLinearOpenCompactCLHash_TableData *mytableData =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_TableData *) tableData;\n"
"	intintHash_CompressLCGData compressFuncData =\n"
"	    mytableData->compressFuncData;\n"
"	unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"	unsigned long int iteration = 0;\n"
"	for (;;) {\n"
"		index =\n"
"		    ((1 * iteration +\n"
"		      c) %\n"
"		     ((__global intintLCGLinearOpenCompactCLHash_TableData *)\n"
"		      tableData)->numBuckets);\n"
"		if ((atomic_cmpxchg\n"
"		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"		      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"			break;\n"
"		} else if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			break;\n"
"		} else if ((index == c && iteration > 0)) {\n"
"			exitCode = HASH_EXIT_CODE_CYCLE;\n"
"			break;\n"
"		}\n"
"		iteration++;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(__global char\n"
"							    *tableData,\n"
"							    unsigned int\n"
"							    numEntries,\n"
"							    __global int *keys,\n"
"							    __global int\n"
"							    *values) {\n"
"	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		__global intintLCGLinearOpenCompactCLHash_TableData *mytableData\n"
"		    =\n"
"		    (__global intintLCGLinearOpenCompactCLHash_TableData *)\n"
"		    tableData;\n"
"		intintHash_CompressLCGData compressFuncData =\n"
"		    mytableData->compressFuncData;\n"
"		unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"		unsigned long int iteration = 0;\n"
"		for (;;) {\n"
"			index =\n"
"			    ((1 * iteration +\n"
"			      c) %\n"
"			     ((__global\n"
"			       intintLCGLinearOpenCompactCLHash_TableData *)\n"
"			      tableData)->numBuckets);\n"
"			if ((atomic_cmpxchg\n"
"			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"			      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"				exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"				break;\n"
"			} else if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"				break;\n"
"			} else if ((index == c && iteration > 0)) {\n"
"				exitCode = HASH_EXIT_CODE_CYCLE;\n"
"				break;\n"
"			}\n"
"			iteration++;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_QuerySingle(__global char *tableData,\n"
"						 int key,\n"
"						 __global int *valueOutput) {\n"
"	return intintLCGLinearOpenCompactCLHash_InnerQuerySingle(tableData, key,\n"
"								 valueOutput);\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_Query(__global char *tableData,\n"
"					   size_t numKeys, __global int *keys,\n"
"					   __global int *valuesOutput) {\n"
"	return intintLCGLinearOpenCompactCLHash_InnerQuery(tableData, numKeys,\n"
"							   keys, valuesOutput);\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InsertSingle(__global char *tableData,\n"
"						  int key, int value) {\n"
"	return intintLCGLinearOpenCompactCLHash_InnerInsertSingle(tableData,\n"
"								  key, value);\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_Insert(__global char *tableData,\n"
"					    size_t numEntries,\n"
"					    __global int *keys,\n"
"					    __global int *values) {\n"
"	return intintLCGLinearOpenCompactCLHash_InnerInsert(tableData,\n"
"							    numEntries, keys,\n"
"							    values);\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(__global char\n"
"							     *tableData,\n"
"							     int key,\n"
"							     int value) {\n"
"	return\n"
"	    intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite\n"
"	    (tableData, key, value);\n"
"}\n"
"int intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(__global char *tableData,\n"
"						       size_t numEntries,\n"
"						       __global int *keys,\n"
"						       __global int *values) {\n"
"	return\n"
"	    intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(tableData,\n"
"								    numEntries,\n"
"								    keys,\n"
"								    values);\n"
"}\n"
"__kernel void intintLCGLinearOpenCompactCLHash_RangeQuerySingle(__global char\n"
"								*tableData,\n"
"								unsigned int\n"
"								numQueries,\n"
"								__global int\n"
"								*keys,\n"
"								__global int\n"
"								*valuesOutput) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintLCGLinearOpenCompactCLHash_InnerQuerySingle(tableData, keys[i],\n"
"							  valuesOutput + i);\n"
"}\n"
"__kernel void intintLCGLinearOpenCompactCLHash_RangeQuery(__global char\n"
"							  *tableData,\n"
"							  unsigned int\n"
"							  numQueries,\n"
"							  unsigned int numKeys,\n"
"							  __global int *keys,\n"
"							  __global int\n"
"							  *valuesOutput) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintLCGLinearOpenCompactCLHash_InnerQuery(tableData, numKeys,\n"
"						    keys + (i * numKeys),\n"
"						    valuesOutput +\n"
"						    (i * numKeys));\n"
"}\n"
"__kernel void intintLCGLinearOpenCompactCLHash_RangeInsertSingle(__global char\n"
"								 *tableData,\n"
"								 unsigned int\n"
"								 numInsertions,\n"
"								 __global int\n"
"								 *keys,\n"
"								 __global int\n"
"								 *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGLinearOpenCompactCLHash_InnerInsertSingle(tableData, keys[i],\n"
"							   values[i]);\n"
"}\n"
"__kernel void intintLCGLinearOpenCompactCLHash_RangeInsert(__global char\n"
"							   *tableData,\n"
"							   unsigned int\n"
"							   numInsertions,\n"
"							   unsigned int\n"
"							   numEntries,\n"
"							   __global int *keys,\n"
"							   __global int\n"
"							   *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGLinearOpenCompactCLHash_InnerInsert(tableData, numEntries,\n"
"						     keys + (i * numEntries),\n"
"						     values + (i * numEntries));\n"
"}\n"
"__kernel void\n"
"intintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite(__global char\n"
"							      *tableData,\n"
"							      unsigned int\n"
"							      numInsertions,\n"
"							      __global int\n"
"							      *keys,\n"
"							      __global int\n"
"							      *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(tableData,\n"
"								      keys[i],\n"
"								      values\n"
"								      [i]);\n"
"}\n"
"__kernel void intintLCGLinearOpenCompactCLHash_RangeInsertNoOverwrite(__global\n"
"								      char\n"
"								      *tableData,\n"
"								      unsigned\n"
"								      int\n"
"								      numInsertions,\n"
"								      unsigned\n"
"								      int\n"
"								      numEntries,\n"
"								      __global\n"
"								      int *keys,\n"
"								      __global\n"
"								      int\n"
"								      *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(tableData,\n"
"								numEntries,\n"
"								keys +\n"
"								(i *\n"
"								 numEntries),\n"
"								values +\n"
"								(i *\n"
"								 numEntries));\n"
"}\n"
"\n"
"typedef struct intintLCGQuadraticOpenCompactCLHash_TableData {\n"
"	int hashID;\n"
"	unsigned int numBuckets;\n"
"	intintHash_CompressLCGData compressFuncData;\n"
"} intintLCGQuadraticOpenCompactCLHash_TableData;\n"
"typedef struct intintLCGQuadraticOpenCompactCLHash_Bucket {\n"
"	int key;\n"
"	int value;\n"
"} intintLCGQuadraticOpenCompactCLHash_Bucket;\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(__global char\n"
"							 *tableData, int key,\n"
"							 __global int\n"
"							 *valueOutput) {\n"
"	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	__global intintLCGQuadraticOpenCompactCLHash_TableData *mytableData =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"	    tableData;\n"
"	intintHash_CompressLCGData compressFuncData =\n"
"	    mytableData->compressFuncData;\n"
"	unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"	unsigned long int iteration = 0;\n"
"	for (;;) {\n"
"		index =\n"
"		    ((1 * iteration * iteration + 0 * iteration +\n"
"		      c) %\n"
"		     ((__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"		      tableData)->numBuckets);\n"
"		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"			break;\n"
"		} else if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			break;\n"
"		} else\n"
"		    if ((iteration >\n"
"			 ((__global\n"
"			   intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"			  tableData)->numBuckets)) {\n"
"			exitCode = HASH_EXIT_CODE_CYCLE;\n"
"			break;\n"
"		}\n"
"		iteration++;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"		*valueOutput = buckets[index].value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		return HASH_EXIT_CODE_KEY_DNE;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerQuery(__global char *tableData,\n"
"						   unsigned int numKeys,\n"
"						   __global int *keys,\n"
"						   __global int *valuesOutput) {\n"
"	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];\n"
"	int key;\n"
"	__global int *valueOutput;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	for (i = 0; i < numKeys; i++) {\n"
"		key = keys[i];\n"
"		valueOutput = &valuesOutput[i];\n"
"		__global intintLCGQuadraticOpenCompactCLHash_TableData\n"
"		    *mytableData =\n"
"		    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"		    tableData;\n"
"		intintHash_CompressLCGData compressFuncData =\n"
"		    mytableData->compressFuncData;\n"
"		unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"		unsigned long int iteration = 0;\n"
"		for (;;) {\n"
"			index =\n"
"			    ((1 * iteration * iteration + 0 * iteration +\n"
"			      c) %\n"
"			     ((__global\n"
"			       intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"			      tableData)->numBuckets);\n"
"			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {\n"
"				exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"				break;\n"
"			} else if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"				break;\n"
"			} else\n"
"			    if ((iteration >\n"
"				 ((__global\n"
"				   intintLCGQuadraticOpenCompactCLHash_TableData\n"
"				   *) tableData)->numBuckets)) {\n"
"				exitCode = HASH_EXIT_CODE_CYCLE;\n"
"				break;\n"
"			}\n"
"			iteration++;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"			*valueOutput = buckets[index].value;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			resultExitCode = HASH_EXIT_CODE_KEY_DNE;\n"
"			break;\n"
"		default:\n"
"			return exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(__global char\n"
"							  *tableData, int key,\n"
"							  int value) {\n"
"	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	__global intintLCGQuadraticOpenCompactCLHash_TableData *mytableData =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"	    tableData;\n"
"	intintHash_CompressLCGData compressFuncData =\n"
"	    mytableData->compressFuncData;\n"
"	unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"	unsigned long int iteration = 0;\n"
"	for (;;) {\n"
"		index =\n"
"		    ((1 * iteration * iteration + 0 * iteration +\n"
"		      c) %\n"
"		     ((__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"		      tableData)->numBuckets);\n"
"		if ((atomic_cmpxchg\n"
"		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"		      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"			break;\n"
"		} else if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			break;\n"
"		} else\n"
"		    if ((iteration >\n"
"			 ((__global\n"
"			   intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"			  tableData)->numBuckets)) {\n"
"			exitCode = HASH_EXIT_CODE_CYCLE;\n"
"			break;\n"
"		}\n"
"		iteration++;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsert(__global char *tableData,\n"
"						    unsigned int numEntries,\n"
"						    __global int *keys,\n"
"						    __global int *values) {\n"
"	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		__global intintLCGQuadraticOpenCompactCLHash_TableData\n"
"		    *mytableData =\n"
"		    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"		    tableData;\n"
"		intintHash_CompressLCGData compressFuncData =\n"
"		    mytableData->compressFuncData;\n"
"		unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"		unsigned long int iteration = 0;\n"
"		for (;;) {\n"
"			index =\n"
"			    ((1 * iteration * iteration + 0 * iteration +\n"
"			      c) %\n"
"			     ((__global\n"
"			       intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"			      tableData)->numBuckets);\n"
"			if ((atomic_cmpxchg\n"
"			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"			      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"				exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"				break;\n"
"			} else if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"				break;\n"
"			} else\n"
"			    if ((iteration >\n"
"				 ((__global\n"
"				   intintLCGQuadraticOpenCompactCLHash_TableData\n"
"				   *) tableData)->numBuckets)) {\n"
"				exitCode = HASH_EXIT_CODE_CYCLE;\n"
"				break;\n"
"			}\n"
"			iteration++;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global\n"
"								     char\n"
"								     *tableData,\n"
"								     int key,\n"
"								     int value) \n"
"{\n"
"	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];\n"
"	int index;\n"
"	int exitCode;\n"
"	__global intintLCGQuadraticOpenCompactCLHash_TableData *mytableData =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"	    tableData;\n"
"	intintHash_CompressLCGData compressFuncData =\n"
"	    mytableData->compressFuncData;\n"
"	unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"	unsigned long int iteration = 0;\n"
"	for (;;) {\n"
"		index =\n"
"		    ((1 * iteration * iteration + 0 * iteration +\n"
"		      c) %\n"
"		     ((__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"		      tableData)->numBuckets);\n"
"		if ((atomic_cmpxchg\n"
"		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"		      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"			exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"			break;\n"
"		} else if (key == buckets[index].key) {\n"
"			exitCode = HASH_SEARCH_CODE_MATCH;\n"
"			break;\n"
"		} else\n"
"		    if ((iteration >\n"
"			 ((__global\n"
"			   intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"			  tableData)->numBuckets)) {\n"
"			exitCode = HASH_EXIT_CODE_CYCLE;\n"
"			break;\n"
"		}\n"
"		iteration++;\n"
"	}\n"
"	switch (exitCode) {\n"
"	case HASH_SEARCH_CODE_MATCH:\n"
"	case HASH_SEARCH_CODE_MISMATCH:\n"
"		return HASH_EXIT_CODE_OVERWRITE;\n"
"	case HASH_SEARCH_CODE_EMPTY:\n"
"		buckets[index].value = value;\n"
"		return HASH_EXIT_CODE_NORMAL;\n"
"	default:\n"
"		return exitCode;\n"
"	}\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(__global char\n"
"							       *tableData,\n"
"							       unsigned int\n"
"							       numEntries,\n"
"							       __global int\n"
"							       *keys,\n"
"							       __global int\n"
"							       *values) {\n"
"	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =\n"
"	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &\n"
"	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];\n"
"	int resultExitCode = HASH_EXIT_CODE_NORMAL;\n"
"	int key;\n"
"	int index;\n"
"	int exitCode;\n"
"	uint i;;\n"
"	for (i = 0; i < numEntries; i++) {\n"
"		key = keys[i];\n"
"		__global intintLCGQuadraticOpenCompactCLHash_TableData\n"
"		    *mytableData =\n"
"		    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"		    tableData;\n"
"		intintHash_CompressLCGData compressFuncData =\n"
"		    mytableData->compressFuncData;\n"
"		unsigned int c = intintHash_CompressLCG(compressFuncData, key);\n"
"		unsigned long int iteration = 0;\n"
"		for (;;) {\n"
"			index =\n"
"			    ((1 * iteration * iteration + 0 * iteration +\n"
"			      c) %\n"
"			     ((__global\n"
"			       intintLCGQuadraticOpenCompactCLHash_TableData *)\n"
"			      tableData)->numBuckets);\n"
"			if ((atomic_cmpxchg\n"
"			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,\n"
"			      key)) == HASH_BUCKET_STATUS_EMPTY) {\n"
"				exitCode = HASH_SEARCH_CODE_EMPTY;\n"
"				break;\n"
"			} else if (key == buckets[index].key) {\n"
"				exitCode = HASH_SEARCH_CODE_MATCH;\n"
"				break;\n"
"			} else\n"
"			    if ((iteration >\n"
"				 ((__global\n"
"				   intintLCGQuadraticOpenCompactCLHash_TableData\n"
"				   *) tableData)->numBuckets)) {\n"
"				exitCode = HASH_EXIT_CODE_CYCLE;\n"
"				break;\n"
"			}\n"
"			iteration++;\n"
"		}\n"
"		switch (exitCode) {\n"
"		case HASH_SEARCH_CODE_MATCH:\n"
"		case HASH_SEARCH_CODE_MISMATCH:\n"
"			resultExitCode = HASH_EXIT_CODE_OVERWRITE;\n"
"			break;\n"
"		case HASH_SEARCH_CODE_EMPTY:\n"
"			buckets[index].value = values[i];\n"
"			break;\n"
"		default:\n"
"			resultExitCode = exitCode;\n"
"		}\n"
"	}\n"
"	return resultExitCode;\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_QuerySingle(__global char *tableData,\n"
"						    int key,\n"
"						    __global int *valueOutput) {\n"
"	return intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(tableData,\n"
"								    key,\n"
"								    valueOutput);\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_Query(__global char *tableData,\n"
"					      size_t numKeys,\n"
"					      __global int *keys,\n"
"					      __global int *valuesOutput) {\n"
"	return intintLCGQuadraticOpenCompactCLHash_InnerQuery(tableData,\n"
"							      numKeys, keys,\n"
"							      valuesOutput);\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InsertSingle(__global char *tableData,\n"
"						     int key, int value) {\n"
"	return intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(tableData,\n"
"								     key,\n"
"								     value);\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_Insert(__global char *tableData,\n"
"					       size_t numEntries,\n"
"					       __global int *keys,\n"
"					       __global int *values) {\n"
"	return intintLCGQuadraticOpenCompactCLHash_InnerInsert(tableData,\n"
"							       numEntries, keys,\n"
"							       values);\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(__global char\n"
"								*tableData,\n"
"								int key,\n"
"								int value) {\n"
"	return\n"
"	    intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite\n"
"	    (tableData, key, value);\n"
"}\n"
"int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(__global char\n"
"							  *tableData,\n"
"							  size_t numEntries,\n"
"							  __global int *keys,\n"
"							  __global int *values) \n"
"{\n"
"	return\n"
"	    intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite\n"
"	    (tableData, numEntries, keys, values);\n"
"}\n"
"__kernel void intintLCGQuadraticOpenCompactCLHash_RangeQuerySingle(__global char\n"
"								   *tableData,\n"
"								   unsigned int\n"
"								   numQueries,\n"
"								   __global int\n"
"								   *keys,\n"
"								   __global int\n"
"								   *valuesOutput) \n"
"{\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(tableData, keys[i],\n"
"							     valuesOutput + i);\n"
"}\n"
"__kernel void intintLCGQuadraticOpenCompactCLHash_RangeQuery(__global char\n"
"							     *tableData,\n"
"							     unsigned int\n"
"							     numQueries,\n"
"							     unsigned int\n"
"							     numKeys,\n"
"							     __global int *keys,\n"
"							     __global int\n"
"							     *valuesOutput) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numQueries) {\n"
"		return;\n"
"	}\n"
"	intintLCGQuadraticOpenCompactCLHash_InnerQuery(tableData, numKeys,\n"
"						       keys + (i * numKeys),\n"
"						       valuesOutput +\n"
"						       (i * numKeys));\n"
"}\n"
"__kernel void intintLCGQuadraticOpenCompactCLHash_RangeInsertSingle(__global\n"
"								    char\n"
"								    *tableData,\n"
"								    unsigned int\n"
"								    numInsertions,\n"
"								    __global int\n"
"								    *keys,\n"
"								    __global int\n"
"								    *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(tableData,\n"
"							      keys[i],\n"
"							      values[i]);\n"
"}\n"
"__kernel void intintLCGQuadraticOpenCompactCLHash_RangeInsert(__global char\n"
"							      *tableData,\n"
"							      unsigned int\n"
"							      numInsertions,\n"
"							      unsigned int\n"
"							      numEntries,\n"
"							      __global int\n"
"							      *keys,\n"
"							      __global int\n"
"							      *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGQuadraticOpenCompactCLHash_InnerInsert(tableData, numEntries,\n"
"							keys + (i * numEntries),\n"
"							values +\n"
"							(i * numEntries));\n"
"}\n"
"__kernel void\n"
"intintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite(__global char\n"
"								 *tableData,\n"
"								 unsigned int\n"
"								 numInsertions,\n"
"								 __global int\n"
"								 *keys,\n"
"								 __global int\n"
"								 *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite\n"
"	    (tableData, keys[i], values[i]);\n"
"}\n"
"__kernel void\n"
"intintLCGQuadraticOpenCompactCLHash_RangeInsertNoOverwrite(__global char\n"
"							   *tableData,\n"
"							   unsigned int\n"
"							   numInsertions,\n"
"							   unsigned int\n"
"							   numEntries,\n"
"							   __global int *keys,\n"
"							   __global int\n"
"							   *values) {\n"
"	uint i = get_global_id(0);\n"
"	if (i >= numInsertions) {\n"
"		return;\n"
"	}\n"
"	intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(tableData,\n"
"								   numEntries,\n"
"								   keys +\n"
"								   (i *\n"
"								    numEntries),\n"
"								   values +\n"
"								   (i *\n"
"								    numEntries));\n"
"}\n"
"__kernel void intintHash_RangeQuery(__global char *tableData,\n"
"				    unsigned int numQueries,\n"
"				    unsigned int numKeys, __global int *keys,\n"
"				    __global int *valuesOutput) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_RangeQuery(tableData,\n"
"							      numQueries,\n"
"							      numKeys, keys,\n"
"							      valuesOutput);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return intintIdentitySentinelPerfectCLHash_RangeQuery(tableData,\n"
"								      numQueries,\n"
"								      numKeys,\n"
"								      keys,\n"
"								      valuesOutput);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGLinearOpenCompactCLHash_RangeQuery(tableData,\n"
"								   numQueries,\n"
"								   numKeys,\n"
"								   keys,\n"
"								   valuesOutput);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGQuadraticOpenCompactCLHash_RangeQuery(tableData,\n"
"								      numQueries,\n"
"								      numKeys,\n"
"								      keys,\n"
"								      valuesOutput);\n"
"	}\n"
"}\n"
"__kernel void intintHash_RangeQuerySingle(__global char *tableData,\n"
"					  unsigned int numQueries,\n"
"					  __global int *keys,\n"
"					  __global int *valueOutput) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_RangeQuerySingle(tableData,\n"
"								    numQueries,\n"
"								    keys,\n"
"								    valueOutput);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_RangeQuerySingle\n"
"		    (tableData, numQueries, keys, valueOutput);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGLinearOpenCompactCLHash_RangeQuerySingle(tableData,\n"
"								      numQueries,\n"
"								      keys,\n"
"								      valueOutput);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_RangeQuerySingle\n"
"		    (tableData, numQueries, keys, valueOutput);\n"
"	}\n"
"}\n"
"__kernel void intintHash_RangeInsert(__global char *tableData,\n"
"				     unsigned int numInsertions,\n"
"				     unsigned int numEntries,\n"
"				     __global int *keys, __global int *values) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_RangeInsert(tableData,\n"
"							       numInsertions,\n"
"							       numEntries, keys,\n"
"							       values);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_RangeInsert(tableData,\n"
"								    numInsertions,\n"
"								    numEntries,\n"
"								    keys,\n"
"								    values);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGLinearOpenCompactCLHash_RangeInsert(tableData,\n"
"								    numInsertions,\n"
"								    numEntries,\n"
"								    keys,\n"
"								    values);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_RangeInsert(tableData,\n"
"								    numInsertions,\n"
"								    numEntries,\n"
"								    keys,\n"
"								    values);\n"
"	}\n"
"}\n"
"__kernel void intintHash_RangeInsertSingle(__global char *tableData,\n"
"					   unsigned int numInsertions,\n"
"					   __global int *keys,\n"
"					   __global int *values) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_RangeInsertSingle(tableData,\n"
"								     numInsertions,\n"
"								     keys,\n"
"								     values);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_RangeInsertSingle\n"
"		    (tableData, numInsertions, keys, values);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGLinearOpenCompactCLHash_RangeInsertSingle\n"
"		    (tableData, numInsertions, keys, values);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_RangeInsertSingle\n"
"		    (tableData, numInsertions, keys, values);\n"
"	}\n"
"}\n"
"__kernel void intintHash_RangeInsertNoOverwrite(__global char *tableData,\n"
"						unsigned int numInsertions,\n"
"						unsigned int numEntries,\n"
"						__global int *keys,\n"
"						__global int *values) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentityPerfectCLHash_RangeInsertNoOverwrite\n"
"		    (tableData, numInsertions, numEntries, keys, values);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_RangeInsertNoOverwrite\n"
"		    (tableData, numInsertions, numEntries, keys, values);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGLinearOpenCompactCLHash_RangeInsertNoOverwrite\n"
"		    (tableData, numInsertions, numEntries, keys, values);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_RangeInsertNoOverwrite\n"
"		    (tableData, numInsertions, numEntries, keys, values);\n"
"	}\n"
"}\n"
"__kernel void intintHash_RangeInsertSingleNoOverwrite(__global char *tableData,\n"
"						      unsigned int\n"
"						      numInsertions,\n"
"						      __global int *keys,\n"
"						      __global int *values) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite\n"
"		    (tableData, numInsertions, keys, values);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite\n"
"		    (tableData, numInsertions, keys, values);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite\n"
"		    (tableData, numInsertions, keys, values);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite\n"
"		    (tableData, numInsertions, keys, values);\n"
"	}\n"
"}\n"
"int intintHash_Query(__global char *tableData, unsigned int numKeys,\n"
"		     __global int *keys, __global int *valuesOutput) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_Query(tableData, numKeys,\n"
"							 keys, valuesOutput);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return intintIdentitySentinelPerfectCLHash_Query(tableData,\n"
"								 numKeys, keys,\n"
"								 valuesOutput);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGLinearOpenCompactCLHash_Query(tableData,\n"
"							      numKeys, keys,\n"
"							      valuesOutput);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGQuadraticOpenCompactCLHash_Query(tableData,\n"
"								 numKeys, keys,\n"
"								 valuesOutput);\n"
"	}\n"
"	return HASH_EXIT_CODE_ERROR;\n"
"}\n"
"int intintHash_QuerySingle(__global char *tableData, int key,\n"
"			   __global int *valueOutput) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_QuerySingle(tableData, key,\n"
"							       valueOutput);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_QuerySingle(tableData,\n"
"								    key,\n"
"								    valueOutput);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGLinearOpenCompactCLHash_QuerySingle(tableData,\n"
"								    key,\n"
"								    valueOutput);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_QuerySingle(tableData,\n"
"								    key,\n"
"								    valueOutput);\n"
"	}\n"
"	return HASH_EXIT_CODE_ERROR;\n"
"}\n"
"int intintHash_Insert(__global char *tableData, unsigned int numEntries,\n"
"		      __global int *keys, __global int *values) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_Insert(tableData, numEntries,\n"
"							  keys, values);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return intintIdentitySentinelPerfectCLHash_Insert(tableData,\n"
"								  numEntries,\n"
"								  keys, values);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGLinearOpenCompactCLHash_Insert(tableData,\n"
"							       numEntries, keys,\n"
"							       values);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGQuadraticOpenCompactCLHash_Insert(tableData,\n"
"								  numEntries,\n"
"								  keys, values);\n"
"	}\n"
"	return HASH_EXIT_CODE_ERROR;\n"
"}\n"
"int intintHash_InsertSingle(__global char *tableData, int key, int value) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_InsertSingle(tableData, key,\n"
"								value);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_InsertSingle(tableData,\n"
"								     key,\n"
"								     value);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return intintLCGLinearOpenCompactCLHash_InsertSingle(tableData,\n"
"								     key,\n"
"								     value);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_InsertSingle(tableData,\n"
"								     key,\n"
"								     value);\n"
"	}\n"
"	return HASH_EXIT_CODE_ERROR;\n"
"}\n"
"int intintHash_InsertNoOverwrite(__global char *tableData,\n"
"				 unsigned int numEntries, __global int *keys,\n"
"				 __global int *values) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return intintIdentityPerfectCLHash_InsertNoOverwrite(tableData,\n"
"								     numEntries,\n"
"								     keys,\n"
"								     values);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_InsertNoOverwrite\n"
"		    (tableData, numEntries, keys, values);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGLinearOpenCompactCLHash_InsertNoOverwrite\n"
"		    (tableData, numEntries, keys, values);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite\n"
"		    (tableData, numEntries, keys, values);\n"
"	}\n"
"	return HASH_EXIT_CODE_ERROR;\n"
"}\n"
"int intintHash_InsertSingleNoOverwrite(__global char *tableData, int key,\n"
"				       int value) {\n"
"	switch (((__global int *)tableData)[0]) {\n"
"	case IDENTITY_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentityPerfectCLHash_InsertSingleNoOverwrite\n"
"		    (tableData, key, value);\n"
"	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:\n"
"		return\n"
"		    intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite\n"
"		    (tableData, key, value);\n"
"	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite\n"
"		    (tableData, key, value);\n"
"	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:\n"
"		return\n"
"		    intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite\n"
"		    (tableData, key, value);\n"
"	}\n"
"	return HASH_EXIT_CODE_ERROR;\n"
"}\n"
;
