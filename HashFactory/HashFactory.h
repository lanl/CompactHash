
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
 * @file   HashFactory.h
 * @author Peter Ahrens
 * @date   Thu Jun 6 2013 
 */
//

#ifndef HASHFACTORY_H
#define HASHFACTORY_H
//
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
//
#ifndef UINT_TYPE
#define UINT_TYPE
typedef unsigned int uint;
#endif
//
#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#else
typedef void *cl_mem;
typedef int cl_context;
typedef int cl_command_queue;
typedef int cl_program;
typedef int cl_device_id;
typedef int cl_int;
int clRetainContext(int context);
int clRetainCommandQueue(int command_queue);
int clGetContextInfo(int context, int param, size_t size, void *value,
		     size_t * size_ret);
int clReleaseContext(int context);
int clReleaseCommandQueue(int command_queue);
int clReleaseProgram(int program);
int clRetainProgram(int program);
int clRetainKernel(int kernel);
cl_mem clCreateBuffer(int context, int flags, size_t size, void *value,
		      int *size_ret);
int clEnqueueWriteBuffer(int command_queue, void *buffer, int blocking_write,
			 size_t offset, size_t cb, const void *ptr,
			 uint nevents, const int *wait_list, int *event);
int clEnqueueReadBuffer(int command_queue, void *buffer, int blocking_write,
			size_t offset, size_t cb, const void *ptr, uint nevents,
			const int *wait_list, int *event);
int clCreateKernel(int program, const char *kernel_name, int *errcode_ret);
int clReleaseKernel(int kernel);
int clReleaseMemObject(void *memobj);
int clSetKernelArg(int kernel, uint arg_index, size_t arg_size,
		   const void *arg_value);
int clGetKernelWorkGroupInfo(int kernel, int device, int param_name,
			     size_t size, void *value, size_t * size_ret);
int clEnqueueNDRangeKernel(int command_queue, int kernel, uint work_dim,
			   const size_t * offset, const size_t * size,
			   const size_t * local_size, uint nevents,
			   const int *wait_list, int *event);
int clFinish(int command_queue);
#endif
//
#include "CLHash_Utilities.h"
//
//
#include "HashFactory.hm"
#include "HashFactory.cm"
//
#define HASH_REPORT_NEVER /**/ 0
#define HASH_REPORT_CYCLE /**/ 1
#define HASH_REPORT_END /****/ 2
//
#define HASH_EXIT_CODE_NORMAL /****************/ -1
#define HASH_EXIT_CODE_ERROR /*****************/ -2
#define HASH_EXIT_CODE_OVERWRITE /*************/ -3
#define HASH_EXIT_CODE_KEY_DNE /***************/ -4
#define HASH_EXIT_CODE_CYCLE /*****************/ -5
#define HASH_EXIT_CODE_MAX_ENTRIES_EXCEEDED /**/ -6
#define HASH_EXIT_CODE_BUCKET_INDEX_OOB /******/ -7
//
#define HASH_SEARCH_CODE_MATCH /*****/ 0
#define HASH_SEARCH_CODE_MISMATCH /**/ 1
#define HASH_SEARCH_CODE_EMPTY /*****/ 2
//
#define HASH_NUM_HASHES /********/ 12
#define HASH_NUM_C_HASHES /*******/ 4
#define HASH_NUM_CL_HASHES /******/ 4
#define HASH_NUM_OPENMP_HASHES /**/ 4
//
#define IDENTITY_PERFECT_HASH_ID /*******************/ 1
#define IDENTITY_SENTINEL_PERFECT_HASH_ID /**********/ 2
#define LCG_LINEAR_OPEN_COMPACT_HASH_ID /************/ 4
#define LCG_QUADRATIC_OPEN_COMPACT_HASH_ID /*********/ 8
//
#define IDENTITY_PERFECT_CL_HASH_ID /****************/ 16
#define IDENTITY_SENTINEL_PERFECT_CL_HASH_ID /*******/ 32
#define LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID /*********/ 64
#define LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID /******/ 128
//
#define IDENTITY_PERFECT_OPENMP_HASH_ID /***********/ 256
#define IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID /***/ 512
#define LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID /*****/ 1024
#define LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID /**/ 2048
//
/**
 * HASH_ALL_C_HASHES hash table types that run on the host.
 */
#define HASH_ALL_C_HASHES (IDENTITY_PERFECT_HASH_ID | IDENTITY_SENTINEL_PERFECT_HASH_ID | LCG_LINEAR_OPEN_COMPACT_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_HASH_ID)
/**
 * HASH_ALL_CL_HASHES hash table types that run in an OpenCL context.
 */
#define HASH_ALL_CL_HASHES (IDENTITY_PERFECT_CL_HASH_ID | IDENTITY_SENTINEL_PERFECT_CL_HASH_ID | LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID)
/**
 * HASH_ALL_OPENMP_HASHES hash table types that run with OpenMP.
 */
#define HASH_ALL_OPENMP_HASHES (IDENTITY_PERFECT_OPENMP_HASH_ID | IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID | LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID)
/**
 * HASH_ALL_HASHES all hash table types.
 */
#define HASH_ALL_HASHES (HASH_ALL_C_HASHES | HASH_ALL_CL_HASHES | HASH_ALL_OPENMP_HASHES)
/**
 * HASH_NOSENTINEL_PERFECT_HASHES perfect hash table types that do not use a sentinel value to mark empty buckets.
 */
#define HASH_NOSENTINEL_PERFECT_HASHES (IDENTITY_PERFECT_HASH_ID | IDENTITY_PERFECT_CL_HASH_ID | IDENTITY_PERFECT_OPENMP_HASH_ID)
/**
 * HASH_SENTINEL_PERFECT_HASHES perfect hash table types that use a sentinel value to mark empty buckets.
 */
#define HASH_SENTINEL_PERFECT_HASHES (IDENTITY_SENTINEL_PERFECT_HASH_ID | IDENTITY_SENTINEL_PERFECT_CL_HASH_ID | IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID)
/**
 * HASH_PERFECT_HASHES perfect hash table types.
 */
#define HASH_PERFECT_HASHES (HASH_NOSENTINEL_PERFECT_HASHES | HASH_SENTINEL_PERFECT_HASHES)
/**
 * HASH_QUADRATIC_COMPACT_HASHES compact hash table types that use a quadratic probe sequence.
 */
#define HASH_QUADRATIC_COMPACT_HASHES (LCG_QUADRATIC_OPEN_COMPACT_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID)
/**
 * HASH_LINEAR_COMPACT_HASHES compact hash table types that use a linear probe sequence.
 */
#define HASH_LINEAR_COMPACT_HASHES (LCG_LINEAR_OPEN_COMPACT_HASH_ID | LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID | LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID)
/**
 * HASH_COMPACT_HASHES compact hash table types that use a linear probe sequence.
 */
#define HASH_COMPACT_HASHES (HASH_QUADRATIC_COMPACT_HASHES | HASH_LINEAR_COMPACT_HASHES)
//
#define HASH_DEFAULT_LOCAL_WORK_SIZE /*******/ 64
#define HASH_DEFAULT_LOAD_FACTOR /*******/ 0.3
#define HASH_MIN_LOAD_FACTOR /*******/ 0.0009
#define HASH_PERFECT_COMPACT_SWITCH_FACTOR /*******/ 20
#define HASH_LCG_A /*******/ 2147483629
#define HASH_LCG_C /*******/ 2147483587
#define HASH_LCG_M /*******/ 2147483647
//
#define HASH_BUCKET_STATUS_EMPTY /**/ -1
#define HASH_BUCKET_STATUS_FULL /***/ -2
#define HASH_BUCKET_STATUS_LOCK /***/ -3

/**
 * Hash_ExitCodeString will return a string representation of the given exit
 * code.
 * 
 * @param exitCode
 * 
 * @return A string representation of that exit code.
 */
char *Hash_ExitCodeString(int exitCode);

/**
 * Hash_ExitCodeDebug will print a string representation of the given exit code
 * if it is not EXIT_CODE_NORMAL.
 * 
 * @param exitCode
 */
void Hash_ExitCodeDebug(int exitCode);

/**
 * Hash_SetReportLevel sets a static report level variable in hash.c. It should 
 * be called before hash tables are created. 
 *
 * @param level The level of data collection desired.
 */
void Hash_SetReportLevel(int level);

/**
 * Hash_GetReportLevel gets this variable.
 *
 * @return The current report level.
 *         Special Values: Meaning
 *         HASH_REPORT_NEVER: The hash table will not collect data. This is the
 *                       default. A call to Hash_Report will return an empty 
 *                       string.
 *         HASH_REPORT_END: The hash table will collect data and calls to
 *                     Hash_Report will return summary information.
 *         HASH_REPORT_CYCLE: The hash table will collect data and calls to 
 *                       Hash_Report will return information related to the
 *                       last important call.
 */
int Hash_GetReportLevel();

const char *Hash_GetKernelSourceString();
int smallestProthPrimeAbove(int N);
int largestProthPrimeUnder(int N);

typedef struct intintHash_Table_ intintHash_Table;
typedef struct intintCLHash_Table_ intintCLHash_Table;
typedef struct intintHash_Factory_ intintHash_Factory;
typedef struct intintCLHash_Factory_ intintCLHash_Factory;
intintHash_Factory *intintHash_CreateFactory(int HashTypes, int *emptyValue,
					     size_t localWorkSize,
					     cl_context * context,
					     cl_command_queue * queue);
int intintHash_DestroyFactory(intintHash_Factory * factory);
intintHash_Table *intintHash_CreateTable(intintHash_Factory * factory,
					 int hashTypes, size_t keyRange,
					 size_t numEntries, float loadFactor);
int intintHash_EmptyTable(intintHash_Table * table);
int intintHash_DestroyTable(intintHash_Table * table);
cl_mem intintHash_GetTableDataBuffer(intintHash_Table * table);
cl_mem *intintHash_GetTableDataBufferPtr(intintHash_Table * table);
int intintHash_GetTableType(intintHash_Table * table);
int intintHash_Query(intintHash_Table * table, size_t numKeys, int *keys,
		     int *valuesOutput);
int intintHash_QuerySingle(intintHash_Table * table, int key, int *valueOutput);
int intintHash_Insert(intintHash_Table * table, size_t numEntries, int *keys,
		      int *values);
int intintHash_InsertSingle(intintHash_Table * table, int key, int value);
int intintHash_InsertNoOverwrite(intintHash_Table * table, size_t numEntries,
				 int *keys, int *values);
int intintHash_InsertSingleNoOverwrite(intintHash_Table * table, int key,
				       int value);
int intintHash_BufferQuery(intintHash_Table * table, size_t numKeys,
			   cl_mem keys, cl_mem valuesOutput);
int intintHash_BufferInsert(intintHash_Table * table, size_t numEntries,
			    cl_mem keys, cl_mem values);
int intintHash_BufferInsertNoOverwrite(intintHash_Table * table,
				       size_t numEntries, cl_mem keys,
				       cl_mem values);
static inline unsigned int intintHash_CompressIdentity(char data, int hashCode) {
	return hashCode;
}

typedef struct intintHash_CompressLCGData {
	long unsigned int a;
	long unsigned int c;
	unsigned int m;
	unsigned int n;
} intintHash_CompressLCGData;
static inline unsigned int intintHash_CompressLCG(intintHash_CompressLCGData
						  compressLCGData,
						  int hashCode) {
	return ((compressLCGData.a * hashCode +
		 compressLCGData.c) % compressLCGData.m) % compressLCGData.n;
}

int intintIdentityPerfectHash_CreateFactory(intintHash_Factory * factory,
					    int hashIndex);
int intintIdentityPerfectHash_DestroyFactory(intintHash_Factory * factory,
					     int hashIndex);
intintHash_Table *intintIdentityPerfectHash_CreateTable(intintHash_Factory *
							factory, int hashIndex,
							size_t keyRange,
							size_t numEntries,
							float loadFactor);
int intintIdentityPerfectHash_InitTable(intintHash_Table * table, va_list args);
int intintIdentityPerfectHash_DestroyTable(intintHash_Table * table);
char *intintIdentityPerfectHash_Report(intintHash_Table * table);
int intintIdentityPerfectHash_EmptyTable(intintHash_Table * table);
int intintIdentityPerfectHash_Query(intintHash_Table * table, size_t numKeys,
				    int *keys, int *valuesOutput);
int intintIdentityPerfectHash_QuerySingle(intintHash_Table * table, int key,
					  int *valueOutput);
int intintIdentityPerfectHash_Insert(intintHash_Table * table,
				     size_t numEntries, int *keys, int *values);
int intintIdentityPerfectHash_InsertSingle(intintHash_Table * table, int key,
					   int value);
int intintIdentityPerfectHash_InsertNoOverwrite(intintHash_Table * table,
						size_t numEntries, int *keys,
						int *values);
int intintIdentityPerfectHash_InsertSingleNoOverwrite(intintHash_Table * table,
						      int key, int value);
int intintIdentityPerfectCLHash_CreateFactory(intintHash_Factory * factory,
					      int hashIndex);
int intintIdentityPerfectCLHash_DestroyFactory(intintHash_Factory * factory,
					       int hashIndex);
intintHash_Table *intintIdentityPerfectCLHash_CreateTable(intintHash_Factory *
							  factory,
							  int hashIndex,
							  size_t keyRange,
							  size_t numEntries,
							  float loadFactor);
int intintIdentityPerfectCLHash_InitTable(intintHash_Table * table,
					  va_list args);
int intintIdentityPerfectCLHash_DestroyTable(intintHash_Table * table);
char *intintIdentityPerfectCLHash_Report(intintHash_Table * table);
int intintIdentityPerfectCLHash_EmptyTable(intintHash_Table * table);
int intintIdentityPerfectCLHash_Query(intintHash_Table * table, size_t numKeys,
				      int *keys, int *valuesOutput);
int intintIdentityPerfectCLHash_QuerySingle(intintHash_Table * table, int key,
					    int *valueOutput);
int intintIdentityPerfectCLHash_Insert(intintHash_Table * table,
				       size_t numEntries, int *keys,
				       int *values);
int intintIdentityPerfectCLHash_InsertSingle(intintHash_Table * table, int key,
					     int value);
int intintIdentityPerfectCLHash_InsertNoOverwrite(intintHash_Table * table,
						  size_t numEntries, int *keys,
						  int *values);
int intintIdentityPerfectCLHash_InsertSingleNoOverwrite(intintHash_Table *
							table, int key,
							int value);
int intintIdentityPerfectCLHash_BufferQuery(intintHash_Table * table,
					    size_t numKeys, cl_mem keysBuffer,
					    cl_mem valuesOutputBuffer);
int intintIdentityPerfectCLHash_BufferInsert(intintHash_Table * table,
					     size_t numEntries,
					     cl_mem keysBuffer,
					     cl_mem valuesBuffer);
int intintIdentityPerfectCLHash_BufferInsertNoOverwrite(intintHash_Table *
							table,
							size_t numEntries,
							cl_mem keysBuffer,
							cl_mem valuesBuffer);
int intintIdentitySentinelPerfectHash_CreateFactory(intintHash_Factory *
						    factory, int hashIndex);
int intintIdentitySentinelPerfectHash_DestroyFactory(intintHash_Factory *
						     factory, int hashIndex);
intintHash_Table
    *intintIdentitySentinelPerfectHash_CreateTable(intintHash_Factory * factory,
						   int hashIndex,
						   size_t keyRange,
						   size_t numEntries,
						   float loadFactor);
int intintIdentitySentinelPerfectHash_InitTable(intintHash_Table * table,
						va_list args);
int intintIdentitySentinelPerfectHash_DestroyTable(intintHash_Table * table);
char *intintIdentitySentinelPerfectHash_Report(intintHash_Table * table);
int intintIdentitySentinelPerfectHash_EmptyTable(intintHash_Table * table);
int intintIdentitySentinelPerfectHash_Query(intintHash_Table * table,
					    size_t numKeys, int *keys,
					    int *valuesOutput);
int intintIdentitySentinelPerfectHash_QuerySingle(intintHash_Table * table,
						  int key, int *valueOutput);
int intintIdentitySentinelPerfectHash_Insert(intintHash_Table * table,
					     size_t numEntries, int *keys,
					     int *values);
int intintIdentitySentinelPerfectHash_InsertSingle(intintHash_Table * table,
						   int key, int value);
int intintIdentitySentinelPerfectHash_InsertNoOverwrite(intintHash_Table *
							table,
							size_t numEntries,
							int *keys, int *values);
int intintIdentitySentinelPerfectHash_InsertSingleNoOverwrite(intintHash_Table *
							      table, int key,
							      int value);
int intintIdentitySentinelPerfectCLHash_CreateFactory(intintHash_Factory *
						      factory, int hashIndex);
int intintIdentitySentinelPerfectCLHash_DestroyFactory(intintHash_Factory *
						       factory, int hashIndex);
intintHash_Table
    *intintIdentitySentinelPerfectCLHash_CreateTable(intintHash_Factory *
						     factory, int hashIndex,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor);
int intintIdentitySentinelPerfectCLHash_InitTable(intintHash_Table * table,
						  va_list args);
int intintIdentitySentinelPerfectCLHash_DestroyTable(intintHash_Table * table);
char *intintIdentitySentinelPerfectCLHash_Report(intintHash_Table * table);
int intintIdentitySentinelPerfectCLHash_EmptyTable(intintHash_Table * table);
int intintIdentitySentinelPerfectCLHash_Query(intintHash_Table * table,
					      size_t numKeys, int *keys,
					      int *valuesOutput);
int intintIdentitySentinelPerfectCLHash_QuerySingle(intintHash_Table * table,
						    int key, int *valueOutput);
int intintIdentitySentinelPerfectCLHash_Insert(intintHash_Table * table,
					       size_t numEntries, int *keys,
					       int *values);
int intintIdentitySentinelPerfectCLHash_InsertSingle(intintHash_Table * table,
						     int key, int value);
int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(intintHash_Table *
							  table,
							  size_t numEntries,
							  int *keys,
							  int *values);
int intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(intintHash_Table
								* table,
								int key,
								int value);
int intintIdentitySentinelPerfectCLHash_BufferQuery(intintHash_Table * table,
						    size_t numKeys,
						    cl_mem keysBuffer,
						    cl_mem valuesOutputBuffer);
int intintIdentitySentinelPerfectCLHash_BufferInsert(intintHash_Table * table,
						     size_t numEntries,
						     cl_mem keysBuffer,
						     cl_mem valuesBuffer);
int intintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite(intintHash_Table
								* table,
								size_t
								numEntries,
								cl_mem
								keysBuffer,
								cl_mem
								valuesBuffer);
int intintLCGLinearOpenCompactHash_CreateFactory(intintHash_Factory * factory,
						 int hashIndex);
int intintLCGLinearOpenCompactHash_DestroyFactory(intintHash_Factory * factory,
						  int hashIndex);
intintHash_Table *intintLCGLinearOpenCompactHash_CreateTable(intintHash_Factory
							     * factory,
							     int hashIndex,
							     size_t keyRange,
							     size_t numEntries,
							     float loadFactor);
int intintLCGLinearOpenCompactHash_InitTable(intintHash_Table * table,
					     va_list args);
int intintLCGLinearOpenCompactHash_DestroyTable(intintHash_Table * table);
char *intintLCGLinearOpenCompactHash_Report(intintHash_Table * table);
int intintLCGLinearOpenCompactHash_EmptyTable(intintHash_Table * table);
int intintLCGLinearOpenCompactHash_Query(intintHash_Table * table,
					 size_t numKeys, int *keys,
					 int *valuesOutput);
int intintLCGLinearOpenCompactHash_QuerySingle(intintHash_Table * table,
					       int key, int *valueOutput);
int intintLCGLinearOpenCompactHash_Insert(intintHash_Table * table,
					  size_t numEntries, int *keys,
					  int *values);
int intintLCGLinearOpenCompactHash_InsertSingle(intintHash_Table * table,
						int key, int value);
int intintLCGLinearOpenCompactHash_InsertNoOverwrite(intintHash_Table * table,
						     size_t numEntries,
						     int *keys, int *values);
int intintLCGLinearOpenCompactHash_InsertSingleNoOverwrite(intintHash_Table *
							   table, int key,
							   int value);
int intintLCGLinearOpenCompactCLHash_CreateFactory(intintHash_Factory * factory,
						   int hashIndex);
int intintLCGLinearOpenCompactCLHash_DestroyFactory(intintHash_Factory *
						    factory, int hashIndex);
intintHash_Table
    *intintLCGLinearOpenCompactCLHash_CreateTable(intintHash_Factory * factory,
						  int hashIndex,
						  size_t keyRange,
						  size_t numEntries,
						  float loadFactor);
int intintLCGLinearOpenCompactCLHash_InitTable(intintHash_Table * table,
					       va_list args);
int intintLCGLinearOpenCompactCLHash_DestroyTable(intintHash_Table * table);
char *intintLCGLinearOpenCompactCLHash_Report(intintHash_Table * table);
int intintLCGLinearOpenCompactCLHash_EmptyTable(intintHash_Table * table);
int intintLCGLinearOpenCompactCLHash_Query(intintHash_Table * table,
					   size_t numKeys, int *keys,
					   int *valuesOutput);
int intintLCGLinearOpenCompactCLHash_QuerySingle(intintHash_Table * table,
						 int key, int *valueOutput);
int intintLCGLinearOpenCompactCLHash_Insert(intintHash_Table * table,
					    size_t numEntries, int *keys,
					    int *values);
int intintLCGLinearOpenCompactCLHash_InsertSingle(intintHash_Table * table,
						  int key, int value);
int intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(intintHash_Table * table,
						       size_t numEntries,
						       int *keys, int *values);
int intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(intintHash_Table *
							     table, int key,
							     int value);
int intintLCGLinearOpenCompactCLHash_BufferQuery(intintHash_Table * table,
						 size_t numKeys,
						 cl_mem keysBuffer,
						 cl_mem valuesOutputBuffer);
int intintLCGLinearOpenCompactCLHash_BufferInsert(intintHash_Table * table,
						  size_t numEntries,
						  cl_mem keysBuffer,
						  cl_mem valuesBuffer);
int intintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite(intintHash_Table *
							     table,
							     size_t numEntries,
							     cl_mem keysBuffer,
							     cl_mem
							     valuesBuffer);
int intintLCGQuadraticOpenCompactHash_CreateFactory(intintHash_Factory *
						    factory, int hashIndex);
int intintLCGQuadraticOpenCompactHash_DestroyFactory(intintHash_Factory *
						     factory, int hashIndex);
intintHash_Table
    *intintLCGQuadraticOpenCompactHash_CreateTable(intintHash_Factory * factory,
						   int hashIndex,
						   size_t keyRange,
						   size_t numEntries,
						   float loadFactor);
int intintLCGQuadraticOpenCompactHash_InitTable(intintHash_Table * table,
						va_list args);
int intintLCGQuadraticOpenCompactHash_DestroyTable(intintHash_Table * table);
char *intintLCGQuadraticOpenCompactHash_Report(intintHash_Table * table);
int intintLCGQuadraticOpenCompactHash_EmptyTable(intintHash_Table * table);
int intintLCGQuadraticOpenCompactHash_Query(intintHash_Table * table,
					    size_t numKeys, int *keys,
					    int *valuesOutput);
int intintLCGQuadraticOpenCompactHash_QuerySingle(intintHash_Table * table,
						  int key, int *valueOutput);
int intintLCGQuadraticOpenCompactHash_Insert(intintHash_Table * table,
					     size_t numEntries, int *keys,
					     int *values);
int intintLCGQuadraticOpenCompactHash_InsertSingle(intintHash_Table * table,
						   int key, int value);
int intintLCGQuadraticOpenCompactHash_InsertNoOverwrite(intintHash_Table *
							table,
							size_t numEntries,
							int *keys, int *values);
int intintLCGQuadraticOpenCompactHash_InsertSingleNoOverwrite(intintHash_Table *
							      table, int key,
							      int value);
int intintLCGQuadraticOpenCompactCLHash_CreateFactory(intintHash_Factory *
						      factory, int hashIndex);
int intintLCGQuadraticOpenCompactCLHash_DestroyFactory(intintHash_Factory *
						       factory, int hashIndex);
intintHash_Table
    *intintLCGQuadraticOpenCompactCLHash_CreateTable(intintHash_Factory *
						     factory, int hashIndex,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor);
int intintLCGQuadraticOpenCompactCLHash_InitTable(intintHash_Table * table,
						  va_list args);
int intintLCGQuadraticOpenCompactCLHash_DestroyTable(intintHash_Table * table);
char *intintLCGQuadraticOpenCompactCLHash_Report(intintHash_Table * table);
int intintLCGQuadraticOpenCompactCLHash_EmptyTable(intintHash_Table * table);
int intintLCGQuadraticOpenCompactCLHash_Query(intintHash_Table * table,
					      size_t numKeys, int *keys,
					      int *valuesOutput);
int intintLCGQuadraticOpenCompactCLHash_QuerySingle(intintHash_Table * table,
						    int key, int *valueOutput);
int intintLCGQuadraticOpenCompactCLHash_Insert(intintHash_Table * table,
					       size_t numEntries, int *keys,
					       int *values);
int intintLCGQuadraticOpenCompactCLHash_InsertSingle(intintHash_Table * table,
						     int key, int value);
int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(intintHash_Table *
							  table,
							  size_t numEntries,
							  int *keys,
							  int *values);
int intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(intintHash_Table
								* table,
								int key,
								int value);
int intintLCGQuadraticOpenCompactCLHash_BufferQuery(intintHash_Table * table,
						    size_t numKeys,
						    cl_mem keysBuffer,
						    cl_mem valuesOutputBuffer);
int intintLCGQuadraticOpenCompactCLHash_BufferInsert(intintHash_Table * table,
						     size_t numEntries,
						     cl_mem keysBuffer,
						     cl_mem valuesBuffer);
int intintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite(intintHash_Table
								* table,
								size_t
								numEntries,
								cl_mem
								keysBuffer,
								cl_mem
								valuesBuffer);

#ifdef __cplusplus
}
#endif				/* __cplusplus */

#endif				/* HASHFACTORY_H */
