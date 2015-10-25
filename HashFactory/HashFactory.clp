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

#define GLOBAL __global

int intintIdentityPerfectCLHash_InsertSingle(__global char *tableData,
                                                  int key, int value);
int intintIdentityPerfectCLHash_InnerInsertSingle(__global char *tableData,
						  int key, int value);
int intintHash_InsertSingle(__global char *tableData, int key, int value);
int intintIdentityPerfectCLHash_InnerQuery(__global char *tableData,
					   unsigned int numKeys,
					   __global int *keys,
					   __global int *valuesOutput);
int intintIdentityPerfectCLHash_InnerQuerySingle(__global char *tableData,
						 int key,
						 __global int *valueOutput);
int intintIdentityPerfectCLHash_InnerInsert(__global char *tableData,
					    unsigned int numEntries,
					    __global int *keys,
					    __global int *values);
int intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(__global char
							     *tableData,
							     int key,
							     int value);
int intintIdentityPerfectCLHash_InnerInsertNoOverwrite(__global char *tableData,
						       unsigned int numEntries,
						       __global int *keys,
						       __global int *values);
int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,
					    __global int *valueOutput);
int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,
					    __global int *valueOutput);
int intintIdentityPerfectCLHash_Query(__global char *tableData, size_t numKeys,
				      __global int *keys,
				      __global int *valuesOutput);
int intintIdentityPerfectCLHash_Insert(__global char *tableData,
				       size_t numEntries, __global int *keys,
				       __global int *values);
int intintIdentityPerfectCLHash_InsertSingleNoOverwrite(__global char
							*tableData, int key,
							int value);
int intintIdentityPerfectCLHash_InsertNoOverwrite(__global char *tableData,
						  size_t numEntries,
						  __global int *keys,
						  __global int *values);
int intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(__global char *tableData,
							       unsigned int numEntries,
							       __global int *keys,
							       __global int *values);
int intintIdentitySentinelPerfectCLHash_InnerQuerySingle(__global char *tableData, int key,
							 __global int *valueOutput);
int intintIdentitySentinelPerfectCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput);
int intintIdentitySentinelPerfectCLHash_InnerInsertSingle(__global char *tableData, int key,
							  int value);
int intintIdentitySentinelPerfectCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values);
int intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite(__global char *tableData,
								     int key,
								     int value);
int intintIdentitySentinelPerfectCLHash_QuerySingle(__global char *tableData,
						    int key,
						    __global int *valueOutput);
int intintIdentitySentinelPerfectCLHash_Query(__global char *tableData,
					      size_t numKeys,
					      __global int *keys,
					      __global int *valuesOutput);
int intintIdentitySentinelPerfectCLHash_InsertSingle(__global char *tableData,
						     int key, int value);
int intintIdentitySentinelPerfectCLHash_Insert(__global char *tableData,
					       size_t numEntries,
					       __global int *keys,
					       __global int *values);
int intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(__global char
								*tableData,
								int key,
								int value);
int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(__global char *tableData,
							  size_t numEntries,
							  __global int *keys,
							  __global int *values);
int intintLCGLinearOpenCompactCLHash_InnerQuerySingle(__global char *tableData,
						      int key,
						      __global int *valueOutput);
int intintLCGLinearOpenCompactCLHash_QuerySingle(__global char *tableData,
						 int key,
						 __global int *valueOutput);
int intintLCGLinearOpenCompactCLHash_Query(__global char *tableData,
					   size_t numKeys, __global int *keys,
					   __global int *valuesOutput);
int intintLCGLinearOpenCompactCLHash_InsertSingle(__global char *tableData,
						  int key, int value);
int intintLCGLinearOpenCompactCLHash_Insert(__global char *tableData,
					    size_t numEntries,
					    __global int *keys,
					    __global int *values);
int intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(__global char
							     *tableData,
							     int key,
							     int value);
int intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(__global char *tableData,
						       size_t numEntries,
						       __global int *keys,
						       __global int *values);
int intintLCGLinearOpenCompactCLHash_InnerQuery(__global char *tableData,
						unsigned int numKeys,
						__global int *keys,
						__global int *valuesOutput);
int intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(__global char *tableData,
							    unsigned int numEntries,
							    __global int *keys,
							    __global int *values);
int intintLCGLinearOpenCompactCLHash_InnerInsertSingle(__global char *tableData,
						       int key, int value);
int intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global char *tableData,
								  int key, int value);
int intintLCGLinearOpenCompactCLHash_InnerInsert(__global char *tableData,
						 unsigned int numEntries,
						 __global int *keys,
						 __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(__global char *tableData, int key,
							 __global int *valueOutput);
int intintLCGQuadraticOpenCompactCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(__global char *tableData, int key,
							  int value);
int intintLCGQuadraticOpenCompactCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global char *tableData,
								     int key, int value);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(__global char *tableData,
							       unsigned int numEntries,
							       __global int *keys,
							       __global int *values);
int intintLCGQuadraticOpenCompactCLHash_QuerySingle(__global char *tableData,
						    int key,
						    __global int *valueOutput);
int intintLCGQuadraticOpenCompactCLHash_Query(__global char *tableData,
					      size_t numKeys,
					      __global int *keys,
					      __global int *valuesOutput);
int intintLCGQuadraticOpenCompactCLHash_InsertSingle(__global char *tableData,
						     int key, int value);
int intintLCGQuadraticOpenCompactCLHash_Insert(__global char *tableData,
					       size_t numEntries,
					       __global int *keys,
					       __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(__global char
								*tableData,
								int key,
								int value);
int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(__global char *tableData,
							  size_t numEntries,
							  __global int *keys,
							  __global int *values);
int intintHash_Query(__global char *tableData, unsigned int numKeys,
		     __global int *keys, __global int *valuesOutput);
int intintHash_QuerySingle(__global char *tableData, int key,
			   __global int *valueOutput);
int intintHash_Insert(__global char *tableData, unsigned int numEntries,
		      __global int *keys, __global int *values);
int intintHash_InsertNoOverwrite(__global char *tableData,
				 unsigned int numEntries, __global int *keys,
				 __global int *values);
int intintHash_InsertSingleNoOverwrite(__global char *tableData, int key,
				       int value);
#define DELAY(X) X
DELAY(#define HASH_REPORT_NEVER /**/ 0)
DELAY(#define HASH_REPORT_CYCLE /**/ 1)
DELAY(#define HASH_REPORT_END /****/ 2)
//
DELAY(#define HASH_EXIT_CODE_NORMAL /****************/ -1)
DELAY(#define HASH_EXIT_CODE_ERROR /*****************/ -2)
DELAY(#define HASH_EXIT_CODE_OVERWRITE /*************/ -3)
DELAY(#define HASH_EXIT_CODE_KEY_DNE /***************/ -4)
DELAY(#define HASH_EXIT_CODE_CYCLE /*****************/ -5)
DELAY(#define HASH_EXIT_CODE_MAX_ENTRIES_EXCEEDED /**/ -6)
DELAY(#define HASH_EXIT_CODE_BUCKET_INDEX_OOB /******/ -7)
//
DELAY(#define HASH_SEARCH_CODE_MATCH /*****/ 0)
DELAY(#define HASH_SEARCH_CODE_MISMATCH /**/ 1)
DELAY(#define HASH_SEARCH_CODE_EMPTY /*****/ 2)
//
DELAY(#define IDENTITY_PERFECT_CL_HASH_ID /****************/ 16)
DELAY(#define IDENTITY_SENTINEL_PERFECT_CL_HASH_ID /*******/ 32)
DELAY(#define LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID /*********/ 64)
DELAY(#define LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID /******/ 128)
//
DELAY(#define HASH_BUCKET_STATUS_EMPTY /**/ -1)
DELAY(#define HASH_BUCKET_STATUS_FULL /***/ -2)
DELAY(#define HASH_BUCKET_STATUS_LOCK /***/ -3)
HASH_KERN_DEFINE(intint, int, int)
