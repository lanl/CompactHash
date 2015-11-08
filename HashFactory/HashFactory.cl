
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
int intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(__global char
							       *tableData,
							       unsigned int
							       numEntries,
							       __global int
							       *keys,
							       __global int
							       *values);
int intintIdentitySentinelPerfectCLHash_InnerQuerySingle(__global char
							 *tableData, int key,
							 __global int
							 *valueOutput);
int intintIdentitySentinelPerfectCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput);
int intintIdentitySentinelPerfectCLHash_InnerInsertSingle(__global char
							  *tableData, int key,
							  int value);
int intintIdentitySentinelPerfectCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values);
int intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite(__global
								     char
								     *tableData,
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
int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(__global char
							  *tableData,
							  size_t numEntries,
							  __global int *keys,
							  __global int *values);
int intintLCGLinearOpenCompactCLHash_InnerQuerySingle(__global char *tableData,
						      int key,
						      __global int
						      *valueOutput);
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
int intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(__global char
							    *tableData,
							    unsigned int
							    numEntries,
							    __global int *keys,
							    __global int
							    *values);
int intintLCGLinearOpenCompactCLHash_InnerInsertSingle(__global char *tableData,
						       int key, int value);
int intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global char
								  *tableData,
								  int key,
								  int value);
int intintLCGLinearOpenCompactCLHash_InnerInsert(__global char *tableData,
						 unsigned int numEntries,
						 __global int *keys,
						 __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(__global char
							 *tableData, int key,
							 __global int
							 *valueOutput);
int intintLCGQuadraticOpenCompactCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(__global char
							  *tableData, int key,
							  int value);
int intintLCGQuadraticOpenCompactCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global
								     char
								     *tableData,
								     int key,
								     int value);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(__global char
							       *tableData,
							       unsigned int
							       numEntries,
							       __global int
							       *keys,
							       __global int
							       *values);
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
int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(__global char
							  *tableData,
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
#define IDENTITY_PERFECT_CL_HASH_ID /****************/ 16
#define IDENTITY_SENTINEL_PERFECT_CL_HASH_ID /*******/ 32
#define LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID /*********/ 64
#define LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID /******/ 128
//
#define HASH_BUCKET_STATUS_EMPTY /**/ -1
#define HASH_BUCKET_STATUS_FULL /***/ -2
#define HASH_BUCKET_STATUS_LOCK /***/ -3
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

typedef struct intintIdentityPerfectCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
} intintIdentityPerfectCLHash_TableData;
typedef struct intintIdentityPerfectCLHash_Bucket {
	int key;
	int value;
} intintIdentityPerfectCLHash_Bucket;
int intintIdentityPerfectCLHash_InnerQuerySingle(__global char *tableData,
						 int key,
						 __global int *valueOutput) {
	__global intintIdentityPerfectCLHash_Bucket *buckets =
	    (__global intintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((__global
					  intintIdentityPerfectCLHash_TableData
					  *) tableData)->compressFuncData, key);
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
int intintIdentityPerfectCLHash_InnerQuery(__global char *tableData,
					   unsigned int numKeys,
					   __global int *keys,
					   __global int *valuesOutput) {
	__global intintIdentityPerfectCLHash_Bucket *buckets =
	    (__global intintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];
	int key;
	__global int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    intintHash_CompressIdentity(((__global
						  intintIdentityPerfectCLHash_TableData
						  *) tableData)->
						compressFuncData, key);
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
int intintIdentityPerfectCLHash_InnerInsertSingle(__global char *tableData,
						  int key, int value) {
	__global intintIdentityPerfectCLHash_Bucket *buckets =
	    (__global intintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((__global
					  intintIdentityPerfectCLHash_TableData
					  *) tableData)->compressFuncData, key);
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
int intintIdentityPerfectCLHash_InnerInsert(__global char *tableData,
					    unsigned int numEntries,
					    __global int *keys,
					    __global int *values) {
	__global intintIdentityPerfectCLHash_Bucket *buckets =
	    (__global intintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((__global
						  intintIdentityPerfectCLHash_TableData
						  *) tableData)->
						compressFuncData, key);
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
int intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(__global char
							     *tableData,
							     int key,
							     int value) {
	__global intintIdentityPerfectCLHash_Bucket *buckets =
	    (__global intintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((__global
					  intintIdentityPerfectCLHash_TableData
					  *) tableData)->compressFuncData, key);
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
int intintIdentityPerfectCLHash_InnerInsertNoOverwrite(__global char *tableData,
						       unsigned int numEntries,
						       __global int *keys,
						       __global int *values) {
	__global intintIdentityPerfectCLHash_Bucket *buckets =
	    (__global intintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentityPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((__global
						  intintIdentityPerfectCLHash_TableData
						  *) tableData)->
						compressFuncData, key);
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
int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,
					    __global int *valueOutput) {
	return intintIdentityPerfectCLHash_InnerQuerySingle(tableData, key,
							    valueOutput);
}
int intintIdentityPerfectCLHash_Query(__global char *tableData, size_t numKeys,
				      __global int *keys,
				      __global int *valuesOutput) {
	return intintIdentityPerfectCLHash_InnerQuery(tableData, numKeys, keys,
						      valuesOutput);
}
int intintIdentityPerfectCLHash_InsertSingle(__global char *tableData, int key,
					     int value) {
	return intintIdentityPerfectCLHash_InnerInsertSingle(tableData, key,
							     value);
}
int intintIdentityPerfectCLHash_Insert(__global char *tableData,
				       size_t numEntries, __global int *keys,
				       __global int *values) {
	return intintIdentityPerfectCLHash_InnerInsert(tableData, numEntries,
						       keys, values);
}
int intintIdentityPerfectCLHash_InsertSingleNoOverwrite(__global char
							*tableData, int key,
							int value) {
	return
	    intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(tableData,
								     key,
								     value);
}
int intintIdentityPerfectCLHash_InsertNoOverwrite(__global char *tableData,
						  size_t numEntries,
						  __global int *keys,
						  __global int *values) {
	return intintIdentityPerfectCLHash_InnerInsertNoOverwrite(tableData,
								  numEntries,
								  keys, values);
}
__kernel void intintIdentityPerfectCLHash_RangeQuerySingle(__global char
							   *tableData,
							   unsigned int
							   numQueries,
							   __global int *keys,
							   __global int
							   *valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintIdentityPerfectCLHash_InnerQuerySingle(tableData, keys[i],
						     valuesOutput + i);
}
__kernel void intintIdentityPerfectCLHash_RangeQuery(__global char *tableData,
						     unsigned int numQueries,
						     unsigned int numKeys,
						     __global int *keys,
						     __global int
						     *valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintIdentityPerfectCLHash_InnerQuery(tableData, numKeys,
					       keys + (i * numKeys),
					       valuesOutput + (i * numKeys));
}
__kernel void intintIdentityPerfectCLHash_RangeInsertSingle(__global char
							    *tableData,
							    unsigned int
							    numInsertions,
							    __global int *keys,
							    __global int
							    *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentityPerfectCLHash_InnerInsertSingle(tableData, keys[i],
						      values[i]);
}
__kernel void intintIdentityPerfectCLHash_RangeInsert(__global char *tableData,
						      unsigned int
						      numInsertions,
						      unsigned int numEntries,
						      __global int *keys,
						      __global int *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentityPerfectCLHash_InnerInsert(tableData, numEntries,
						keys + (i * numEntries),
						values + (i * numEntries));
}
__kernel void intintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite(__global
								       char
								       *tableData,
								       unsigned
								       int
								       numInsertions,
								       __global
								       int
								       *keys,
								       __global
								       int
								       *values) 
{
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(tableData,
								 keys[i],
								 values[i]);
}
__kernel void intintIdentityPerfectCLHash_RangeInsertNoOverwrite(__global char
								 *tableData,
								 unsigned int
								 numInsertions,
								 unsigned int
								 numEntries,
								 __global int
								 *keys,
								 __global int
								 *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentityPerfectCLHash_InnerInsertNoOverwrite(tableData,
							   numEntries,
							   keys +
							   (i * numEntries),
							   values +
							   (i * numEntries));
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
int intintIdentitySentinelPerfectCLHash_InnerQuerySingle(__global char
							 *tableData, int key,
							 __global int
							 *valueOutput) {
	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((__global
					  intintIdentitySentinelPerfectCLHash_TableData
					  *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((__global intintIdentitySentinelPerfectCLHash_TableData *)
	     tableData)->emptyValue) {
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
int intintIdentitySentinelPerfectCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput) {
	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];
	int key;
	__global int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    intintHash_CompressIdentity(((__global
						  intintIdentitySentinelPerfectCLHash_TableData
						  *) tableData)->
						compressFuncData, key);
		if (buckets[index].value !=
		    ((__global intintIdentitySentinelPerfectCLHash_TableData *)
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
int intintIdentitySentinelPerfectCLHash_InnerInsertSingle(__global char
							  *tableData, int key,
							  int value) {
	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((__global
					  intintIdentitySentinelPerfectCLHash_TableData
					  *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((__global intintIdentitySentinelPerfectCLHash_TableData *)
	     tableData)->emptyValue) {
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
int intintIdentitySentinelPerfectCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values) {
	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((__global
						  intintIdentitySentinelPerfectCLHash_TableData
						  *) tableData)->
						compressFuncData, key);
		if (buckets[index].value !=
		    ((__global intintIdentitySentinelPerfectCLHash_TableData *)
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
int intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite(__global
								     char
								     *tableData,
								     int key,
								     int value) 
{
	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];
	int index;
	int exitCode;
	index =
	    intintHash_CompressIdentity(((__global
					  intintIdentitySentinelPerfectCLHash_TableData
					  *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((__global intintIdentitySentinelPerfectCLHash_TableData *)
	     tableData)->emptyValue) {
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
int intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(__global char
							       *tableData,
							       unsigned int
							       numEntries,
							       __global int
							       *keys,
							       __global int
							       *values) {
	__global intintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global intintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    intintHash_CompressIdentity(((__global
						  intintIdentitySentinelPerfectCLHash_TableData
						  *) tableData)->
						compressFuncData, key);
		if (buckets[index].value !=
		    ((__global intintIdentitySentinelPerfectCLHash_TableData *)
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
int intintIdentitySentinelPerfectCLHash_QuerySingle(__global char *tableData,
						    int key,
						    __global int *valueOutput) {
	return intintIdentitySentinelPerfectCLHash_InnerQuerySingle(tableData,
								    key,
								    valueOutput);
}
int intintIdentitySentinelPerfectCLHash_Query(__global char *tableData,
					      size_t numKeys,
					      __global int *keys,
					      __global int *valuesOutput) {
	return intintIdentitySentinelPerfectCLHash_InnerQuery(tableData,
							      numKeys, keys,
							      valuesOutput);
}
int intintIdentitySentinelPerfectCLHash_InsertSingle(__global char *tableData,
						     int key, int value) {
	return intintIdentitySentinelPerfectCLHash_InnerInsertSingle(tableData,
								     key,
								     value);
}
int intintIdentitySentinelPerfectCLHash_Insert(__global char *tableData,
					       size_t numEntries,
					       __global int *keys,
					       __global int *values) {
	return intintIdentitySentinelPerfectCLHash_InnerInsert(tableData,
							       numEntries, keys,
							       values);
}
int intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(__global char
								*tableData,
								int key,
								int value) {
	return
	    intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite
	    (tableData, key, value);
}
int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(__global char
							  *tableData,
							  size_t numEntries,
							  __global int *keys,
							  __global int *values) 
{
	return
	    intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite
	    (tableData, numEntries, keys, values);
}
__kernel void intintIdentitySentinelPerfectCLHash_RangeQuerySingle(__global char
								   *tableData,
								   unsigned int
								   numQueries,
								   __global int
								   *keys,
								   __global int
								   *valuesOutput) 
{
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintIdentitySentinelPerfectCLHash_InnerQuerySingle(tableData, keys[i],
							     valuesOutput + i);
}
__kernel void intintIdentitySentinelPerfectCLHash_RangeQuery(__global char
							     *tableData,
							     unsigned int
							     numQueries,
							     unsigned int
							     numKeys,
							     __global int *keys,
							     __global int
							     *valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintIdentitySentinelPerfectCLHash_InnerQuery(tableData, numKeys,
						       keys + (i * numKeys),
						       valuesOutput +
						       (i * numKeys));
}
__kernel void intintIdentitySentinelPerfectCLHash_RangeInsertSingle(__global
								    char
								    *tableData,
								    unsigned int
								    numInsertions,
								    __global int
								    *keys,
								    __global int
								    *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentitySentinelPerfectCLHash_InnerInsertSingle(tableData,
							      keys[i],
							      values[i]);
}
__kernel void intintIdentitySentinelPerfectCLHash_RangeInsert(__global char
							      *tableData,
							      unsigned int
							      numInsertions,
							      unsigned int
							      numEntries,
							      __global int
							      *keys,
							      __global int
							      *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentitySentinelPerfectCLHash_InnerInsert(tableData, numEntries,
							keys + (i * numEntries),
							values +
							(i * numEntries));
}
__kernel void
intintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite(__global char
								 *tableData,
								 unsigned int
								 numInsertions,
								 __global int
								 *keys,
								 __global int
								 *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite
	    (tableData, keys[i], values[i]);
}
__kernel void
intintIdentitySentinelPerfectCLHash_RangeInsertNoOverwrite(__global char
							   *tableData,
							   unsigned int
							   numInsertions,
							   unsigned int
							   numEntries,
							   __global int *keys,
							   __global int
							   *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(tableData,
								   numEntries,
								   keys +
								   (i *
								    numEntries),
								   values +
								   (i *
								    numEntries));
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
int intintLCGLinearOpenCompactCLHash_InnerQuerySingle(__global char *tableData,
						      int key,
						      __global int
						      *valueOutput) {
	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];
	int index;
	int exitCode;
	__global intintLCGLinearOpenCompactCLHash_TableData *mytableData =
	    (__global intintLCGLinearOpenCompactCLHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((__global intintLCGLinearOpenCompactCLHash_TableData *)
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
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int intintLCGLinearOpenCompactCLHash_InnerQuery(__global char *tableData,
						unsigned int numKeys,
						__global int *keys,
						__global int *valuesOutput) {
	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];
	int key;
	__global int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		__global intintLCGLinearOpenCompactCLHash_TableData *mytableData
		    =
		    (__global intintLCGLinearOpenCompactCLHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((__global
			       intintLCGLinearOpenCompactCLHash_TableData *)
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
int intintLCGLinearOpenCompactCLHash_InnerInsertSingle(__global char *tableData,
						       int key, int value) {
	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];
	int index;
	int exitCode;
	__global intintLCGLinearOpenCompactCLHash_TableData *mytableData =
	    (__global intintLCGLinearOpenCompactCLHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((__global intintLCGLinearOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
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
int intintLCGLinearOpenCompactCLHash_InnerInsert(__global char *tableData,
						 unsigned int numEntries,
						 __global int *keys,
						 __global int *values) {
	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global intintLCGLinearOpenCompactCLHash_TableData *mytableData
		    =
		    (__global intintLCGLinearOpenCompactCLHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((__global
			       intintLCGLinearOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
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
int intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global char
								  *tableData,
								  int key,
								  int value) {
	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];
	int index;
	int exitCode;
	__global intintLCGLinearOpenCompactCLHash_TableData *mytableData =
	    (__global intintLCGLinearOpenCompactCLHash_TableData *) tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((__global intintLCGLinearOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
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
int intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(__global char
							    *tableData,
							    unsigned int
							    numEntries,
							    __global int *keys,
							    __global int
							    *values) {
	__global intintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global intintLCGLinearOpenCompactCLHash_TableData *mytableData
		    =
		    (__global intintLCGLinearOpenCompactCLHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((__global
			       intintLCGLinearOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
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
int intintLCGLinearOpenCompactCLHash_QuerySingle(__global char *tableData,
						 int key,
						 __global int *valueOutput) {
	return intintLCGLinearOpenCompactCLHash_InnerQuerySingle(tableData, key,
								 valueOutput);
}
int intintLCGLinearOpenCompactCLHash_Query(__global char *tableData,
					   size_t numKeys, __global int *keys,
					   __global int *valuesOutput) {
	return intintLCGLinearOpenCompactCLHash_InnerQuery(tableData, numKeys,
							   keys, valuesOutput);
}
int intintLCGLinearOpenCompactCLHash_InsertSingle(__global char *tableData,
						  int key, int value) {
	return intintLCGLinearOpenCompactCLHash_InnerInsertSingle(tableData,
								  key, value);
}
int intintLCGLinearOpenCompactCLHash_Insert(__global char *tableData,
					    size_t numEntries,
					    __global int *keys,
					    __global int *values) {
	return intintLCGLinearOpenCompactCLHash_InnerInsert(tableData,
							    numEntries, keys,
							    values);
}
int intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(__global char
							     *tableData,
							     int key,
							     int value) {
	return
	    intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite
	    (tableData, key, value);
}
int intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(__global char *tableData,
						       size_t numEntries,
						       __global int *keys,
						       __global int *values) {
	return
	    intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(tableData,
								    numEntries,
								    keys,
								    values);
}
__kernel void intintLCGLinearOpenCompactCLHash_RangeQuerySingle(__global char
								*tableData,
								unsigned int
								numQueries,
								__global int
								*keys,
								__global int
								*valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintLCGLinearOpenCompactCLHash_InnerQuerySingle(tableData, keys[i],
							  valuesOutput + i);
}
__kernel void intintLCGLinearOpenCompactCLHash_RangeQuery(__global char
							  *tableData,
							  unsigned int
							  numQueries,
							  unsigned int numKeys,
							  __global int *keys,
							  __global int
							  *valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintLCGLinearOpenCompactCLHash_InnerQuery(tableData, numKeys,
						    keys + (i * numKeys),
						    valuesOutput +
						    (i * numKeys));
}
__kernel void intintLCGLinearOpenCompactCLHash_RangeInsertSingle(__global char
								 *tableData,
								 unsigned int
								 numInsertions,
								 __global int
								 *keys,
								 __global int
								 *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGLinearOpenCompactCLHash_InnerInsertSingle(tableData, keys[i],
							   values[i]);
}
__kernel void intintLCGLinearOpenCompactCLHash_RangeInsert(__global char
							   *tableData,
							   unsigned int
							   numInsertions,
							   unsigned int
							   numEntries,
							   __global int *keys,
							   __global int
							   *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGLinearOpenCompactCLHash_InnerInsert(tableData, numEntries,
						     keys + (i * numEntries),
						     values + (i * numEntries));
}
__kernel void
intintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite(__global char
							      *tableData,
							      unsigned int
							      numInsertions,
							      __global int
							      *keys,
							      __global int
							      *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(tableData,
								      keys[i],
								      values
								      [i]);
}
__kernel void intintLCGLinearOpenCompactCLHash_RangeInsertNoOverwrite(__global
								      char
								      *tableData,
								      unsigned
								      int
								      numInsertions,
								      unsigned
								      int
								      numEntries,
								      __global
								      int *keys,
								      __global
								      int
								      *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(tableData,
								numEntries,
								keys +
								(i *
								 numEntries),
								values +
								(i *
								 numEntries));
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
int intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(__global char
							 *tableData, int key,
							 __global int
							 *valueOutput) {
	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];
	int index;
	int exitCode;
	__global intintLCGQuadraticOpenCompactCLHash_TableData *mytableData =
	    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)
	    tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((__global intintLCGQuadraticOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((__global
			   intintLCGQuadraticOpenCompactCLHash_TableData *)
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
int intintLCGQuadraticOpenCompactCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput) {
	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];
	int key;
	__global int *valueOutput;
	int index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		__global intintLCGQuadraticOpenCompactCLHash_TableData
		    *mytableData =
		    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((__global
			       intintLCGQuadraticOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((__global
				   intintLCGQuadraticOpenCompactCLHash_TableData
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
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(__global char
							  *tableData, int key,
							  int value) {
	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];
	int index;
	int exitCode;
	__global intintLCGQuadraticOpenCompactCLHash_TableData *mytableData =
	    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)
	    tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((__global intintLCGQuadraticOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((__global
			   intintLCGQuadraticOpenCompactCLHash_TableData *)
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
int intintLCGQuadraticOpenCompactCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values) {
	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global intintLCGQuadraticOpenCompactCLHash_TableData
		    *mytableData =
		    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((__global
			       intintLCGQuadraticOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((__global
				   intintLCGQuadraticOpenCompactCLHash_TableData
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
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global
								     char
								     *tableData,
								     int key,
								     int value) 
{
	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];
	int index;
	int exitCode;
	__global intintLCGQuadraticOpenCompactCLHash_TableData *mytableData =
	    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)
	    tableData;
	intintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = intintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((__global intintLCGQuadraticOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((__global
			   intintLCGQuadraticOpenCompactCLHash_TableData *)
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
int intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(__global char
							       *tableData,
							       unsigned int
							       numEntries,
							       __global int
							       *keys,
							       __global int
							       *values) {
	__global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global intintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	int key;
	int index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global intintLCGQuadraticOpenCompactCLHash_TableData
		    *mytableData =
		    (__global intintLCGQuadraticOpenCompactCLHash_TableData *)
		    tableData;
		intintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c = intintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((__global
			       intintLCGQuadraticOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((__global
				   intintLCGQuadraticOpenCompactCLHash_TableData
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
int intintLCGQuadraticOpenCompactCLHash_QuerySingle(__global char *tableData,
						    int key,
						    __global int *valueOutput) {
	return intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(tableData,
								    key,
								    valueOutput);
}
int intintLCGQuadraticOpenCompactCLHash_Query(__global char *tableData,
					      size_t numKeys,
					      __global int *keys,
					      __global int *valuesOutput) {
	return intintLCGQuadraticOpenCompactCLHash_InnerQuery(tableData,
							      numKeys, keys,
							      valuesOutput);
}
int intintLCGQuadraticOpenCompactCLHash_InsertSingle(__global char *tableData,
						     int key, int value) {
	return intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(tableData,
								     key,
								     value);
}
int intintLCGQuadraticOpenCompactCLHash_Insert(__global char *tableData,
					       size_t numEntries,
					       __global int *keys,
					       __global int *values) {
	return intintLCGQuadraticOpenCompactCLHash_InnerInsert(tableData,
							       numEntries, keys,
							       values);
}
int intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(__global char
								*tableData,
								int key,
								int value) {
	return
	    intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite
	    (tableData, key, value);
}
int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(__global char
							  *tableData,
							  size_t numEntries,
							  __global int *keys,
							  __global int *values) 
{
	return
	    intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite
	    (tableData, numEntries, keys, values);
}
__kernel void intintLCGQuadraticOpenCompactCLHash_RangeQuerySingle(__global char
								   *tableData,
								   unsigned int
								   numQueries,
								   __global int
								   *keys,
								   __global int
								   *valuesOutput) 
{
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(tableData, keys[i],
							     valuesOutput + i);
}
__kernel void intintLCGQuadraticOpenCompactCLHash_RangeQuery(__global char
							     *tableData,
							     unsigned int
							     numQueries,
							     unsigned int
							     numKeys,
							     __global int *keys,
							     __global int
							     *valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	intintLCGQuadraticOpenCompactCLHash_InnerQuery(tableData, numKeys,
						       keys + (i * numKeys),
						       valuesOutput +
						       (i * numKeys));
}
__kernel void intintLCGQuadraticOpenCompactCLHash_RangeInsertSingle(__global
								    char
								    *tableData,
								    unsigned int
								    numInsertions,
								    __global int
								    *keys,
								    __global int
								    *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(tableData,
							      keys[i],
							      values[i]);
}
__kernel void intintLCGQuadraticOpenCompactCLHash_RangeInsert(__global char
							      *tableData,
							      unsigned int
							      numInsertions,
							      unsigned int
							      numEntries,
							      __global int
							      *keys,
							      __global int
							      *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGQuadraticOpenCompactCLHash_InnerInsert(tableData, numEntries,
							keys + (i * numEntries),
							values +
							(i * numEntries));
}
__kernel void
intintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite(__global char
								 *tableData,
								 unsigned int
								 numInsertions,
								 __global int
								 *keys,
								 __global int
								 *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite
	    (tableData, keys[i], values[i]);
}
__kernel void
intintLCGQuadraticOpenCompactCLHash_RangeInsertNoOverwrite(__global char
							   *tableData,
							   unsigned int
							   numInsertions,
							   unsigned int
							   numEntries,
							   __global int *keys,
							   __global int
							   *values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(tableData,
								   numEntries,
								   keys +
								   (i *
								    numEntries),
								   values +
								   (i *
								    numEntries));
}
__kernel void intintHash_RangeQuery(__global char *tableData,
				    unsigned int numQueries,
				    unsigned int numKeys, __global int *keys,
				    __global int *valuesOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_RangeQuery(tableData,
							      numQueries,
							      numKeys, keys,
							      valuesOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return intintIdentitySentinelPerfectCLHash_RangeQuery(tableData,
								      numQueries,
								      numKeys,
								      keys,
								      valuesOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGLinearOpenCompactCLHash_RangeQuery(tableData,
								   numQueries,
								   numKeys,
								   keys,
								   valuesOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGQuadraticOpenCompactCLHash_RangeQuery(tableData,
								      numQueries,
								      numKeys,
								      keys,
								      valuesOutput);
	}
}
__kernel void intintHash_RangeQuerySingle(__global char *tableData,
					  unsigned int numQueries,
					  __global int *keys,
					  __global int *valueOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_RangeQuerySingle(tableData,
								    numQueries,
								    keys,
								    valueOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_RangeQuerySingle
		    (tableData, numQueries, keys, valueOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGLinearOpenCompactCLHash_RangeQuerySingle(tableData,
								      numQueries,
								      keys,
								      valueOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_RangeQuerySingle
		    (tableData, numQueries, keys, valueOutput);
	}
}
__kernel void intintHash_RangeInsert(__global char *tableData,
				     unsigned int numInsertions,
				     unsigned int numEntries,
				     __global int *keys, __global int *values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_RangeInsert(tableData,
							       numInsertions,
							       numEntries, keys,
							       values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_RangeInsert(tableData,
								    numInsertions,
								    numEntries,
								    keys,
								    values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGLinearOpenCompactCLHash_RangeInsert(tableData,
								    numInsertions,
								    numEntries,
								    keys,
								    values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_RangeInsert(tableData,
								    numInsertions,
								    numEntries,
								    keys,
								    values);
	}
}
__kernel void intintHash_RangeInsertSingle(__global char *tableData,
					   unsigned int numInsertions,
					   __global int *keys,
					   __global int *values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_RangeInsertSingle(tableData,
								     numInsertions,
								     keys,
								     values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_RangeInsertSingle
		    (tableData, numInsertions, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGLinearOpenCompactCLHash_RangeInsertSingle
		    (tableData, numInsertions, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_RangeInsertSingle
		    (tableData, numInsertions, keys, values);
	}
}
__kernel void intintHash_RangeInsertNoOverwrite(__global char *tableData,
						unsigned int numInsertions,
						unsigned int numEntries,
						__global int *keys,
						__global int *values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    intintIdentityPerfectCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGLinearOpenCompactCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	}
}
__kernel void intintHash_RangeInsertSingleNoOverwrite(__global char *tableData,
						      unsigned int
						      numInsertions,
						      __global int *keys,
						      __global int *values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    intintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	}
}
int intintHash_Query(__global char *tableData, unsigned int numKeys,
		     __global int *keys, __global int *valuesOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_Query(tableData, numKeys,
							 keys, valuesOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return intintIdentitySentinelPerfectCLHash_Query(tableData,
								 numKeys, keys,
								 valuesOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGLinearOpenCompactCLHash_Query(tableData,
							      numKeys, keys,
							      valuesOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGQuadraticOpenCompactCLHash_Query(tableData,
								 numKeys, keys,
								 valuesOutput);
	}
	return HASH_EXIT_CODE_ERROR;
}
int intintHash_QuerySingle(__global char *tableData, int key,
			   __global int *valueOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_QuerySingle(tableData, key,
							       valueOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_QuerySingle(tableData,
								    key,
								    valueOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGLinearOpenCompactCLHash_QuerySingle(tableData,
								    key,
								    valueOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_QuerySingle(tableData,
								    key,
								    valueOutput);
	}
	return HASH_EXIT_CODE_ERROR;
}
int intintHash_Insert(__global char *tableData, unsigned int numEntries,
		      __global int *keys, __global int *values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_Insert(tableData, numEntries,
							  keys, values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return intintIdentitySentinelPerfectCLHash_Insert(tableData,
								  numEntries,
								  keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGLinearOpenCompactCLHash_Insert(tableData,
							       numEntries, keys,
							       values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGQuadraticOpenCompactCLHash_Insert(tableData,
								  numEntries,
								  keys, values);
	}
	return HASH_EXIT_CODE_ERROR;
}
int intintHash_InsertSingle(__global char *tableData, int key, int value) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_InsertSingle(tableData, key,
								value);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_InsertSingle(tableData,
								     key,
								     value);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return intintLCGLinearOpenCompactCLHash_InsertSingle(tableData,
								     key,
								     value);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_InsertSingle(tableData,
								     key,
								     value);
	}
	return HASH_EXIT_CODE_ERROR;
}
int intintHash_InsertNoOverwrite(__global char *tableData,
				 unsigned int numEntries, __global int *keys,
				 __global int *values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return intintIdentityPerfectCLHash_InsertNoOverwrite(tableData,
								     numEntries,
								     keys,
								     values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_InsertNoOverwrite
		    (tableData, numEntries, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGLinearOpenCompactCLHash_InsertNoOverwrite
		    (tableData, numEntries, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite
		    (tableData, numEntries, keys, values);
	}
	return HASH_EXIT_CODE_ERROR;
}
int intintHash_InsertSingleNoOverwrite(__global char *tableData, int key,
				       int value) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    intintIdentityPerfectCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	}
	return HASH_EXIT_CODE_ERROR;
}
