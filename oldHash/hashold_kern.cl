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

/*
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Peter Ahrens            peter.ahrens@lanl.gov, ptrahrens@gmail.com
 *           Sara Hartse             sara@lanl.gov, sara.hartse@gmail.com
 *           Rebecka Tumblin         rtumblin@lanl.gov, rebeckatumblin@gmail.com
 */

/* neigh2d_kern.cl */

#ifndef  __HAVE_CL_DOUBLE_DEFINED__
#define __HAVE_CL_DOUBLE_DEFINED__

#ifdef HAVE_CL_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double  real;
#else
typedef float   real;
#endif

#endif

__constant ulong prime = 4294967291;
__constant uint hash_jump_prime = 41;

void write_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const uint  giX,
            const ulong hashkey,
   __global       int   *hash)
{
   uint hashloc;
   int icount = 0;
   uint jump;
   int old_key;
   int MaxTries = 1000;
#ifndef __APPLE_CC__
   switch (hash_method) {
   case -1:
#endif
      hash[hashkey] = giX;
#ifndef __APPLE_CC__
      break;
   case 0:
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (int icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc++;
         hashloc %= hash_table_size;
     
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }    

      if (icount < MaxTries) hash[2*hashloc+1] = giX; 
      break;
   case 1:
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (int icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc+=(icount*icount);
         hashloc %= hash_table_size;
     
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }    

      if (icount < MaxTries) hash[2*hashloc+1] = giX; 
      break;
   case 2:
      jump = 1+hashkey%hash_jump_prime;
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (int icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc += (icount*jump);
         hashloc %= hash_table_size;

         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }

      if (icount < MaxTries) hash[2*hashloc+1] = giX;
      break;
   }
#endif
}

int read_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const ulong hashkey,
   __global const int   *hash)
{
   int hashval = -1;
   uint hashloc;
   int icount = 0;
   uint jump;

#ifndef __APPLE_CC__
   switch (hash_method) {
   case -1:
#endif
      return(hash[hashkey]);
#ifndef __APPLE_CC__
      break;
   case 0:
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc++,hashloc %= hash_table_size);
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   case 1:
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc %= hash_table_size){
         icount++;
      }
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   case 2:
      jump = 1+hashkey%hash_jump_prime;
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc %= hash_table_size){
         icount++;
      }
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   }
#endif
   return(hashval);
}
