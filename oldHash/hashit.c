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

#include "hashit.h"
//#include "timer/timer.h"

static ulong AA;
static ulong BB;
static ulong prime=4294967291;
static uint hashtablesize;
static uint hash_stride;
static uint hash_ncells;
static uint write_hash_collisions;
static uint read_hash_collisions;
static double write_hash_collisions_runsum = 0.0;
static double read_hash_collisions_runsum = 0.0;
static uint write_hash_collisions_count = 0;
static uint read_hash_collisions_count = 0;
static uint hash_report_level = 0;
static uint hash_queries;
static uint hash_method = 1;
static uint hash_jump_prime = 41;
static double hash_mult = 3.0;
static int do_compact_hash = 0;

float mem_opt_factor;

void compact_hash_delete(int *hash){
      free(hash);
}

int (*read_hash)(ulong, int *);
void (*write_hash)(uint, ulong, int *);

#ifdef _OPENMP
#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void (*write_hash_openmp)(uint, ulong, int *);
#else
void (*write_hash_openmp)(uint, ulong, int *, omp_lock_t * lock);
#endif
#endif

int get_hash_method(void) {
  return(hash_method);
}

long long get_hashtablesize(void) {
  return(hashtablesize);
}

int *compact_hash_init(int ncells, uint isize, uint jsize, uint report_level){

   hash_ncells = 0;
   write_hash_collisions = 0;
   read_hash_collisions = 0;
   hash_queries = 0;
   hash_report_level = report_level;
   hash_stride = isize;
   int *hash = NULL;

   uint compact_hash_size = (uint)((double)ncells*hash_mult);
   compact_hash_size = ncells;
   uint perfect_hash_size = (uint)(isize*jsize);
   float hash_mem_factor = 20.0;
   float hash_mem_ratio = (double)perfect_hash_size/(double)compact_hash_size;
   if (mem_opt_factor != 1.0) hash_mem_factor /= (mem_opt_factor*0.2); 
   do_compact_hash = (hash_mem_ratio < hash_mem_factor) ? 0 : 1;
   do_compact_hash = 1;
   if (hash_report_level >= 2) printf("DEBUG do_compact_hash %d hash_mem_ratio %f hash_mem_factor %f mem_opt_factor %f perfect_hash_size %u compact_hash_size %u\n",do_compact_hash,hash_mem_ratio,hash_mem_factor,mem_opt_factor,perfect_hash_size,compact_hash_size);

   if (do_compact_hash) {
      hashtablesize = compact_hash_size;
      AA = (ulong)(1.0+(double)(prime-1)*drand48());
      BB = (ulong)(0.0+(double)(prime-1)*drand48());
      if (AA > prime-1 || BB > prime-1) exit(0);
      if (hash_report_level > 1) printf("Factors AA %lu BB %lu\n",AA,BB);

      hash = (int *)malloc(2*hashtablesize*sizeof(int));
      for (uint ii = 0; ii<2*hashtablesize; ii+=2){
         hash[ii] = -1;
      }
   } else {
      hashtablesize = perfect_hash_size;

      hash = (int *)malloc(hashtablesize*sizeof(int));
      for (uint ii = 0; ii<hashtablesize; ii++){
         hash[ii] = -1;
      }
   }

   if (hash_report_level >= 2) {
      printf("Hash table size %u perfect hash table size %u memory savings %u by percentage %lf\n",
        hashtablesize,isize*jsize,isize*jsize-hashtablesize,
        (double)hashtablesize/(double)(isize*jsize));
   }

   if (! do_compact_hash) {
      read_hash  = read_hash_perfect;
      write_hash = write_hash_perfect;
   } else if (hash_method == 0){
      if (hash_report_level == 0){
         read_hash  = read_hash_linear;
         write_hash = write_hash_linear;
      } else if (hash_report_level == 1){
         read_hash  = read_hash_linear_report_level_1;
         write_hash = write_hash_linear_report_level_1;
      } else if (hash_report_level == 2){
         read_hash  = read_hash_linear_report_level_2;
         write_hash = write_hash_linear_report_level_2;
      } else if (hash_report_level == 3){
         read_hash  = read_hash_linear_report_level_3;
         write_hash = write_hash_linear_report_level_3;
      }
   } else if (hash_method == 1){
      if (hash_report_level == 0){
         read_hash  = read_hash_quadratic;
         write_hash = write_hash_quadratic;
      } else if (hash_report_level == 1){
         read_hash  = read_hash_quadratic_report_level_1;
         write_hash = write_hash_quadratic_report_level_1;
      } else if (hash_report_level == 2){
         read_hash  = read_hash_quadratic_report_level_2;
         write_hash = write_hash_quadratic_report_level_2;
      } else if (hash_report_level == 3){
         read_hash  = read_hash_quadratic_report_level_3;
         write_hash = write_hash_quadratic_report_level_3;
      }
   } else if (hash_method == 2){
      if (hash_report_level == 0){
         read_hash  = read_hash_primejump;
         write_hash = write_hash_primejump;
      } else if (hash_report_level == 1){
         read_hash  = read_hash_primejump_report_level_1;
         write_hash = write_hash_primejump_report_level_1;
      } else if (hash_report_level == 2){
         read_hash  = read_hash_primejump_report_level_2;
         write_hash = write_hash_primejump_report_level_2;
      } else if (hash_report_level == 3){
         read_hash  = read_hash_primejump_report_level_3;
         write_hash = write_hash_primejump_report_level_3;
      }
   }

   return(hash);
}

void write_hash_perfect(uint ic, ulong hashkey, int *hash){
   hash[hashkey] = ic;
}

void write_hash_linear(uint ic, ulong hashkey, int *hash){
   uint hashloc;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize);

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_linear_report_level_1(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
      write_hash_collisions++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_linear_report_level_2(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
      write_hash_collisions++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_linear_report_level_3(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
      int hashloctmp = hashloc+1;
      hashloctmp = hashloctmp%hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize) {
      icount++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic_report_level_1(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic_report_level_2(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic_report_level_3(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
      int hashloctmp = hashloc+icount*icount;
      hashloctmp = hashloctmp%hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize) {
      icount++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump_report_level_1(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump_report_level_2(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump_report_level_3(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;

   hash_ncells++;
   hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
      int hashloctmp = hashloc+1;
      hashloctmp = hashloctmp%hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

#ifdef _OPENMP

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_perfect_openmp(uint ic, ulong hashkey, int *hash){
#else
void write_hash_perfect_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   hash[hashkey] = ic;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_linear_openmp(uint ic, ulong hashkey, int *hash){
#else
void write_hash_linear_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc++;
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      hashloc++;
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_linear_openmp_report_level_1(uint ic, ulong hashkey, int *hash){
#else
void write_hash_linear_openmp_report_level_1(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc++;
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
      icount++;
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      hashloc++;
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
      icount++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_linear_openmp_report_level_2(uint ic, ulong hashkey, int *hash){
#else
void write_hash_linear_openmp_report_level_2(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc++;
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
      icount++;
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      hashloc++;
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
      icount++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_linear_openmp_report_level_3(uint ic, ulong hashkey, int *hash){
#else
void write_hash_linear_openmp_report_level_3(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;;
   printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc++;
      hashloc %= hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
      icount++;
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      hashloc++;
      hashloc = hashloc%hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
      omp_set_lock(&(lock[hashloc]));
      icount++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_quadratic_openmp(uint ic, ulong hashkey, int *hash){
#else
void write_hash_quadratic_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*icount);
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*icount);
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_quadratic_openmp_report_level_1(uint ic, ulong hashkey, int *hash){
#else
void write_hash_quadratic_openmp_report_level_1(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*icount);
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*icount);
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_quadratic_openmp_report_level_2(uint ic, ulong hashkey, int *hash){
#else
void write_hash_quadratic_openmp_report_level_2(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*icount);
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*icount);
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_quadratic_openmp_report_level_3(uint ic, ulong hashkey, int *hash){
#else
void write_hash_quadratic_openmp_report_level_3(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*icount);
      hashloc %= hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*icount);
      hashloc = hashloc%hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_primejump_openmp(uint ic, ulong hashkey, int *hash){
#else
void write_hash_primejump_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint jump = 1+hashkey%hash_jump_prime;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*jump);
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*jump);
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_primejump_openmp_report_level_1(uint ic, ulong hashkey, int *hash){
#else
void write_hash_primejump_openmp_report_level_1(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint jump = 1+hashkey%hash_jump_prime;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*jump);
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*jump);
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_primejump_openmp_report_level_2(uint ic, ulong hashkey, int *hash){
#else
void write_hash_primejump_openmp_report_level_2(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint jump = 1+hashkey%hash_jump_prime;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*jump);
      hashloc %= hashtablesize;

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*jump);
      hashloc = hashloc%hashtablesize;
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
void write_hash_primejump_openmp_report_level_3(uint ic, ulong hashkey, int *hash){
#else
void write_hash_primejump_openmp_report_level_3(uint ic, ulong hashkey, int *hash, omp_lock_t *lock){
#endif
   int icount = 0;
   uint jump = 1+hashkey%hash_jump_prime;
   uint hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
   int MaxTries = 1000;

   int old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   //printf("old_key is %d\n",old_key);

   for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
      hashloc+=(icount*jump);
      hashloc %= hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);

      old_key = __sync_val_compare_and_swap(&hash[2*hashloc], -1, hashkey); 
   }

   if (icount < MaxTries) hash[2*hashloc+1] = ic;
#else
   omp_set_lock(&(lock[hashloc]));
   while (hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey){
      omp_unset_lock(&(lock[hashloc]));
      icount++;
      hashloc+=(icount*jump);
      hashloc = hashloc%hashtablesize;
      printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
      omp_set_lock(&(lock[hashloc]));
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
   omp_unset_lock(&(lock[hashloc]));
#endif

#pragma omp atomic
   write_hash_collisions += icount;;
#pragma omp atomic
   hash_ncells++;
}
#endif

int read_hash_perfect(ulong hashkey, int *hash){
   return(hash[hashkey]);
}

int read_hash_linear(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize);

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_linear_report_level_1(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
      icount++;
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_linear_report_level_2(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
      icount++;
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_linear_report_level_3(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
      icount++;
      uint hashloctmp = hashloc+1;
      hashloctmp = hashloctmp%hashtablesize;
      printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic_report_level_1(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic_report_level_2(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic_report_level_3(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
      uint hashloctmp = hashloc+1;
      hashloctmp = hashloctmp%hashtablesize;
      printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump_report_level_1(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;
   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump_report_level_2(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump_report_level_3(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;

   hash_queries++;
   hashloc = (hashkey*AA+BB)%prime%hashtablesize;
   printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
      uint hashloctmp = hashloc+1;
      hashloctmp = hashloctmp%hashtablesize;
      printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_dev_hash(int hash_method, ulong hashtablesize, ulong AA, ulong BB, ulong hashkey, int *hash){
   //int hash_report_level = 3;
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;
   if (hash_method == -1) {
      return(hash[hashkey]);
   }
   if (hash_method == 0) {
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_queries++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
            uint hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      }
   } else if (hash_method == 1) {
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_queries++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
            uint hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      }
   } else if (hash_method == 2) {
      uint jump = 1+hashkey%hash_jump_prime;
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_queries++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
            uint hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      }
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

void write_hash_collision_report(void){
   if (! do_compact_hash) return;
   if (hash_report_level == 1) {
      write_hash_collisions_runsum += (double)write_hash_collisions/(double)hash_ncells;
      write_hash_collisions_count++;
   } else if (hash_report_level >= 2) {
      printf("Write hash collision report -- collisions per cell %lf, collisions %d cells %d\n",(double)write_hash_collisions/(double)hash_ncells,write_hash_collisions,hash_ncells);
   }
}

void read_hash_collision_report(void){
   if (! do_compact_hash) return;
   if (hash_report_level == 1) {
      read_hash_collisions_runsum += (double)read_hash_collisions/(double)hash_queries;
      read_hash_collisions_count++;
   } else if (hash_report_level >= 2) {
      printf("Read hash collision report -- collisions per cell %lf, collisions %d cells %d\n",(double)read_hash_collisions/(double)hash_queries,read_hash_collisions,hash_queries);
      hash_queries = 0;
      read_hash_collisions = 0;
   }
}

void final_hash_collision_report(void){
   if (hash_report_level >= 1 && read_hash_collisions_count > 0) {
      printf("Final hash collision report -- write/read collisions per cell %lf/%lf\n",write_hash_collisions_runsum/(double)write_hash_collisions_count,read_hash_collisions_runsum/(double)read_hash_collisions_count);
   }
}

