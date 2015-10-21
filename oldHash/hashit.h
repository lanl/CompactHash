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

#ifndef _HASH_H
#define _HASH_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef unsigned int uint;
typedef unsigned long ulong;

int *compact_hash_init(int ncells, uint isize, uint jsize, uint report_level);

#ifdef _OPENMP
   #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
      int *compact_hash_init_openmp(int ncells, uint isize, uint jsize, uint report_level);
   #else
      int *compact_hash_init_openmp(int ncells, uint isize, uint jsize, uint report_level, omp_lock_t **lock);
   #endif
#endif

int get_hash_method(void);
long long get_hash_table_size(void);

void write_hash_perfect(uint ic, ulong hashkey, int *hash);
void write_hash_linear(uint ic, ulong hashkey, int *hash);
void write_hash_linear_report_level_1(uint ic, ulong hashkey, int *hash);
void write_hash_linear_report_level_2(uint ic, ulong hashkey, int *hash);
void write_hash_linear_report_level_3(uint ic, ulong hashkey, int *hash);
void write_hash_quadratic(uint ic, ulong hashkey, int *hash);
void write_hash_quadratic_report_level_1(uint ic, ulong hashkey, int *hash);
void write_hash_quadratic_report_level_2(uint ic, ulong hashkey, int *hash);
void write_hash_quadratic_report_level_3(uint ic, ulong hashkey, int *hash);
void write_hash_primejump(uint ic, ulong hashkey, int *hash);
void write_hash_primejump_report_level_1(uint ic, ulong hashkey, int *hash);
void write_hash_primejump_report_level_2(uint ic, ulong hashkey, int *hash);
void write_hash_primejump_report_level_3(uint ic, ulong hashkey, int *hash);
extern void (*write_hash)(uint ic, ulong hashkey, int *hash); // declared in hash.c

#ifdef _OPENMP
   #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
      void write_hash_perfect_openmp(uint ic, ulong hashkey, int *hash);
      void write_hash_linear_openmp(uint ic, ulong hashkey, int *hash);
      void write_hash_linear_openmp_report_level_1(uint ic, ulong hashkey, int *hash);
      void write_hash_linear_openmp_report_level_2(uint ic, ulong hashkey, int *hash);
      void write_hash_linear_openmp_report_level_3(uint ic, ulong hashkey, int *hash);
      void write_hash_quadratic_openmp(uint ic, ulong hashkey, int *hash);
      void write_hash_quadratic_openmp_report_level_1(uint ic, ulong hashkey, int *hash);
      void write_hash_quadratic_openmp_report_level_2(uint ic, ulong hashkey, int *hash);
      void write_hash_quadratic_openmp_report_level_3(uint ic, ulong hashkey, int *hash);
      void write_hash_primejump_openmp(uint ic, ulong hashkey, int *hash);
      void write_hash_primejump_openmp_report_level_1(uint ic, ulong hashkey, int *hash);
      void write_hash_primejump_openmp_report_level_2(uint ic, ulong hashkey, int *hash);
      void write_hash_primejump_openmp_report_level_3(uint ic, ulong hashkey, int *hash);
      extern void (*write_hash_openmp)(uint ic, ulong hashkey, int *hash); // declared in hash.c
   #else
      void write_hash_perfect_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_linear_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_linear_openmp_report_level_1(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_linear_openmp_report_level_2(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_linear_openmp_report_level_3(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_quadratic_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_quadratic_openmp_report_level_1(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_quadratic_openmp_report_level_2(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_quadratic_openmp_report_level_3(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_primejump_openmp(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_primejump_openmp_report_level_1(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_primejump_openmp_report_level_2(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      void write_hash_primejump_openmp_report_level_3(uint ic, ulong hashkey, int *hash, omp_lock_t *lock);
      extern void (*write_hash_openmp)(uint ic, ulong hashkey, int *hash, omp_lock_t *lock); // declared in hash.c
   #endif
#endif

int read_hash_perfect(ulong hashkey, int *hash);
int read_hash_linear(ulong hashkey, int *hash);
int read_hash_linear_report_level_1(ulong hashkey, int *hash);
int read_hash_linear_report_level_2(ulong hashkey, int *hash);
int read_hash_linear_report_level_3(ulong hashkey, int *hash);
int read_hash_quadratic(ulong hashkey, int *hash);
int read_hash_quadratic_report_level_1(ulong hashkey, int *hash);
int read_hash_quadratic_report_level_2(ulong hashkey, int *hash);
int read_hash_quadratic_report_level_3(ulong hashkey, int *hash);
int read_hash_primejump(ulong hashkey, int *hash);
int read_hash_primejump_report_level_1(ulong hashkey, int *hash);
int read_hash_primejump_report_level_2(ulong hashkey, int *hash);
int read_hash_primejump_report_level_3(ulong hashkey, int *hash);
extern int (*read_hash)(ulong hashkey, int *hash); // declared in hash.c

void compact_hash_delete(int *hash);
#ifdef _OPENMP
   #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
      void compact_hash_delete_openmp(int *hash);
   #else
      void compact_hash_delete_openmp(int *hash, omp_lock_t *lock);
   #endif
#endif

void read_hash_collision_report(void);
void final_hash_collision_report(void);

#ifdef __cplusplus
}
#endif

#endif // _HASH_H
