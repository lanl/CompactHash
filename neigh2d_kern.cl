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

struct neighbor2d {
  int left;
  int right;
  int bottom;
  int top;
};

__kernel void init_kern(const uint size,
               __global int       *temp) {
	const uint idx = get_global_id(0);

  if (idx >= size) return;

	temp[idx] = -1;
}

__kernel void hash_kern(const real  min_val,
                        const real  min_diff,
                        const uint  length,
               __global const real *arr,
               __global int        *temp){
	const uint idx = get_global_id(0);
	
  if(idx >= length) return;

  temp[(uint)((arr[idx]-min_val)/min_diff)] = idx;
}

#define hashval(j,i) hash[(j)*imaxsize+(i)]
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

__kernel void neighbors2d_hashwrite_kern(const uint isize,
                                         const uint mesh_size,
                                         const uint levmx,
                                __global const int *levtable,
                                __global const int *i,
                                __global const int *j,
                                __global const int *level,
                                __global int       *hash) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];

  for (   int jjj = jjmin; jjj < jjmax; jjj++) {
    for (int iii = iimin; iii < iimax; iii++) {
      hashval(jjj, iii) = giX;
    }
  }
}

__kernel void neighbors2d_hashwrite_opt_1_kern(const uint isize,
                                               const uint mesh_size,
                                               const uint levmx,
                                      __global const int *levtable,
                                      __global const int *i,
                                      __global const int *j,
                                      __global const int *level,
                                      __global int       *hash){
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];
  
  for (int iii = iimin; iii < iimax; iii++) {
    hashval(jjmin,   iii) = giX;
    hashval(jjmax-1, iii) = giX;
  }
  for (int jjj = jjmin+1; jjj < jjmax-1; jjj++) {
    hashval(jjj, iimin) = giX;
    hashval(jjj, iimax-1) = giX;
  }
}

__kernel void neighbors2d_hashwrite_opt_2_kern(const uint isize,
                                               const uint mesh_size,
                                               const uint levmx,
                                      __global const int *levtable,
                                      __global const int *i,
                                      __global const int *j,
                                      __global const int *level,
                                      __global int       *hash) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;
  int levmult = levtable[levmx - lev];

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];
  
  if(lev == levmx) {
    hashval(jj,ii) = giX;
    return;
  }

  jj *= levmult;
  ii *= levmult;
  hashval(jj,ii) = giX;

  ii += levmult/2;
  hashval(jj,ii) = giX;
  if(levmult > 2) {
    ii += levmult/2 - 1;
    hashval(jj,ii) = giX;
    ii -= levmult/2 - 1;
  }
  ii -= levmult/2;
  jj += levmult/2;
  hashval(jj,ii) = giX;
  ii += levmult - 1;
  hashval(jj,ii) = giX;

  if(levmult > 2) {
    ii -= levmult - 1;
    jj += levmult/2 - 1;
    hashval(jj,ii) = giX;
    ii += levmult/2;
    hashval(jj,ii) = giX;
  }
}

__kernel void neighbors2d_hashwrite_opt_3_kern(const uint isize,
                                               const uint mesh_size,
                                               const uint levmx,
                                      __global const int *levtable,
                                      __global const int *i,
                                      __global const int *j,
                                      __global const int *level,
                                      __global int       *hash) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;
  int levmult = levtable[levmx - lev];

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];
  
  jj *= levmult;
  ii *= levmult;
  hashval(jj,ii) = giX;
}

__kernel void neighbors2d_hasholdlibwrite_opt_3_kern(const int   isize,
                                                     const int   levmx,
                                                     const int   imaxsize,
                                                     const int   hash_method,
                                                     const ulong hash_table_size,
                                                     const ulong AA,
                                                     const ulong BB,
                                            __global const int  *levtable,
                                            __global const int  *level,
                                            __global const int  *i,
                                            __global const int  *j,
                                            __global int        *hash){
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];
  int levmult = levtable[levmx - lev];

  jj *= levmult;
  ii *= levmult;

  write_hash(hash_method, hash_table_size, AA, BB, giX, jj*imaxsize+ii, hash);
}

__kernel void neighbors2d_hashread_kern(const int          isize,
                                        const uint         mesh_size,
                                        const int          levmx,
                               __global const int         *levtable,
                               __global const int         *i,
                               __global const int         *j,
                               __global const int         *level,
                               __global const int         *hash,
                               __global struct neighbor2d *neigh2d) {
  const uint giX  = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];

  int ii = i[giX];
  int jj = j[giX];
  int lev = level[giX];
  int levmult = levtable[levmx-lev];

  int nlftval = hashval(      jj   *levmult               , MAX(  ii   *levmult-1, 0         ));
  int nrhtval = hashval(      jj   *levmult               , MIN( (ii+1)*levmult,   imaxsize-1));
  int nbotval = hashval(MAX(  jj   *levmult-1, 0)         ,       ii   *levmult               );
  int ntopval = hashval(MIN( (jj+1)*levmult,   jmaxsize-1),       ii   *levmult               );

  neigh2d[giX].left = nlftval;
  neigh2d[giX].right = nrhtval;
  neigh2d[giX].bottom = nbotval;
  neigh2d[giX].top = ntopval;
}

__kernel void neighbors2d_hashread_opt_3_kern(const int          isize,
                                              const uint         mesh_size,
                                              const int          levmx,
                                     __global const int         *levtable,
                                     __global const int         *i,
                                     __global const int         *j,
                                     __global const int         *level,
                                     __global const int         *hash,
                                     __global struct neighbor2d *neigh2d) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  int imax = mesh_size;
  int jmax = mesh_size; 

  int ii = i[giX];
  int jj = j[giX];
  int lev = level[giX];
  int levmult = levtable[levmx-lev];
 
  int iicur = ii*levmult;
  int iilft = MAX( (ii-1)*levmult, 0         );
  int iirht = MIN( (ii+1)*levmult, imaxsize-1);
  int jjcur = jj*levmult;
  int jjbot = MAX( (jj-1)*levmult, 0         );
  int jjtop = MIN( (jj+1)*levmult, jmaxsize-1);

  int nlftval = -1;
  int nrhtval = -1;
  int nbotval = -1;
  int ntopval = -1;

  // Taking care of boundary cells
  // Force each boundary cell to point to itself on its boundary direction
  if (iicur <    1  ) nlftval = giX;
  if (jjcur <    1  ) nbotval = giX;
  if ((ii + 1)*levmult >(imax)*levtable[levmx]-1) nrhtval = giX;
  if ((jj + 1)*levmult >(jmax)*levtable[levmx]-1) ntopval = giX;
  // Need to check for finer neighbor first
  // Right and top neighbor don't change for finer, so drop through to same size
  // Left and bottom need to be half of same size index for finer test
     
  if (lev != levmx) {
    int iilftfiner = iicur-(iicur-iilft)/2;
    int jjbotfiner = jjcur-(jjcur-jjbot)/2;
    if (nlftval < 0) nlftval = hashval(jjcur, iilftfiner);
    if (nbotval < 0) nbotval = hashval(jjbotfiner,iicur);
  }

  // same size neighbor
  if (nlftval < 0) nlftval = hashval(jjcur,iilft);
  if (nrhtval < 0) nrhtval = hashval(jjcur,iirht);
  if (nbotval < 0) nbotval = hashval(jjbot,iicur);
  if (ntopval < 0) ntopval = hashval(jjtop,iicur);
  // coarser neighbor
  if (lev != 0){
    if (nlftval < 0) {
      iilft -= iicur-iilft;
      int jjlft = (jj/2)*2*levmult;
      nlftval = hashval(jjlft,iilft);
    }
    if (nrhtval < 0) {
      int jjrht = (jj/2)*2*levmult;
      nrhtval = hashval(jjrht,iirht);
    }
    if (nbotval < 0) {
      jjbot -= jjcur-jjbot;
      int iibot = (ii/2)*2*levmult;
      nbotval = hashval(jjbot,iibot);
    }
    if (ntopval < 0) {
      int iitop = (ii/2)*2*levmult;
      ntopval = hashval(jjtop,iitop);
    }
  }

  neigh2d[giX].left   = nlftval;
  neigh2d[giX].right  = nrhtval;
  neigh2d[giX].bottom = nbotval;
  neigh2d[giX].top    = ntopval;
}

__kernel void neighbors2d_hasholdlibread_opt_3_kern(const int          isize,
                                                    const uint         mesh_size,
                                                    const int          levmx,
                                           __global const int         *levtable,
                                           __global const int         *i,
                                           __global const int         *j,
                                           __global const int         *level,
                                           __global const int         *hash,
                                           __global struct neighbor2d *neigh2d,
                                                    const int          hash_method,
                                                    const ulong        hash_table_size,
                                                    const ulong        AA,
                                                    const ulong        BB) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  int imax = mesh_size;
  int jmax = mesh_size; 

  int ii = i[giX];
  int jj = j[giX];
  int lev = level[giX];
  int levmult = levtable[levmx-lev];
 
  int iicur = ii*levmult;
  int iilft = MAX( (ii-1)*levmult, 0         );
  int iirht = MIN( (ii+1)*levmult, imaxsize-1);
  int jjcur = jj*levmult;
  int jjbot = MAX( (jj-1)*levmult, 0         );
  int jjtop = MIN( (jj+1)*levmult, jmaxsize-1);

  int nlftval = -1;
  int nrhtval = -1;
  int nbotval = -1;
  int ntopval = -1;

  // Taking care of boundary cells
  // Force each boundary cell to point to itself on its boundary direction
  if (iicur <    1  ) nlftval = giX;
  if (jjcur <    1  ) nbotval = giX;
  if ((ii + 1)*levmult >(imax)*levtable[levmx]-1) nrhtval = giX;
  if ((jj + 1)*levmult >(jmax)*levtable[levmx]-1) ntopval = giX;
  // Need to check for finer neighbor first
  // Right and top neighbor don't change for finer, so drop through to same size
  // Left and bottom need to be half of same size index for finer test
     
  if (lev != levmx) {
    int iilftfiner = iicur-(iicur-iilft)/2;
    int jjbotfiner = jjcur-(jjcur-jjbot)/2;
    if (nlftval < 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur *imaxsize+iilftfiner, hash); 
    if (nbotval < 0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*imaxsize+iicur,   hash);
  }

  // same size neighbor
  if (nlftval < 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*imaxsize+iilft, hash); 
  if (nrhtval < 0) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*imaxsize+iirht, hash);
  if (nbotval < 0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbot*imaxsize+iicur, hash);
  if (ntopval < 0) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*imaxsize+iicur, hash);
  // coarser neighbor
  if (lev != 0){
    if (nlftval < 0) {
      iilft -= iicur-iilft;
      int jjlft = (jj/2)*2*levmult;
      nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjlft*imaxsize+iilft, hash);
    }
    if (nrhtval < 0) {
      int jjrht = (jj/2)*2*levmult;
      nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjrht*imaxsize+iirht, hash);
    }
    if (nbotval < 0) {
      jjbot -= jjcur-jjbot;
      int iibot = (ii/2)*2*levmult;
      nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbot*imaxsize+iibot, hash);
    }
    if (ntopval < 0) {
      int iitop = (ii/2)*2*levmult;
      ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*imaxsize+iitop, hash);
    }
  }

  neigh2d[giX].left   = nlftval;
  neigh2d[giX].right  = nrhtval;
  neigh2d[giX].bottom = nbotval;
  neigh2d[giX].top    = ntopval;
}

__kernel void neighbors2d_hashlibwrite_kern(const uint isize,
                                            const uint mesh_size,
                                            const uint levmx,
                                   __global const int *levtable,
                                   __global const int *i,
                                   __global const int *j,
                                   __global const int *level,
                                   __global char      *hashTableData) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];

  for (   int jjj = jjmin; jjj < jjmax; jjj++) {
    for (int iii = iimin; iii < iimax; iii++) {
      intintHash_InsertSingle(hashTableData, jjj*imaxsize+iii, giX);
    }
  }
}

__kernel void neighbors2d_hashlibwrite_opt_1_kern(const uint isize,
                                                  const uint mesh_size,
                                                  const uint levmx,
                                         __global const int *levtable,
                                         __global const int *i,
                                         __global const int *j,
                                         __global const int *level,
                                         __global char      *hashTableData) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];

  /*
  if (lev == levmx) {
    intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
  }else{
    int width = levtable[levmx-lev] - 1;
    for(int x = 0; x < width; x++){
      intintHash_InsertSingle(hashTableData, (iimin + x    ) + (jjmin        ) * imaxsize, giX);
      intintHash_InsertSingle(hashTableData, (iimin + x + 1) + (jjmin + width) * imaxsize, giX);
      intintHash_InsertSingle(hashTableData, (iimin + width) + (jjmin + x    ) * imaxsize, giX);
      intintHash_InsertSingle(hashTableData, (iimin        ) + (jjmin + x + 1) * imaxsize, giX);
    }
  }
  return;
  */

  for (int iii = iimin; iii < iimax; iii++) {
    intintHash_InsertSingle(hashTableData, jjmin*imaxsize+iii, giX);
    intintHash_InsertSingle(hashTableData, (jjmax-1)*imaxsize+iii, giX);
  }
  for (int jjj = jjmin+1; jjj < jjmax-1; jjj++) {
    intintHash_InsertSingle(hashTableData, jjj*imaxsize+iimin, giX);
    intintHash_InsertSingle(hashTableData, jjj*imaxsize+iimax-1, giX);
  }
  return;
}

__kernel void neighbors2d_hashlibwrite_opt_2_kern(const uint isize,
                                                  const uint mesh_size,
                                                  const uint levmx,
                                         __global const int *levtable,
                                         __global const int *i,
                                         __global const int *j,
                                         __global const int *level,
                                         __global char      *hashTableData) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;
  int levmult = levtable[levmx - lev];

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];
  
  if(lev == levmx) {
    intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
    return;
  }

  jj *= levmult;
  ii *= levmult;
  intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);

  ii += levmult/2;
  intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
  if(levmult > 2) {
    ii += levmult/2 - 1;
    intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
    ii -= levmult/2 - 1;
  }
  ii -= levmult/2;
  jj += levmult/2;
  intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
  ii += levmult - 1;
  intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);

  if(levmult > 2) {
    ii -= levmult - 1;
    jj += levmult/2 - 1;
    intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
    ii += levmult/2;
    intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
  }
}

__kernel void neighbors2d_hashlibwrite_opt_3_kern(const uint isize,
                                                  const uint mesh_size,
                                                  const uint levmx,
                                         __global const int *levtable,
                                         __global const int *i,
                                         __global const int *j,
                                         __global const int *level,
                                         __global char      *hashTableData) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];

  int lev = level[giX];
  int ii = i[giX];
  int jj = j[giX];

  int levdiff = levmx - lev;
  int levmult = levtable[levmx - lev];

  int iimin =  ii   *levtable[levdiff];
  int iimax = (ii+1)*levtable[levdiff];
  int jjmin =  jj   *levtable[levdiff];
  int jjmax = (jj+1)*levtable[levdiff];
  
  jj *= levmult;
  ii *= levmult;
  intintHash_InsertSingle(hashTableData, jj * imaxsize + ii, giX);
}

__kernel void neighbors2d_hashlibread_opt_3_kern(const int          isize,
                                                 const uint         mesh_size,
                                                 const int          levmx,
                                        __global const int         *levtable,
                                        __global const int         *i,
                                        __global const int         *j,
                                        __global const int         *level,
                                        __global char              *hashTableData,
                                        __global struct neighbor2d *neigh2d) {
  const uint giX = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  int imax = mesh_size;
  int jmax = mesh_size; 

  int ii = i[giX];
  int jj = j[giX];
  int lev = level[giX];
  int levmult = levtable[levmx-lev];
 
  int iicur = ii*levmult;
  int iilft = MAX( (ii-1)*levmult, 0         );
  int iirht = MIN( (ii+1)*levmult, imaxsize-1);
  int jjcur = jj*levmult;
  int jjbot = MAX( (jj-1)*levmult, 0         );
  int jjtop = MIN( (jj+1)*levmult, jmaxsize-1);

  neigh2d[giX].left = -1;
  neigh2d[giX].right = -1;
  neigh2d[giX].bottom = -1;
  neigh2d[giX].top = -1;
  // Taking care of boundary cells
  // Force each boundary cell to point to itself on its boundary direction
  if (iicur <    1  ) neigh2d[giX].left = giX;
  if (jjcur <    1  ) neigh2d[giX].bottom = giX;
  if ((ii + 1)*levmult >(imax)*levtable[levmx]-1) neigh2d[giX].right = giX;
  if ((jj + 1)*levmult >(jmax)*levtable[levmx]-1) neigh2d[giX].top = giX;
  // Need to check for finer neighbor first
  // Right and top neighbor don't change for finer, so drop through to same size
  // Left and bottom need to be half of same size index for finer test
     
  if (lev != levmx) {
    int iilftfiner = iicur-(iicur-iilft)/2;
    int jjbotfiner = jjcur-(jjcur-jjbot)/2;
    if (neigh2d[giX].left == -1) intintHash_QuerySingle(hashTableData, iilftfiner + jjcur * imaxsize, &neigh2d[giX].left);
    if (neigh2d[giX].bottom == -1) intintHash_QuerySingle(hashTableData, iicur + jjbotfiner * imaxsize, &neigh2d[giX].bottom);
  }

  // same size neighbor
  if (neigh2d[giX].left == -1) intintHash_QuerySingle(hashTableData, iilft + jjcur * imaxsize, &neigh2d[giX].left);
  if (neigh2d[giX].right == -1) intintHash_QuerySingle(hashTableData, iirht + jjcur * imaxsize, &neigh2d[giX].right);
  if (neigh2d[giX].bottom == -1) intintHash_QuerySingle(hashTableData, iicur + jjbot * imaxsize, &neigh2d[giX].bottom);
  if (neigh2d[giX].top == -1) intintHash_QuerySingle(hashTableData, iicur + jjtop * imaxsize, &neigh2d[giX].top);
  // coarser neighbor
  if (lev != 0){
    if (neigh2d[giX].left == -1) {
      iilft -= iicur-iilft;
      int jjlft = (jj/2)*2*levmult;
      intintHash_QuerySingle(hashTableData, iilft + jjlft * imaxsize, &neigh2d[giX].left);
    }
    if (neigh2d[giX].right == -1) {
      int jjrht = (jj/2)*2*levmult;
      intintHash_QuerySingle(hashTableData, iirht + jjrht * imaxsize, &neigh2d[giX].right);
    }
    if (neigh2d[giX].bottom == -1) {
      jjbot -= jjcur-jjbot;
      int iibot = (ii/2)*2*levmult;
      intintHash_QuerySingle(hashTableData, iibot + jjbot * imaxsize, &neigh2d[giX].bottom);
    }
    if (neigh2d[giX].top == -1) {
      int iitop = (ii/2)*2*levmult;
      intintHash_QuerySingle(hashTableData, iitop + jjtop * imaxsize, &neigh2d[giX].top);
    }
  }
}

__kernel void neighbors2d_hashlibread_kern(const int          isize,
                                           const uint         mesh_size,
                                           const int          levmx,
                                  __global const int         *levtable,
                                  __global const int         *i,
                                  __global const int         *j,
                                  __global const int         *level,
                                  __global char              *hashTableData,
                                  __global struct neighbor2d *neigh2d) {
  const uint giX  = get_global_id(0);

  if (giX >= isize) return;

  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];

  int ii = i[giX];
  int jj = j[giX];
  int lev = level[giX];
  int levmult = levtable[levmx-lev];

  intintHash_QuerySingle(hashTableData, (MAX(  ii   *levmult-1, 0         )) + (      jj   *levmult               ) * imaxsize, &neigh2d[giX].left);
  intintHash_QuerySingle(hashTableData, (MIN( (ii+1)*levmult,   imaxsize-1)) + (      jj   *levmult               ) * imaxsize, &neigh2d[giX].right);
  intintHash_QuerySingle(hashTableData, (      ii   *levmult               ) + (MAX(  jj   *levmult-1, 0)         ) * imaxsize, &neigh2d[giX].bottom);
  intintHash_QuerySingle(hashTableData, (      ii   *levmult               ) + (MIN( (jj+1)*levmult,   jmaxsize-1)) * imaxsize, &neigh2d[giX].top);
}
