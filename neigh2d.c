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
//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "kdtree/KDTree2d.h"
#include "gpu.h"
#include "HashFactory/HashFactory.h"
#include "simplehash/simplehash.h"
#include "timer/timer.h"

#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#endif

#ifdef HAVE_CL_DOUBLE
typedef double real;
#ifdef HAVE_OPENCL
typedef cl_double cl_real;
typedef cl_double4 cl_real4;
#endif
#define SQRT(x) sqrt(x)
#define ONE 1.0
#define TWO 2.0
#else
typedef float real;
#ifdef HAVE_OPENCL
typedef cl_float cl_real;
typedef cl_float4 cl_real4;
#endif
#define SQRT(x) sqrtf(x)
#define ONE 1.0f
#define TWO 2.0f
#endif

#define SQ(x) (( (x)*(x) ))

#define HASH_TYPE HASH_PERFECT_HASHES
#define HASH_LOAD_FACTOR 0.3333333
int TILE_SIZE = 128;
int TILE_SIZE_INIT = 128;

#ifndef DETAILED_TIMING
#define DETAILED_TIMING 0
#endif

#ifndef WRITE_MEM_USAGE
#define WRITE_MEM_USAGE 0
#endif

#define CHECK 1
#define LONG_RUNS 1
#define PRINT_GOLD_ON_FAILED_TEST 1
#define PRINT_GOLD_ON_RUN 0
#define PRINT_PROBLEM_SET_ON_RUN 0

#ifndef MIN
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))
#define SWAP_PTR(p1,p2,p3) ((p3=p1), (p1=p2), (p2=p3))
#endif

int powerOfFour(int n) {
  int result = 1;
  int i;
  for(i = 0; i < n; i++) {
    result *= 4;
  }
  return result;
}

void swap_double(double** a, double** b) {
  double* c = *a;
  *a = *b;
  *b = c;
}

void swap_real(real** a, real** b) {
  real* c = *a;
  *a = *b;
  *b = c;
}

void swap_int(int** a, int** b) {
  int* c = *a;
  *a = *b;
  *b = c;
}

struct neighbor2d {
   int left;
   int right;
   int bottom;
   int top;
};

int device_type = VENDOR_UNKNOWN;

#define END_CHARACTER    -1

enum hash_type
{
  BRUTE            =  0,
  KDTREE           =  1,
  HASHCPU          =  2,
  HASHLIBCPU       =  3,
  HASHCPU_1        =  4,
  HASHLIBCPU_1     =  5,
  HASHCPU_2        =  6,
  HASHLIBCPU_2     =  7,
  HASHCPU_3        =  8,
  HASHLIBCPU_3     =  9,
  HASHOLDLIBCPU_3  = 10,
  HASHGPU          = 11,
  HASHGPU_1        = 12,
  HASHGPU_2        = 13,
  HASHGPU_3        = 14,
  HASHLIBGPU       = 15,
  HASHLIBGPU_1     = 16,
  HASHLIBGPU_2     = 17,
  HASHLIBGPU_3     = 18,
  HASHOLDLIBGPU_3  = 19
};

#ifdef HAVE_OPENCL
cl_context context;
cl_command_queue queue;
cl_program program;
cl_kernel neighbors2d_hashwrite_kern, neighbors2d_hashwrite_opt_1_kern, neighbors2d_hashwrite_opt_2_kern, neighbors2d_hashwrite_opt_3_kern, neighbors2d_hasholdlibwrite_opt_3_kern, neighbors2d_hashread_kern, neighbors2d_hashread_opt_3_kern, neighbors2d_hasholdlibread_opt_3_kern, neighbors2d_hashlibwrite_kern, neighbors2d_hashlibread_kern, neighbors2d_hashlibwrite_opt_1_kern, neighbors2d_hashlibwrite_opt_2_kern, neighbors2d_hashlibwrite_opt_3_kern, neighbors2d_hashlibread_opt_3_kern; 
#endif
intintHash_Factory *factory;
intintHash_Factory *CLFactory;

void neighbors2d( uint mesh_size, int levmx, float threshold, int *options, int haveGPU);
struct neighbor2d *neighbors2d_bruteforce( uint ncells, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_kdtree( uint ncells, int mesh_size, double *x, double *y, int *level );
struct neighbor2d *neighbors2d_hashcpu( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hashlibcpu( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hashcpu_opt_1( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hashlibcpu_opt_1( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hashcpu_opt_2( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hashlibcpu_opt_2( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hashcpu_opt_3( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hashlibcpu_opt_3( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
struct neighbor2d *neighbors2d_hasholdlibcpu_opt_3( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level );
cl_mem neighbors2d_hashgpu( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable);
#ifdef HAVE_OPENCL
cl_mem neighbors2d_hashgpu_opt_1( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable);
cl_mem neighbors2d_hashgpu_opt_2( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable );
cl_mem neighbors2d_hashgpu_opt_3( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable );
cl_mem neighbors2d_hasholdlibgpu_opt_3( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable );
cl_mem neighbors2d_hashlibgpu( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable);
cl_mem neighbors2d_hashlibgpu_opt_1( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable);
cl_mem neighbors2d_hashlibgpu_opt_2( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable);
cl_mem neighbors2d_hashlibgpu_opt_3( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable);
cl_mem neighbors2d_hashlibgpu_opt_4( uint ncells, int mesh_size, int levmx, cl_mem i, cl_mem j, cl_mem level, cl_mem levtable);
#endif

int adaptiveMeshConstructorWij(const int n, const int l, int** level_ptr, double** x_ptr, double** y_ptr, int **i_ptr, int **j_ptr, float threshold, int target_ncells);
void genmatrixfree(void **var);
void **genmatrix(int jnum, int inum, size_t elsize);
FILE * fmem; 
#ifdef HAVE_OPENCL
void cl_error_check(int error, char *file, int line){
  if (error != CL_SUCCESS) printf("Error is %d in file %s at line %d\n",error, file, line);
}
#endif
   
#ifndef bool
typedef int bool;
#endif

static bool randomize = false;
int main (int argc, const char * argv[])
{
  CLHash_Init(argv[0]);
  char header[1024];

  int nopts = 20;

  int   inputkey[]   = {BRUTE,        KDTREE,           HASHCPU,             HASHLIBCPU,     HASHCPU_1,    HASHLIBCPU_1,     HASHCPU_2,    HASHLIBCPU_2,     HASHCPU_3,    HASHLIBCPU_3,     HASHOLDLIBCPU_3,     HASHGPU,        HASHGPU_1,    HASHGPU_2,        HASHGPU_3,    HASHLIBGPU,     HASHLIBGPU_1,     HASHLIBGPU_2,     HASHLIBGPU_3,     HASHOLDLIBGPU_3};
  char *inputopt[]   = {"br",         "kd",             "hc",                "hlc",          "hc1",        "hlc1",           "hc2",        "hlc2",           "hc3",        "hlc3",           "holc3",             "hg",           "hg1",        "hg2",            "hg3",        "hlg",          "hlg1",           "hlg2",           "hlg3",           "holg3"};
  char *descriptor[] = {"Brute",      "kdtree",         "Hash CPU",          "Hash Lib CPU", "Hash CPU 1", "Hash Lib CPU 1", "Hash CPU 2", "Hash Lib CPU 2", "Hash CPU 3", "Hash Lib CPU 3", "Hash OldLib CPU 3", "Hash GPU",     "Hash GPU 1", "Hash GPU 2",     "Hash GPU 3", "Hash Lib GPU", "Hash Lib GPU 1", "Hash Lib GPU 2", "Hash Lib GPU 3", "Hash OldLib GPU 3"};

  char add_string[40];

  int *options;
  float threshold = 1.0;   

  uint levmx_min = 1;
  uint levmx_max = 6;
  uint n_min     = 256;
  uint n_max     = 256;

  if(WRITE_MEM_USAGE) fmem = fopen("memory.out", "w");
  if (argc == 1) {
    printf("Use -o to select different tests. -h for more information\n");
    return 0;
  }else{
    for (int i = 0; i < argc; i++) {
      if (strcmp(argv[i], "-h") == 0) {
        printf("Use '-t' to set threshold and '-o' to add tests: br kd hc hlc hc1 hlc1 hc2 hlc2 hc3 hlc3 holc3 hg hg1 hg2 hg3 hlg hlg1 hlg2 hlg3 holg3\n");
        printf("Use -r for randomization of cell data to see performance without cache\n");
        printf("Use -n for min mesh size and -N for max mesh size [ default is 256 min to 256 max ]\n");
        printf("Use -l for min levmx     and -L for max levmx     [ default is 1 min to 5 max ]\n");
        return 0;
      }
      if (strcmp(argv[i], "-r") == 0) {
        randomize = true; 
      } 
      if (strcmp(argv[i], "-t") == 0) {
        i++;
        threshold =(atof(argv[i])); 
      } 
      if (strcmp(argv[i], "-l") == 0) {
        i++;
        levmx_min =(int)(atoi(argv[i])); 
      } 
      if (strcmp(argv[i], "-L") == 0) {
        i++;
        levmx_max =(int)(atoi(argv[i])); 
        levmx_max++;
      } 
      if (strcmp(argv[i], "-n") == 0) {
        i++;
        n_min =(int)(atoi(argv[i])); 
      } 
      if (strcmp(argv[i], "-N") == 0) {
        i++;
        n_max =(int)(atoi(argv[i])); 
      }
      if (strcmp(argv[i], "-o") == 0) {
        i++;
        int option_size =  0; 
        while ( i + option_size < argc && argv[i+option_size][0] != '-') {
          option_size++;
        }
        if(option_size == 0){
          printf("Please add tests: br, kd, hc, hlc, hc1, hlc1, hc2, hlc2, hc3, hlc3, holc3, hg, hg1, hg2, hg3, hg4\n");
          return 0;
        }
        sprintf(header,"Size, \tlevmx: \tncells,  \tcompression,");
        if(WRITE_MEM_USAGE) fprintf(fmem,"Size,   \tncells,   ");
        fflush(stdout);
        options = (int*)malloc((option_size+2)*sizeof(int));
        int j = 0;
        int option;
        while (i < argc && argv[i][0] != '-' ) {
          option = -1;
          for (int iopt = 0; iopt < nopts; iopt++){
            if (strcmp(argv[i], inputopt[iopt]) == 0) { 
              option = inputkey[iopt];
              sprintf(add_string," \t%s,",descriptor[iopt]);
              strcat(header,add_string);
              if(WRITE_MEM_USAGE) fprintf(fmem,"\t%s,",descriptor[iopt]);
              break;
            }
          }
          if (option == -1){
            printf("Invalid Option: %s\n",argv[i]);
            exit(1);
          }
          options[j] = option;
          j++;
          i++;
        }
        options[j] = END_CHARACTER;
        //strcat(header,"\n");
        if(WRITE_MEM_USAGE) fprintf(fmem,"\n");
      }
    }
  }

  int haveGPU = 0;
  for(int z = 0; options[z] != END_CHARACTER; z++) {
    if (options[z] == HASHGPU || options[z] == HASHGPU_1 ||options[z] == HASHGPU_2 || options[z] == HASHGPU_3 || options[z] == HASHOLDLIBGPU_3 || options[z] == HASHLIBGPU || options[z] == HASHLIBGPU_1 || options[z] == HASHLIBGPU_2 || options[z] == HASHLIBGPU_3){ 
      haveGPU = 1;
      break;
    } 
  }

  cl_int error;

  int emptyNeighborValue = -1;
  factory = intintHash_CreateFactory(HASH_ALL_C_HASHES, &emptyNeighborValue, 0, NULL, NULL);
  if(haveGPU){
#ifdef HAVE_OPENCL
    char *bothsources = (char*)malloc(strlen(get_hash_kernel_source_string()) + strlen(Hash_GetKernelSourceString()) + 1);
    strcpy(bothsources, get_hash_kernel_source_string());
    strcat(bothsources, Hash_GetKernelSourceString());
    GPUInit(&context, &queue, &device_type, &program, "neigh2d_kern.cl", bothsources);
    if (device_type == MIC) {
       TILE_SIZE      = 180;
       TILE_SIZE_INIT = 240;
    }
    hash_lib_init(context);
    uint lws = TILE_SIZE;
    CLFactory = intintHash_CreateFactory(HASH_ALL_CL_HASHES, &emptyNeighborValue, lws, &context, &queue);

    neighbors2d_hashwrite_kern = clCreateKernel(program, "neighbors2d_hashwrite_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashwrite_opt_1_kern = clCreateKernel(program, "neighbors2d_hashwrite_opt_1_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashwrite_opt_2_kern = clCreateKernel(program, "neighbors2d_hashwrite_opt_2_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashwrite_opt_3_kern = clCreateKernel(program, "neighbors2d_hashwrite_opt_3_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashread_kern = clCreateKernel(program, "neighbors2d_hashread_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashread_opt_3_kern = clCreateKernel(program, "neighbors2d_hashread_opt_3_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashlibwrite_kern = clCreateKernel(program, "neighbors2d_hashlibwrite_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashlibwrite_opt_1_kern = clCreateKernel(program, "neighbors2d_hashlibwrite_opt_1_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashlibwrite_opt_2_kern = clCreateKernel(program, "neighbors2d_hashlibwrite_opt_2_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashlibwrite_opt_3_kern = clCreateKernel(program, "neighbors2d_hashlibwrite_opt_3_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashlibread_kern = clCreateKernel(program, "neighbors2d_hashlibread_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hashlibread_opt_3_kern = clCreateKernel(program, "neighbors2d_hashlibread_opt_3_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hasholdlibwrite_opt_3_kern = clCreateKernel(program, "neighbors2d_hasholdlibwrite_opt_3_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    neighbors2d_hasholdlibread_opt_3_kern = clCreateKernel(program, "neighbors2d_hasholdlibread_opt_3_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
#else
    printf("Error -- asking for GPU test and there is no GPU\n");
    exit(-1);
#endif
  }
  printf("\n    2D Neighbors Performance Results\n\n");
  printf("%s\n",header);

  //printf("\nThreshold for refinement is %f\n",threshold); 
  //printf("\nWork group size is %d and init workgroup size is %d\n",TILE_SIZE,TILE_SIZE_INIT); 
  for (uint levmx = levmx_min; levmx < levmx_max; levmx++ ){
    //printf("\nMax levels is %d\n",levmx); 
    if(WRITE_MEM_USAGE) fprintf(fmem,"\nMax levels is %d\n",levmx); 
    for( uint i = n_min; i <= n_max; i*=2 ) {
      //if (levmx > 3 && i > 512) continue;
      printf("%d, \t%d: ", i,levmx);
      if(WRITE_MEM_USAGE) fprintf(fmem,"%d,     ", i);
      neighbors2d(i, levmx, threshold, options, haveGPU);
      printf("\n");
      if(WRITE_MEM_USAGE) fprintf(fmem,"\n");
    }
  }
}

void neigh2dTester(int ncells, int *level, struct neighbor2d *neigh2d_gold, struct neighbor2d *neigh2d_test){
  for(int ic = 0; ic < ncells; ic++) {
    if(neigh2d_gold[ic].left   != neigh2d_test[ic].left   ||
       neigh2d_gold[ic].right  != neigh2d_test[ic].right  ||
       neigh2d_gold[ic].bottom != neigh2d_test[ic].bottom ||
       neigh2d_gold[ic].top    != neigh2d_test[ic].top    ) {
      if(PRINT_GOLD_ON_FAILED_TEST){
		    printf("Gold Cell %d: left = %3d \t right = %3d \t bottom = %3d \t top = %3d \t (level = %3d) \n",
               ic, neigh2d_gold[ic].left, neigh2d_gold[ic].right, neigh2d_gold[ic].bottom, neigh2d_gold[ic].top, level[ic]);      
      }
		  printf("Test Cell %d: left = %3d \t right = %3d \t bottom = %3d \t top = %3d \t (level = %3d) \n",
		         ic, neigh2d_test[ic].left, neigh2d_test[ic].right, neigh2d_test[ic].bottom, neigh2d_test[ic].top, level[ic]);      
		  printf("\n");
		}
	}
  free(neigh2d_test);
}

#ifdef HAVE_OPENCL
void neigh2dCLTester(int ncells, int *level, struct neighbor2d *neigh2d_gold, cl_mem neigh2d_buffer){
  struct neighbor2d *neigh2d_test = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  cl_int error = clEnqueueReadBuffer(queue, neigh2d_buffer, CL_TRUE, 0, ncells*sizeof(cl_uint4), neigh2d_test, 0, NULL, NULL);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
  clReleaseMemObject(neigh2d_buffer);
}
#endif

void neighbors2d( uint mesh_size, int levmx, float threshold, int *options, int haveGPU) 
{
  struct neighbor2d *neigh2d_gold, *neigh2d_test;

  int* level = NULL;
  double* x    = NULL;
  double* y    = NULL;
  int* i     = NULL;
  int* j     = NULL;

  int target_ncells = 0;
  int ncells = adaptiveMeshConstructorWij(mesh_size, levmx, &level, &x, &y, &i, &j,threshold, target_ncells);
  printf("\t%8d,", ncells);
  if(WRITE_MEM_USAGE) fprintf(fmem,"\t%8d,", ncells);
  
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  
  size_t jmaxsize = mesh_size*levtable[levmx];
  size_t imaxsize = mesh_size*levtable[levmx];
  float compression = (imaxsize*jmaxsize)/((double)ncells);
  printf("\t%.2f,        ", compression);
  if(PRINT_PROBLEM_SET_ON_RUN){
    for(int ic = 0; ic < ncells; ic++) {
      printf("Gold Cell %d: x = %f \t y = %f \t i = %3d \t j = %3d \t (level = %3d) \n",
             ic, x[ic], y[ic], i[ic], j[ic], level[ic]);      
    }
  }
  
  struct timespec timer;
  double time;
  cpu_timer_start(&timer);
  switch(options[0]){
    case BRUTE: neigh2d_gold = neighbors2d_bruteforce(ncells, i, j, level); break;
    default: neigh2d_gold = neighbors2d_kdtree(ncells, mesh_size, x, y, level); break;
  }
  time = cpu_timer_stop(timer);

  if(PRINT_GOLD_ON_RUN){
    for(int ic = 0; ic < ncells; ic++) {
      printf("Cell %d: left = %3d \t right = %3d \t bottom = %3d \t top = %3d \t (level = %3d) \n",
        ic, neigh2d_gold[ic].left, neigh2d_gold[ic].right, neigh2d_gold[ic].bottom, neigh2d_gold[ic].top, level[ic]);      
    }
  }

#ifdef HAVE_OPENCL
  cl_mem i_buffer;
  cl_mem j_buffer;
  cl_mem level_buffer;
  cl_int error;
  cl_mem levtable_buffer;
  cl_mem neigh2d_buffer;
#endif

  if (haveGPU) {
#ifdef HAVE_OPENCL
    error = 0;
    i_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(int), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    j_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(int), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    level_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(int), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clEnqueueWriteBuffer(queue, i_buffer, CL_TRUE, 0, ncells*sizeof(int), i, 0, NULL, NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clEnqueueWriteBuffer(queue, j_buffer, CL_TRUE, 0, ncells*sizeof(int), j, 0, NULL, NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clEnqueueWriteBuffer(queue, level_buffer, CL_TRUE, 0, ncells*sizeof(int), level, 0, NULL, NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

    int *levtable = (int *)malloc((levmx+1)*sizeof(int));
    for (int lev=0; lev<levmx+1; lev++){
      levtable[lev] = (int)pow(2,lev);
    }
    levtable_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, (levmx+1)*sizeof(int), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clEnqueueWriteBuffer(queue, levtable_buffer, CL_TRUE, 0, (levmx+1)*sizeof(int), levtable, 0, NULL, NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
#endif
  }

  for(int z = 0; options[z] != END_CHARACTER; z++) {
    switch(options[z]) {
      case BRUTE:
        if(z != 0){
          cpu_timer_start(&timer);
          neigh2d_test = neighbors2d_bruteforce(ncells, i, j, level);
          time = cpu_timer_stop(timer);
        }
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        if(z != 0){
          neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        }
        break;

      case KDTREE:
        if(z != 0){
          cpu_timer_start(&timer);
          neigh2d_test = neighbors2d_kdtree(ncells, mesh_size, x, y, level);
          time = cpu_timer_stop(timer);
        }
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        if(z != 0){
          neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        }
        break;

      case HASHCPU:
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashcpu(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

      case HASHLIBCPU:    
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashlibcpu(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

      case HASHCPU_1: 
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashcpu_opt_1(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

      case HASHLIBCPU_1:
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashlibcpu_opt_1(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

      case HASHCPU_2:
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashcpu_opt_2(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

      case HASHLIBCPU_2:
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashlibcpu_opt_2(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

      case HASHCPU_3:
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashcpu_opt_3(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;
              
      case HASHLIBCPU_3:
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hashlibcpu_opt_3(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

      case HASHOLDLIBCPU_3:
        cpu_timer_start(&timer);
        neigh2d_test = neighbors2d_hasholdlibcpu_opt_3(ncells, mesh_size, levmx, i, j, level);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dTester(ncells, level, neigh2d_gold, neigh2d_test);
        break;

#ifdef HAVE_OPENCL
      case HASHGPU:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashgpu(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHGPU_1:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashgpu_opt_1(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHGPU_2:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashgpu_opt_2(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHGPU_3:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashgpu_opt_3(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHLIBGPU:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashlibgpu(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHLIBGPU_1:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashlibgpu_opt_1(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHLIBGPU_2:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashlibgpu_opt_2(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHLIBGPU_3:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hashlibgpu_opt_3(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;

      case HASHOLDLIBGPU_3:
        cpu_timer_start(&timer);
        neigh2d_buffer = neighbors2d_hasholdlibgpu_opt_3(ncells, mesh_size, levmx, i_buffer, j_buffer, level_buffer, levtable_buffer);
        time = cpu_timer_stop(timer);
        printf("\t%.6lf,", time);
        fflush(stdout);
        if(WRITE_MEM_USAGE) fprintf(fmem, "\t%8s,", "x");
        neigh2dCLTester(ncells, level, neigh2d_gold, neigh2d_buffer);
        break;
#endif
    }
  }
  if(haveGPU){
#ifdef HAVE_OPENCL
    clReleaseMemObject(i_buffer);
    clReleaseMemObject(j_buffer);
    clReleaseMemObject(level_buffer);
    clReleaseMemObject(levtable_buffer);
#endif
  }
  free(neigh2d_gold);
  free(levtable);
  free(level);
  free(x);
  free(y);
  free(i);
  free(j);
}


/**
 * neighbors2d_bruteforce calculates neighbors through a brute force technique
 */
struct neighbor2d *neighbors2d_bruteforce( uint ncells, int *i, int *j, int *level )
{
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));

  for (uint index1 = 0; index1 < ncells; index1++) {
    int lev = level[index1];
    int ii = i[index1];
    int jj = j[index1];

    int left  = index1;
    for (uint index2 = 0; index2 < ncells; index2++) {
      if (abs(level[index2] - lev) > 1) continue;
      if ((level[index2] == lev && i[index2] ==  ii-1    && j[index2] == jj  ) ||
          (level[index2] <  lev && i[index2] == (ii-1)/2 && j[index2] == jj/2) ||
          (level[index2] >  lev && i[index2] ==  ii*2-1  && j[index2] == jj*2) ) {
        left = index2;
        break;
      }
    }
    int right = index1;
    for (uint index2 = 0; index2 < ncells; index2++) {
      if (abs(level[index2] - lev) > 1) continue;
      if ((level[index2] == lev && i[index2] ==  ii+1    && j[index2] == jj  ) ||
          (level[index2] <  lev && i[index2] == (ii+1)/2 && j[index2] == jj/2) ||
          (level[index2] >  lev && i[index2] == (ii+1)*2 && j[index2] == jj*2) ) {
        right = index2;
        break;
      }
    }
    int bottom = index1;
    for (uint index2 = 0; index2 < ncells; index2++) {
      if (abs(level[index2] - lev) > 1) continue;
      if ((level[index2] == lev && i[index2] == ii   && j[index2] ==  jj-1)    ||
          (level[index2] <  lev && i[index2] == ii/2 && j[index2] == (jj-1)/2) ||
          (level[index2] >  lev && i[index2] == ii*2 && j[index2] ==  jj*2-1)  ) {
        bottom = index2;
        break;
      }
    }
    int top = index1;
    for (uint index2 = 0; index2 < ncells; index2++) {
      if (abs(level[index2] - lev) > 1) continue;
      if ((level[index2] == lev && i[index2] == ii   && j[index2] ==  jj+1   ) ||
          (level[index2] < lev && i[index2] == ii/2 && j[index2] == (jj+1)/2) ||
          (level[index2] > lev && i[index2] == ii*2 && j[index2] == (jj+1)*2) ) {
        top = index2;
        break;
      }
    }
    neigh2d[index1].left   = left;
    neigh2d[index1].right  = right;
    neigh2d[index1].bottom = bottom;
    neigh2d[index1].top    = top;
  }

  return(neigh2d);
}

/**
 * neighbors2d_kdtree calculates neighbors through a brute force technique
 */
struct neighbor2d *neighbors2d_kdtree( uint ncells, int mesh_size, double *x, double *y, int *level ) 
{
  TKDTree2d tree;

  KDTree_Initialize2d(&tree);

  TBounds2d box;
  for(uint ic = 0; ic < ncells; ic++) {
    double lev_power = pow(2,(double)level[ic]);
    box.min.x = x[ic]-1.0*0.5/lev_power/mesh_size;
    box.max.x = x[ic]+1.0*0.5/lev_power/mesh_size;
    box.min.y = y[ic]-1.0*0.5/lev_power/mesh_size;
    box.max.y = y[ic]+1.0*0.5/lev_power/mesh_size;
    //printf("Adding cell %d : xmin %lf xmax %lf ymin %lf ymax %lf\n",ic,box.min.x,box.max.x,box.min.y,box.max.y);
    KDTree_AddElement2d(&tree, &box);
  }

  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));

  int index_list[20];
  int num;
  for (uint ic = 0; ic < ncells; ic++) {
    neigh2d[ic].left = ic;
    neigh2d[ic].right = ic;
    neigh2d[ic].bottom = ic;
    neigh2d[ic].top = ic;
    double lev_power = pow(2,(double)level[ic]);
    box.min.x = x[ic]-1.1*1.0*0.5/lev_power/mesh_size;
    box.max.x = x[ic]-1.1*1.0*0.5/lev_power/mesh_size;
    box.min.y = y[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    box.max.y = y[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    KDTree_QueryBoxIntersect2d(&tree, &num, &(index_list[0]), &box);
    if (num == 1) neigh2d[ic].left = index_list[0];

    box.min.x = x[ic]+1.1*1.0*0.5/lev_power/mesh_size;
    box.max.x = x[ic]+1.1*1.0*0.5/lev_power/mesh_size;
    box.min.y = y[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    box.max.y = y[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    KDTree_QueryBoxIntersect2d(&tree, &num, &(index_list[0]), &box);
    if (num == 1) neigh2d[ic].right = index_list[0];

    box.min.x = x[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    box.max.x = x[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    box.min.y = y[ic]-1.1*1.0*0.5/lev_power/mesh_size;
    box.max.y = y[ic]-1.1*1.0*0.5/lev_power/mesh_size;
    KDTree_QueryBoxIntersect2d(&tree, &num, &(index_list[0]), &box);
    if (num == 1) neigh2d[ic].bottom = index_list[0];

    box.min.x = x[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    box.max.x = x[ic]-0.5*1.0*0.5/lev_power/mesh_size;
    box.min.y = y[ic]+1.1*1.0*0.5/lev_power/mesh_size;
    box.max.y = y[ic]+1.1*1.0*0.5/lev_power/mesh_size;
    KDTree_QueryBoxIntersect2d(&tree, &num, &(index_list[0]), &box);
    if (num == 1) neigh2d[ic].top = index_list[0];

  }

  KDTree_Destroy2d(&tree);

  return(neigh2d);
}

/**
 * neighbors2d_hashcpu calculates neighbors with a hard-coded perfect hash and no memory optimizations
 */
struct neighbor2d *neighbors2d_hashcpu( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  size_t jmaxsize = mesh_size*levtable[levmx];
  size_t imaxsize = mesh_size*levtable[levmx];

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  int **hash = (int **)genmatrix(jmaxsize, imaxsize, sizeof(int));

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
//for (int jj = 0; jj<jmaxsize; jj++){
//  for (int ii = 0; ii<imaxsize; ii++){
//    hash[jj][ii]=-1;
//  }
//}

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = jmaxsize * imaxsize * sizeof(int); 
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  for(int ic=0; ic<ncells; ic++){
    int lev = level[ic];
    if (lev == levmx) {
      hash[j[ic]][i[ic]] = ic;
    } else {
      for (int jj = j[ic]*levtable[levmx-lev]; jj < (j[ic]+1)*levtable[levmx-lev]; jj++) {
        for (int ii=i[ic]*levtable[levmx-lev]; ii<(i[ic]+1)*levtable[levmx-lev]; ii++) {
          hash[jj][ii] = ic;
        }
      }
    }
  }

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    neigh2d[ic].left   = hash[      jj   *levmult               ][MAX(  ii   *levmult-1, 0         )];
    neigh2d[ic].right  = hash[      jj   *levmult               ][MIN( (ii+1)*levmult,   imaxsize-1)];
    neigh2d[ic].bottom = hash[MAX(  jj   *levmult-1, 0)         ][      ii   *levmult               ];
    neigh2d[ic].top    = hash[MIN( (jj+1)*levmult,   jmaxsize-1)][      ii   *levmult               ];
  }

  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  free(levtable);
  genmatrixfree((void**)hash);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}

/**
 * neighbors2d_hashlibcpu calculates neighbors with the library hash and no memory optimizations
 */
struct neighbor2d *neighbors2d_hashlibcpu( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  size_t jmaxsize = mesh_size*levtable[levmx];
  size_t imaxsize = mesh_size*levtable[levmx];
  int numKeys = imaxsize + (jmaxsize - 1) * imaxsize;

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(factory, HASH_PERFECT_HASHES, imaxsize + (jmaxsize - 1) * imaxsize, numKeys, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_EmptyTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = numKeys * sizeof(int); 
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

   //Write Cells
  for(int ic=0; ic<ncells; ic++){
    int lev = level[ic];
    if (lev == levmx) {
      intintHash_InsertSingleNoOverwrite(hashTable, i[ic] + j[ic] * imaxsize, ic);
    } else {
      for (int jj = j[ic]*levtable[levmx-lev]; jj < (j[ic]+1)*levtable[levmx-lev]; jj++) {
        for (int ii=i[ic]*levtable[levmx-lev]; ii<(i[ic]+1)*levtable[levmx-lev]; ii++) {
          intintHash_InsertSingleNoOverwrite(hashTable, ii + jj * imaxsize, ic);
        }
      }
    }
  }

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    intintHash_QuerySingle(hashTable, (MAX(  ii   *levmult-1, 0         )) + (      jj   *levmult               ) * imaxsize, &neigh2d[ic].left);
    intintHash_QuerySingle(hashTable, (MIN( (ii+1)*levmult,   imaxsize-1)) + (      jj   *levmult               ) * imaxsize, &neigh2d[ic].right);
    intintHash_QuerySingle(hashTable, (      ii   *levmult               ) + (MAX(  jj   *levmult-1, 0)         ) * imaxsize, &neigh2d[ic].bottom);
    intintHash_QuerySingle(hashTable, (      ii   *levmult               ) + (MIN( (jj+1)*levmult,   jmaxsize-1)) * imaxsize, &neigh2d[ic].top);
  }

  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}

/**
 * neighbors2d_hashcpu_opt_1 calculates neighbors with a hard-coded perfect hash and the boundary cell memory optimization
 */
struct neighbor2d *neighbors2d_hashcpu_opt_1(uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int jmaxsize = mesh_size*levtable[levmx];
  int imaxsize = mesh_size*levtable[levmx];

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);
   
  //Create Hash Table
  int **hash = (int **)genmatrix(jmaxsize, imaxsize, sizeof(int));

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);
   
  //Empty Hash Table
//for (int jj = 0; jj<jmaxsize; jj++){
//  for (int ii = 0; ii<imaxsize; ii++){
//    hash[jj][ii]=-1;
//  }
//}

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = jmaxsize * imaxsize * sizeof(int); 
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  for(int ic=0; ic<ncells; ic++){ 
    int lev = level[ic];
    int ii = i[ic];
    int jj = j[ic];
    int levmult = levtable[levmx - lev];
    int cellnumber = ic; //+noffset;
 
    int iimin = ii    *levmult;//-iminsize;
    int iimax = (ii+1)*levmult;//-iminsize;
    int jjmin = jj    *levmult;//-jminsize; 
    int jjmax = (jj+1)*levmult;//-jminsize; 

    for (int iii = iimin; iii < iimax; iii++) {
      hash[jjmin][iii]   = cellnumber;
      hash[jjmax-1][iii] = cellnumber;
    }
    for (int jjj = jjmin+1; jjj < jjmax-1; jjj++) {
      hash[jjj][iimin]   = cellnumber;
      hash[jjj][iimax-1] = cellnumber;
    }
  }

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);
  
  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii  = i[ic];
    int jj  = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    neigh2d[ic].left   = hash[      jj   *levmult               ][MAX(  ii   *levmult-1, 0         )];
    neigh2d[ic].right  = hash[      jj   *levmult               ][MIN( (ii+1)*levmult,   imaxsize-1)];
    neigh2d[ic].bottom = hash[MAX(  jj   *levmult-1, 0)         ][      ii   *levmult               ];
    neigh2d[ic].top    = hash[MIN( (jj+1)*levmult,   jmaxsize-1)][      ii   *levmult               ];
  }

  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  free(levtable);
  genmatrixfree((void**)hash);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}


/**
 * neighbors2d_hashlibcpu_opt_1 calculates neighbors with the library hash and the boundary cell memory optimization
 */
struct neighbor2d *neighbors2d_hashlibcpu_opt_1( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int jmaxsize = mesh_size*levtable[levmx];
  int imaxsize = mesh_size*levtable[levmx]; 
  int numKeys = imaxsize + (jmaxsize - 1) * imaxsize;

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(factory, HASH_PERFECT_HASHES, imaxsize + (jmaxsize - 1) * imaxsize, numKeys, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_SetupTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = numKeys; 
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  int istart;
  int jstart;
  int width;
  for(int ic=0; ic < ncells; ic++){
    if (level[ic] == levmx) {
      intintHash_InsertSingle(hashTable, i[ic] + j[ic] * imaxsize, ic);
    } else {
      width = levtable[levmx-level[ic]] - 1;
      istart = i[ic] * levtable[levmx-level[ic]];
      jstart = j[ic] * levtable[levmx-level[ic]];
      for(int x = 0; x < width; x++){
        intintHash_InsertSingle(hashTable, (istart + x    ) + (jstart        ) * imaxsize, ic);
        intintHash_InsertSingle(hashTable, (istart + x + 1) + (jstart + width) * imaxsize, ic);
        intintHash_InsertSingle(hashTable, (istart + width) + (jstart + x    ) * imaxsize, ic);
        intintHash_InsertSingle(hashTable, (istart        ) + (jstart + x + 1) * imaxsize, ic);
      }
    }
  }

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){ int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    intintHash_QuerySingle(hashTable, (MAX(  ii   *levmult-1, 0         )) + (      jj   *levmult               ) * imaxsize, &neigh2d[ic].left);
    intintHash_QuerySingle(hashTable, (MIN( (ii+1)*levmult,   imaxsize-1)) + (      jj   *levmult               ) * imaxsize, &neigh2d[ic].right);
    intintHash_QuerySingle(hashTable, (      ii   *levmult               ) + (MAX(  jj   *levmult-1, 0)         ) * imaxsize, &neigh2d[ic].bottom);
    intintHash_QuerySingle(hashTable, (      ii   *levmult               ) + (MIN( (jj+1)*levmult,   jmaxsize-1)) * imaxsize, &neigh2d[ic].top);
  }

  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}

/**
 * neighbors2d_hashcpu_opt_2 calculates neighbors with a hard-coded perfect hash and the 7-write, 1 read memory optimization
 */
struct neighbor2d *neighbors2d_hashcpu_opt_2( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int jmaxsize = mesh_size*levtable[levmx];
  int imaxsize = mesh_size*levtable[levmx];

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  int **hash = (int **)genmatrix(jmaxsize, imaxsize, sizeof(int));

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
//for (int jj = 0; jj<jmaxsize; jj++){
//  for (int ii = 0; ii<imaxsize; ii++){
//    hash[jj][ii]=-1;
//  }
//}

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  for(int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    int cellnumber = ic;

    jj = jj*levmult;
    ii = ii*levmult; 

    hash[jj][ii] = cellnumber; // lower left corner

    ii += levmult/2;
    hash[jj][ii] = cellnumber; // lower boundary mid-point
    if(levmult > 2) {
      ii += levmult/2 - 1;
      hash[jj][ii] = cellnumber; // lower right corner
      ii -= levmult/2 - 1;
    }
    ii -= levmult/2;
    jj += levmult/2;
    hash[jj][ii] = cellnumber; // left boundary mid-point
    ii += levmult - 1;
    hash[jj][ii] = cellnumber; // right boundary mid-point

    if(levmult > 2) {
      ii -= levmult - 1;
      jj += levmult/2 - 1;
      hash[jj][ii] = cellnumber; // upper left boundary
      ii += levmult/2;
      hash[jj][ii] = cellnumber; // upper boundary mid-point
    }
  } 

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    neigh2d[ic].left   = hash[      jj   *levmult               ][MAX(  ii   *levmult-1, 0         )];
    neigh2d[ic].right  = hash[      jj   *levmult               ][MIN( (ii+1)*levmult,   imaxsize-1)];
    neigh2d[ic].bottom = hash[MAX(  jj   *levmult-1, 0)         ][      ii   *levmult               ];
    neigh2d[ic].top    = hash[MIN( (jj+1)*levmult,   jmaxsize-1)][      ii   *levmult               ];
  }
  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  free(levtable);
  genmatrixfree((void**)hash);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}

/**
 * neighbors2d_hashlibcpu_opt_2 calculates neighbors with the library hash and the 7-write, 1 read memory optimization
 */
struct neighbor2d *neighbors2d_hashlibcpu_opt_2( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int jmaxsize = mesh_size*levtable[levmx];
  int imaxsize = mesh_size*levtable[levmx];

  int numKeys = imaxsize + (jmaxsize-1) * imaxsize;

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(factory, HASH_PERFECT_HASHES, imaxsize + (jmaxsize - 1) * imaxsize, numKeys, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_SetupTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = numKeys * sizeof(int); 
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  for(int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    int cellnumber = ic;

    jj = jj*levmult;
    ii = ii*levmult; 

    intintHash_InsertSingle(hashTable, jj*imaxsize+ii, cellnumber); // lower left corner

    ii += levmult/2;
    intintHash_InsertSingle(hashTable, jj*imaxsize+ii, cellnumber); // lower boundary mid-point
    if(levmult > 2) {
      ii += levmult/2 - 1;
      intintHash_InsertSingle(hashTable, jj*imaxsize+ii, cellnumber); // lower right corner
      ii -= levmult/2 - 1;
    }
    ii -= levmult/2;
    jj += levmult/2;
    intintHash_InsertSingle(hashTable, jj*imaxsize+ii, cellnumber); // left boundary mid-point
    ii += levmult - 1;
    intintHash_InsertSingle(hashTable, jj*imaxsize+ii, cellnumber); // right boundary mid-point

    if(levmult > 2) {
      ii -= levmult - 1;
      jj += levmult/2 - 1;
      intintHash_InsertSingle(hashTable, jj*imaxsize+ii, cellnumber); // upper left boundary
      ii += levmult/2;
      intintHash_InsertSingle(hashTable, jj*imaxsize+ii, cellnumber); // upper boundary mid-point
    }
  } 

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    intintHash_QuerySingle(hashTable, (MAX(  ii   *levmult-1, 0         )) + (      jj   *levmult               ) * imaxsize, &neigh2d[ic].left);
    intintHash_QuerySingle(hashTable, (MIN( (ii+1)*levmult,   imaxsize-1)) + (      jj   *levmult               ) * imaxsize, &neigh2d[ic].right);
    intintHash_QuerySingle(hashTable, (      ii   *levmult               ) + (MAX(  jj   *levmult-1, 0)         ) * imaxsize, &neigh2d[ic].bottom);
    intintHash_QuerySingle(hashTable, (      ii   *levmult               ) + (MIN( (jj+1)*levmult,   jmaxsize-1)) * imaxsize, &neigh2d[ic].top);
  }

  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}

/**
 * neighbors2d_hashcpu_opt_3 calculates neighbors with a hard-coded perfect hash and the 1-write, 3-read memory optimization
 */
struct neighbor2d *neighbors2d_hashcpu_opt_3(uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int imax = mesh_size;
  int jmax = mesh_size;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int jmaxsize = mesh_size*levtable[levmx];
  int imaxsize = mesh_size*levtable[levmx];

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  int **hash = (int **)genmatrix(jmaxsize, imaxsize, sizeof(int));

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  for (int jj = 0; jj<jmaxsize; jj++){
    for (int ii = 0; ii<imaxsize; ii++){
      hash[jj][ii]=-1;
    }
  }

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  for(int ic=0; ic<ncells; ic++){ 
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    int ii = i[ic]*levmult;
    int jj = j[ic]*levmult;
      
    hash[jj][ii] = ic;
  } 

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);
 
  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic]; 
    int lev = level[ic]; 
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
    if (iicur <    1  ) nlftval = ic;
    if (jjcur <    1  ) nbotval = ic;
    if ((ii + 1)*levmult >(imax)*levtable[levmx]-1) nrhtval = ic;
    if ((jj + 1)*levmult >(jmax)*levtable[levmx]-1) ntopval = ic;
    // Need to check for finer neighbor first
    // Right and top neighbor don't change for finer, so drop through to same size
    // Left and bottom need to be half of same size index for finer test
    if (lev != levmx) {
      int iilftfiner = iicur-(iicur-iilft)/2;
      int jjbotfiner = jjcur-(jjcur-jjbot)/2;
      if (nlftval < 0) nlftval = hash[jjcur][iilftfiner];
      if (nbotval < 0) nbotval = hash[jjbotfiner][iicur];
    }

    // same size neighbor
    if (nlftval < 0) nlftval = hash[jjcur][iilft];
    if (nrhtval < 0) nrhtval = hash[jjcur][iirht];
    if (nbotval < 0) nbotval = hash[jjbot][iicur];
    if (ntopval < 0) ntopval = hash[jjtop][iicur];
    // coarser neighbor
    if (lev != 0){
      if (nlftval < 0) {
        iilft -= iicur-iilft;
        int jjlft = (jj/2)*2*levmult;
        nlftval = hash[jjlft][iilft];
      }
      if (nrhtval < 0) {
        int jjrht = (jj/2)*2*levmult;
        nrhtval = hash[jjrht][iirht];
      }
      if (nbotval < 0) {
        jjbot -= jjcur-jjbot;
        int iibot = (ii/2)*2*levmult;
        nbotval = hash[jjbot][iibot];
      }
      if (ntopval < 0) {
        int iitop = (ii/2)*2*levmult;
        ntopval = hash[jjtop][iitop];
      }
    }
 
    neigh2d[ic].left   = nlftval;
    neigh2d[ic].right  = nrhtval;
    neigh2d[ic].bottom = nbotval;
    neigh2d[ic].top    = ntopval;
  } 

  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  free(levtable);
  genmatrixfree((void**)hash);
  
  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }
 
  return(neigh2d);
} 

/**
 * neighbors2d_hashlibcpu_opt_3 calculates neighbors with the library hash and the 1-write, 3-read memory optimization
 */
struct neighbor2d *neighbors2d_hashlibcpu_opt_3( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int imax = mesh_size;
  int jmax = mesh_size;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  size_t jmaxsize = mesh_size*levtable[levmx];
  size_t imaxsize = mesh_size*levtable[levmx];
  int numKeys = ncells;

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(factory, HASH_TYPE, imaxsize + (jmaxsize - 1) * imaxsize, ncells, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_EmptyTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes;
    if(intintHash_GetTableType(hashTable) & HASH_COMPACT_HASHES){
      nbytes = numKeys * 2.0 / HASH_LOAD_FACTOR * sizeof(int); 
    }else{
      nbytes = numKeys * sizeof(int);
    }
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  for(int ic=0; ic<ncells; ic++){
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    int ii = i[ic]*levmult;
    int jj = j[ic]*levmult;
    intintHash_InsertSingle(hashTable, jj*imaxsize+ii, ic);
  }

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
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
    if (iicur <    1  ) nlftval = ic;
    if (jjcur <    1  ) nbotval = ic;
    if ((ii + 1)*levmult >(imax)*levtable[levmx]-1) nrhtval = ic;
    if ((jj + 1)*levmult >(jmax)*levtable[levmx]-1) ntopval = ic;
    // Need to check for finer neighbor first
    // Right and top neighbor don't change for finer, so drop through to same size
    // Left and bottom need to be half of same size index for finer test
    if (lev != levmx) {
      int iilftfiner = iicur-(iicur-iilft)/2;
      int jjbotfiner = jjcur-(jjcur-jjbot)/2;
      if (nlftval < 0) intintHash_QuerySingle(hashTable, jjcur*imaxsize+iilftfiner, &nlftval);
      if (nbotval < 0) intintHash_QuerySingle(hashTable, jjbotfiner*imaxsize+iicur, &nbotval);
    }

    // same size neighbor
    if (nlftval < 0) intintHash_QuerySingle(hashTable, jjcur*imaxsize+iilft, &nlftval);
    if (nrhtval < 0) intintHash_QuerySingle(hashTable, jjcur*imaxsize+iirht, &nrhtval);
    if (nbotval < 0) intintHash_QuerySingle(hashTable, jjbot*imaxsize+iicur, &nbotval);
    if (ntopval < 0) intintHash_QuerySingle(hashTable, jjtop*imaxsize+iicur, &ntopval);
    // coarser neighbor
    if (lev != 0){
      if (nlftval < 0) {
        iilft -= iicur-iilft;
        int jjlft = (jj/2)*2*levmult;
        intintHash_QuerySingle(hashTable, jjlft*imaxsize+iilft, &nlftval);
      }
      if (nrhtval < 0) {
        int jjrht = (jj/2)*2*levmult;
        intintHash_QuerySingle(hashTable, jjrht*imaxsize+iirht, &nrhtval);
      }
      if (nbotval < 0) {
        jjbot -= jjcur-jjbot;
        int iibot = (ii/2)*2*levmult;
        intintHash_QuerySingle(hashTable, jjbot*imaxsize+iibot, &nbotval);
      }
      if (ntopval < 0) {
        int iitop = (ii/2)*2*levmult;
        intintHash_QuerySingle(hashTable, jjtop*imaxsize+iitop, &ntopval);
      }
    }

    neigh2d[ic].left   = nlftval;
    neigh2d[ic].right  = nrhtval;
    neigh2d[ic].bottom = nbotval;
    neigh2d[ic].top    = ntopval;
  }

  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}

struct neighbor2d *neighbors2d_hasholdlibcpu_opt_3( uint ncells, int mesh_size, int levmx, int *i, int *j, int *level )
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  struct neighbor2d *neigh2d = (struct neighbor2d *)malloc(ncells*sizeof(struct neighbor2d));
  int imax = mesh_size;
  int jmax = mesh_size;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int jmaxsize = mesh_size*levtable[levmx];
  int imaxsize = mesh_size*levtable[levmx];

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  int *hash = compact_hash_init(smallestProthPrimeAbove(ncells / HASH_LOAD_FACTOR), imaxsize, jmaxsize, 0, 1);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    int nbytes = ncells * 2 / HASH_LOAD_FACTOR * sizeof(int); 
    fprintf(fmem, "\t%8d,", nbytes);
  }

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Write Cells
  for(int ic=0; ic<ncells; ic++){
    int lev = level[ic];
    int levmult = levtable[levmx-lev];
    int ii = i[ic]*levmult;
    int jj = j[ic]*levmult;
    write_hash(ic, jj*imaxsize+ii, hash);
  }

  if(DETAILED_TIMING) writeCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Read Neighbors
  for (int ic=0; ic<ncells; ic++){
    int ii = i[ic];
    int jj = j[ic];
    int lev = level[ic];
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
    if (iicur <    1  ) nlftval = ic;
    if (jjcur <    1  ) nbotval = ic;
    if ((ii + 1)*levmult >(imax)*levtable[levmx]-1) nrhtval = ic;
    if ((jj + 1)*levmult >(jmax)*levtable[levmx]-1) ntopval = ic;
    // Need to check for finer neighbor first
    // Right and top neighbor don't change for finer, so drop through to same size
    // Left and bottom need to be half of same size index for finer test
    if (lev != levmx) {
      int iilftfiner = iicur-(iicur-iilft)/2;
      int jjbotfiner = jjcur-(jjcur-jjbot)/2;
      if (nlftval < 0) nlftval = read_hash(jjcur*imaxsize+iilftfiner, hash);
      if (nbotval < 0) nbotval = read_hash(jjbotfiner*imaxsize+iicur, hash);
    }

    // same size neighbor
    if (nlftval < 0) nlftval = read_hash(jjcur*imaxsize+iilft, hash);
    if (nrhtval < 0) nrhtval = read_hash(jjcur*imaxsize+iirht, hash);
    if (nbotval < 0) nbotval = read_hash(jjbot*imaxsize+iicur, hash);
    if (ntopval < 0) ntopval = read_hash(jjtop*imaxsize+iicur, hash);
    // coarser neighbor
    if (lev != 0){
      if (nlftval < 0) {
        iilft -= iicur-iilft;
        int jjlft = (jj/2)*2*levmult;
        nlftval = read_hash(jjlft*imaxsize+iilft, hash);
      }
      if (nrhtval < 0) {
        int jjrht = (jj/2)*2*levmult;
        nrhtval = read_hash(jjrht*imaxsize+iirht, hash);
      }
      if (nbotval < 0) {
        jjbot -= jjcur-jjbot;
        int iibot = (ii/2)*2*levmult;
        nbotval = read_hash(jjbot*imaxsize+iibot, hash);
      }
      if (ntopval < 0) {
        int iitop = (ii/2)*2*levmult;
        ntopval = read_hash(jjtop*imaxsize+iitop, hash);
      }
    }

    neigh2d[ic].left   = nlftval;
    neigh2d[ic].right  = nrhtval;
    neigh2d[ic].bottom = nbotval;
    neigh2d[ic].top    = ntopval;
  }
  if(DETAILED_TIMING) readCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  read_hash_collision_report();
  compact_hash_delete(hash);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeCPUTime, readCPUTime, cleanupCPUTime, totalCPUTime);
  }

  return(neigh2d);
}

#ifdef HAVE_OPENCL
cl_mem neighbors2d_hashgpu( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer=NULL;
  cl_int error = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  free(levtable);
  uint hash_size = (uint)(imaxsize*jmaxsize);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) GPUTime = 0;

  //Create Hash Table
  hash_buffer = hash_init(hash_size, TILE_SIZE, 0, context, queue, &GPUTime);

  if(DETAILED_TIMING) createGPUTime = (double)(GPUTime)*1.0e-9;
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];
  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];

  error = clSetKernelArg(neighbors2d_hashwrite_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashwrite_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  error = clSetKernelArg(neighbors2d_hashread_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  clEnqueueNDRangeKernel(queue, neighbors2d_hashread_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  clReleaseMemObject(hash_buffer);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
}

cl_mem neighbors2d_hashlibgpu( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer=NULL;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  free(levtable);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(CLFactory, HASH_PERFECT_HASHES, imaxsize + (jmaxsize - 1) * imaxsize, imaxsize + (jmaxsize - 1) * imaxsize, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_SetupTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    int nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8d,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];

  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibwrite_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  error = clSetKernelArg(neighbors2d_hashlibread_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibread_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
}

cl_mem neighbors2d_hashgpu_opt_1( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  free(levtable);
  uint hash_size = (uint)(imaxsize*jmaxsize);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) GPUTime = 0;

  //Create Hash Table
  hash_buffer = hash_init(hash_size, TILE_SIZE, 0, context, queue, &GPUTime);

  if(DETAILED_TIMING) createGPUTime = (double)(GPUTime)*1.0e-9;
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
   
  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    int nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8d,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];
  
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_1_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashwrite_opt_1_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  error = clSetKernelArg(neighbors2d_hashread_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashread_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  clReleaseMemObject(hash_buffer);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
} 

cl_mem neighbors2d_hashlibgpu_opt_1( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer=NULL;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  int numKeys = imaxsize + (jmaxsize - 1) * imaxsize;
  free(levtable);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(CLFactory, HASH_PERFECT_HASHES, imaxsize + (jmaxsize - 1) * imaxsize, numKeys, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_SetupTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    int nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8d,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];

  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_1_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibwrite_opt_1_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  error = clSetKernelArg(neighbors2d_hashlibread_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibread_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
}

cl_mem neighbors2d_hashgpu_opt_2( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  free(levtable);
  uint hash_size = (uint)(imaxsize*jmaxsize);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) GPUTime = 0;

  //Create Hash Table
  hash_buffer = hash_init(hash_size, TILE_SIZE, 0, context, queue, &GPUTime);

  if(DETAILED_TIMING) createGPUTime = (double)(GPUTime)*1.0e-9;
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  
  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    int nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8d,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];
  
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_2_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashwrite_opt_2_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  error = clSetKernelArg(neighbors2d_hashread_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashread_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  clReleaseMemObject(hash_buffer);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);

}

cl_mem neighbors2d_hashlibgpu_opt_2( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer=NULL;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  int numKeys = imaxsize + (jmaxsize - 1) * imaxsize;
  free(levtable);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(CLFactory, HASH_PERFECT_HASHES, imaxsize + (jmaxsize - 1) * imaxsize, ncells, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_SetupTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    int nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8d,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];

  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_2_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibwrite_opt_2_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  error = clSetKernelArg(neighbors2d_hashlibread_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibread_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
}

cl_mem neighbors2d_hashgpu_opt_3( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  free(levtable);
  uint hash_size = (uint)(imaxsize*jmaxsize);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) GPUTime = 0;

  //Create Hash Table
  hash_buffer = hash_init(hash_size, TILE_SIZE, 1, context, queue, &GPUTime);

  if(DETAILED_TIMING) createGPUTime = (double)(GPUTime)*1.0e-9;
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  
  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    int nbytes = jmaxsize * imaxsize * sizeof(int);
    fprintf(fmem, "\t%8d,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];
  
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashwrite_opt_3_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashwrite_opt_3_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 7, sizeof(cl_mem), &hash_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashread_opt_3_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashread_opt_3_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  clReleaseMemObject(hash_buffer);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
}

cl_mem neighbors2d_hashlibgpu_opt_3( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer=NULL;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  int numKeys = ncells;
  free(levtable);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Create Hash Table
  intintHash_Table *hashTable = intintHash_CreateTable(CLFactory, HASH_TYPE, imaxsize + (jmaxsize - 1) * imaxsize, numKeys, HASH_LOAD_FACTOR);

  if(DETAILED_TIMING) createCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  intintHash_SetupTable(hashTable);

  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    size_t nbytes;
    if(intintHash_GetTableType(hashTable) & HASH_COMPACT_HASHES){
      nbytes = numKeys * 2.0 / HASH_LOAD_FACTOR * sizeof(int); 
    }else{
      nbytes = numKeys * sizeof(int);
    }
    fprintf(fmem, "\t%8ld,", nbytes);
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];

  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibwrite_opt_3_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashwrite_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibwrite_opt_3_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 0, sizeof(cl_uint), &ncells);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 1, sizeof(cl_int), &mesh_size);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 2, sizeof(cl_int), &levmx);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 3, sizeof(cl_mem), &levtable_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 4, sizeof(cl_mem), &i_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 5, sizeof(cl_mem), &j_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 6, sizeof(cl_mem), &level_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 7, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(hashTable));
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  error = clSetKernelArg(neighbors2d_hashlibread_opt_3_kern, 8, sizeof(cl_mem), &neighbor2d_buffer);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  cl_event hashread_event;

  error = clEnqueueNDRangeKernel(queue, neighbors2d_hashlibread_opt_3_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  intintHash_DestroyTable(hashTable);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createCPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
}

cl_mem neighbors2d_hasholdlibgpu_opt_3( uint ncells, int mesh_size, int levmx, cl_mem i_buffer, cl_mem j_buffer,
      cl_mem level_buffer, cl_mem levtable_buffer)
{
  struct timespec tSection, tTotal;
  double setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime;
  long startGPUTime, endGPUTime, GPUTime;
  if(DETAILED_TIMING) cpu_timer_start(&tTotal);
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Setup
  cl_mem hash_buffer, neighbor2d_buffer;
  cl_int error = 0;
  long gpu_time = 0;
  int *levtable = (int *)malloc((levmx+1)*sizeof(int));
  for (int lev=0; lev<levmx+1; lev++){
    levtable[lev] = (int)pow(2,lev);
  }
  int imaxsize = mesh_size*levtable[levmx];
  int jmaxsize = mesh_size*levtable[levmx];
  free(levtable);

  if(DETAILED_TIMING) setupCPUTime = cpu_timer_stop(tSection);
  if(DETAILED_TIMING) GPUTime = 0;

  //Create Hash Table
  size_t hash_size;
  double hash_mult = 3.0;
  double mem_opt_factor = 1.0;
  int gpu_do_compact_hash = 0; 
  int gpu_hash_method = 1;
  ulong gpu_hash_table_size;
  ulong gpu_AA=0;
  ulong gpu_BB=0;
  ulong prime=4294967291;
  int hash_report_level = 1;

  uint gpu_compact_hash_size = (uint)((double)ncells*hash_mult);
  uint gpu_perfect_hash_size = (uint)(imaxsize*jmaxsize);
  float gpu_hash_mem_factor = 20.0;
  float gpu_hash_mem_ratio = (double)gpu_perfect_hash_size/(double)gpu_compact_hash_size;
  if (mem_opt_factor != 1.0) gpu_hash_mem_factor /= (mem_opt_factor*0.2);
  gpu_do_compact_hash = (gpu_hash_mem_ratio < gpu_hash_mem_factor) ? 0 : 1; 
  #ifdef __APPLE_CC__
    gpu_do_compact_hash = 0; 
  #endif
  gpu_do_compact_hash = 1; 
   
  if (gpu_do_compact_hash) {
    gpu_hash_method = 1; 
    gpu_hash_table_size = gpu_compact_hash_size;
    hash_size = 2*gpu_compact_hash_size;
    gpu_AA = (ulong)(1.0+(double)(prime-1)*drand48());
    gpu_BB = (ulong)(0.0+(double)(prime-1)*drand48());
    if (gpu_AA > prime-1 || gpu_BB > prime-1) exit(0);
    if (hash_report_level > 1) printf("Factors AA %lu BB %lu\n",gpu_AA,gpu_BB);
  } else {
    gpu_hash_method = -1; 
    gpu_hash_table_size = gpu_perfect_hash_size;
    hash_size = gpu_perfect_hash_size;
  }
  
  hash_buffer = hash_init(hash_size, TILE_SIZE, 1, context, queue, &GPUTime);

  if(DETAILED_TIMING) createGPUTime = (double)(GPUTime)*1.0e-9;
  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Empty Hash Table
  
  if(DETAILED_TIMING) emptyCPUTime = cpu_timer_stop(tSection);

  if(WRITE_MEM_USAGE){
    fprintf(fmem, "\t%8ld,", hash_size*sizeof(int));
  }

  //Write Cells
  size_t global_work_size[1];
  size_t local_work_size[1];

  local_work_size[0] = TILE_SIZE;
  global_work_size[0] = ((ncells+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];
  
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 0,  sizeof(cl_uint),  &ncells), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 1,  sizeof(cl_int),   &levmx), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 2,  sizeof(cl_int),   &imaxsize), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 3,  sizeof(cl_int),   &gpu_hash_method), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 4,  sizeof(cl_ulong), &gpu_hash_table_size), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 5,  sizeof(cl_ulong), &gpu_AA), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 6,  sizeof(cl_ulong), &gpu_BB), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 7,  sizeof(cl_mem),   &levtable_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 8,  sizeof(cl_mem),   &level_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 9,  sizeof(cl_mem),   &i_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 10, sizeof(cl_mem),   &j_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibwrite_opt_3_kern, 11, sizeof(cl_mem),   &hash_buffer), __FILE__, __LINE__);

  cl_event hashwrite_event;

  cl_error_check(clEnqueueNDRangeKernel(queue, neighbors2d_hasholdlibwrite_opt_3_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashwrite_event), __FILE__, __LINE__);
  
  //Read Neighbors
  neighbor2d_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, ncells*sizeof(cl_uint4), NULL, &error);
  if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 0, sizeof(cl_uint), &ncells), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 1, sizeof(cl_int), &mesh_size), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 2, sizeof(cl_int), &levmx), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 3, sizeof(cl_mem), &levtable_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 4, sizeof(cl_mem), &i_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 5, sizeof(cl_mem), &j_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 6, sizeof(cl_mem), &level_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 7, sizeof(cl_mem), &hash_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 8, sizeof(cl_mem), &neighbor2d_buffer), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 9, sizeof(cl_int), &gpu_hash_method), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 10, sizeof(cl_ulong), &gpu_hash_table_size), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 11, sizeof(cl_ulong), &gpu_AA), __FILE__, __LINE__);
  cl_error_check(clSetKernelArg(neighbors2d_hasholdlibread_opt_3_kern, 12, sizeof(cl_ulong), &gpu_BB), __FILE__, __LINE__);
  
  cl_event hashread_event;

  cl_error_check(clEnqueueNDRangeKernel(queue, neighbors2d_hasholdlibread_opt_3_kern, 1, 0, global_work_size, local_work_size, 0, NULL, &hashread_event), __FILE__, __LINE__);

  clWaitForEvents(1,&hashread_event);

  if(DETAILED_TIMING) cpu_timer_start(&tSection);

  //Cleanup
  clReleaseMemObject(hash_buffer);

  if (DETAILED_TIMING){
    cleanupCPUTime = cpu_timer_stop(tSection);
    totalCPUTime = cpu_timer_stop(tTotal);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashwrite_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    writeGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_START, sizeof(startGPUTime), &startGPUTime, NULL);
    clGetEventProfilingInfo(hashread_event, CL_PROFILING_COMMAND_END, sizeof(endGPUTime), &endGPUTime, NULL);
    readGPUTime = (double)(endGPUTime - startGPUTime)*1.0e-9;
    printf("\nDetailed Timing:\nSetup: %lf Create: %lf Empty: %lf Write: %lf Read: %lf Cleanup: %lf Total: %lf\n",setupCPUTime, createGPUTime, emptyCPUTime, writeGPUTime, readGPUTime, cleanupCPUTime, totalCPUTime);
  }

  clReleaseEvent(hashwrite_event);
  clReleaseEvent(hashread_event);

  return(neighbor2d_buffer);
} 
#endif
 
// adaptiveMeshConstructor()
// Inputs: n (width/height of the square mesh), levmax (maximum level of refinement),
//         pointers for the level, x, and y arrays (should be NULL for all three)
// Output: number of cells in the adaptive mesh
//
int adaptiveMeshConstructorWij(const int n, const int levmax, 
         int** level_ptr, double** x_ptr, double** y_ptr, int **i_ptr, int **j_ptr, float threshold, int target_ncells) {
  int ncells = SQ(n);

  // ints used for for() loops later
  int ic, xc, yc, xlc, ylc, nlc;

  //printf("\nBuilding the mesh...\n");

  // Initialize Coarse Mesh
  int*  level = (int*)  malloc(sizeof(int)*ncells);
  double* x   = (double*) malloc(sizeof(double)*ncells);
  double* y   = (double*) malloc(sizeof(double)*ncells);
  int*  i     = (int*)  malloc(sizeof(int)*ncells);
  int*  j     = (int*)  malloc(sizeof(int)*ncells);
  for(yc = 0; yc < n; yc++) {
    for(xc = 0; xc < n; xc++) {
      level[n*yc+xc] = 0;
      x[n*yc+xc]     = (real)(TWO*xc+ONE) / (real)(TWO*n);
      y[n*yc+xc]     = (real)(TWO*yc+ONE) / (real)(TWO*n);
      i[n*yc+xc]     = xc;
      j[n*yc+xc]     = yc;
    }
  }
  //printf("Coarse mesh initialized.\n");

  // Randomly Set Level of Refinement
  //unsigned int iseed = (unsigned int)time(NULL);
  //srand (iseed);
  //srand (0);
  for(int ii = levmax; ii > 0; ii--) {
    float lev_threshold = threshold*(float)ii/(float)levmax;
    for(ic = 0; ic < ncells; ic++) {
      float jj = (100.0*(float)rand() / ((float)RAND_MAX));
      if(jj<lev_threshold && level[ic] == 0) level[ic] = ii;
    }
  }

  //printf("Levels of refinement randomly set.\n");

  // Smooth the Refinement
  int newcount = -1;
  while(newcount != 0) {
    newcount = 0;
    int lev = 0;
    for(ic = 0; ic < ncells; ic++) {
      lev = level[ic];
      lev++;
      // Check bottom neighbor
      if(ic - n >= 0) {
        if(level[ic-n] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
      // Check top neighbor
      if(ic + n < ncells) {
        if(level[ic+n] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
      // Check left neighbor
      if((ic%n)-1 >= 0) {
        if(level[ic-1] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
      // Check right neighbor
      if((ic%n)+1 < n) {
        if(level[ic+1] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
    }
  }

  //printf("\nDEBUG -- ncells %d target_ncells %ld fine mesh size %ld\n",ncells,target_ncells,n*two_to_the(levmax)*n*two_to_the(levmax));
  if (target_ncells > ncells && target_ncells < n*2<<(levmax)*n*2<<(levmax)) {
    int icount = 0;
    int newcount = 0;
    for(ic = 0; ic < ncells; ic++) {newcount += (powerOfFour(level[ic]) - 1);}

    while ( abs((ncells+newcount) - target_ncells) > MAX(5,target_ncells/10000) && icount < 40) {
      icount++;
      //printf("DEBUG -- Adjusting cell count %ld target %ld diff %ld\n",ncells+newcount, target_ncells, abs(ncells+newcount - target_ncells));

      if (ncells+newcount > target_ncells){
        int reduce_count = ((ncells+newcount) - target_ncells);
        //printf("DEBUG -- Too many cells -- need to reduce by %d\n",reduce_count);
        int jcount = 0;
        while (reduce_count > 0 && jcount < ncells) {
          int jj = 1 + (int)((float)ncells*rand() / (RAND_MAX+1.0));
          if(jj>0 && jj<ncells && level[jj] > 0) {
             reduce_count-=4;
          //   printf("DEBUG reducing level for ic %d level %d reduce_count %d jj %d\n",jj,level[jj],reduce_count,jj);
             level[jj]--;
          }
          jcount++;
        }
      } else {
        int increase_count = (target_ncells - (ncells+newcount));
        increase_count /= (levmax*4);
        //printf("DEBUG -- Too few cells -- need to increase by %d\n",increase_count);
        int jcount = 0;
        while (increase_count > 0 && jcount < ncells) {
          int jj = 1 + (int)((float)ncells*rand() / (RAND_MAX+1.0));
          if(jj>0 && jj<ncells && level[jj] < levmax) {
            increase_count-=4;
            level[jj]++;
          }
          jcount++;
        }
      }

      // Smooth the Refinement
      newcount = -1;
      while(newcount != 0) {
        newcount = 0;
        uint lev = 0;
        for(ic = 0; ic < ncells; ic++) {
          lev = level[ic];
          lev++;
          // Check bottom neighbor
          if(ic - n >= 0) {
            if(level[ic-n] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
          // Check top neighbor
          if(ic + n < ncells) {
            if(level[ic+n] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
          // Check left neighbor
          if((ic%n)-1 >= 0) {
            if(level[ic-1] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
          // Check right neighbor
          if((ic%n)+1 < n) {
            if(level[ic+1] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
        }
      } // while(newcount != 0) {
      newcount = 0;
      for(ic = 0; ic < ncells; ic++) {newcount += (powerOfFour(level[ic]) - 1);}
    } // while ( abs(ncells+newcount - target_ncells) > 10 && icount < 10) {

  } //if (target_ncells > 0) {

  //printf("Refinement smoothed.\n");
  int small_cells = 0;
  for(ic = 0; ic < ncells; ic++) {
    if (level[ic] == levmax) {
      small_cells++;
    }
  }
  //printf("%8d small cells, ", small_cells);
  
  // Allocate Space for the Adaptive Mesh
  newcount = 0;
  for(ic = 0; ic < ncells; ic++) {newcount += (powerOfFour(level[ic]) - 1);}
  int*  level_temp = (int*)  malloc(sizeof(int)*(ncells+newcount));
  double* x_temp   = (double*) malloc(sizeof(double)*(ncells+newcount));
  double* y_temp   = (double*) malloc(sizeof(double)*(ncells+newcount));
  int*  i_temp     = (int*)  malloc(sizeof(int)*(ncells+newcount));
  int*  j_temp     = (int*)  malloc(sizeof(int)*(ncells+newcount));

  // Set the Adaptive Mesh
  int offset = 0;
  for(yc = 0; yc < n; yc++) {
    for(xc = 0; xc < n; xc++) {
      ic = n*yc + xc;
      nlc = (int) sqrt( (real) powerOfFour(level[ic]) );
      for(ylc = 0; ylc < nlc; ylc++) {
        for(xlc = 0; xlc < nlc; xlc++) {
          level_temp[ic + offset + (nlc*ylc + xlc)] = level[ic];
          x_temp[ic + offset + (nlc*ylc + xlc)] = x[ic]-(ONE / (real)(TWO*n))
                               + ((real)(TWO*xlc+ONE) / (real)(n*nlc*TWO));
          y_temp[ic + offset + (nlc*ylc + xlc)] = y[ic]-(ONE / (real)(TWO*n))
                               + ((real)(TWO*ylc+ONE) / (real)(n*nlc*TWO));
          i_temp[ic + offset + (nlc*ylc + xlc)] = i[ic]*pow(2,level[ic]) + xlc;
          j_temp[ic + offset + (nlc*ylc + xlc)] = j[ic]*pow(2,level[ic]) + ylc;
        }         
      }
      offset += powerOfFour(level[ic])-1;
    }
  }
  //printf("Adaptive mesh built.\n");

  // Swap pointers and free memory used by Coarse Mesh
  swap_int(&level, &level_temp);
  swap_double(&x, &x_temp);
  swap_double(&y, &y_temp);
  swap_int(&i, &i_temp);
  swap_int(&j, &j_temp);
  free(level_temp);
  free(x_temp);
  free(y_temp);
  free(i_temp);
  free(j_temp);

  //printf("Old ncells: %d", ncells);
  // Update ncells
  ncells += newcount;
  //printf("\tNew ncells: %d\n", ncells);

  if (randomize) {
    // Randomize the order of the arrays

    int* random = (int*) malloc(sizeof(int)*ncells);
    int* temp1 = (int*) malloc(sizeof(int)*ncells);
    real* temp2 = (real*) malloc(sizeof(real)*ncells*2);
    int* temp3 = (int*) malloc(sizeof(int)*ncells*2);
    // XXX Want better randomization? XXX
    // XXX Why is the time between printf() statements the longest part? XXX
    //printf("Shuffling");
    //fflush(stdout);
    for(ic = 0; ic < ncells; ic++) {random[ic] = ic;}
    //iseed = (unsigned int)time(NULL);
    //srand (iseed);
    srand(0);
    nlc = 0;
    for(int ii = 0; ii < 7; ii++) {
      for(ic = 0; ic < ncells; ic++) {
    
        int jj = (int)( (real)ncells*((real)rand() / (real)(RAND_MAX+ONE) ) );
        // occasionally jj will be ncells and random ratio is 1.0
        if (jj >= ncells) jj=ncells-1;
        nlc = random[jj];
        random[jj] = random[ic];
        random[ic] = nlc;
         if (random[ic] >= ncells) {
           printf("DEBUG -- ic %d file %s line %d\n",ic,__FILE__,__LINE__);
           exit(0);
         }
      }
      //printf(".");
      //fflush(stdout);
    }
    //printf("\n");

    for(ic = 0; ic < ncells; ic++) {
      temp1[ic] = level[random[ic]];
      temp2[2*ic] = x[random[ic]];
      temp2[2*ic+1] = y[random[ic]];
      temp3[2*ic] = i[random[ic]];
      temp3[2*ic+1] = j[random[ic]];
    }
    for(ic = 0; ic < ncells; ic++) {
      level[ic] = temp1[ic];
      x[ic]     = temp2[2*ic];
      y[ic]     = temp2[2*ic+1];
      i[ic]     = temp3[2*ic];
      j[ic]     = temp3[2*ic+1];
    }

    free(temp1);
    free(temp2);
    free(temp3);
    free(random);
    //printf("Adaptive mesh randomized.\n");
  } // End of if randomize

  *level_ptr = level;
  *x_ptr = x;
  *y_ptr = y;
  *i_ptr = i;
  *j_ptr = j;

  //printf("Adaptive mesh construction complete.\n");
  return ncells;
}

void *genvector(int inum, size_t elsize)
{
  void *out;
  size_t mem_size;

  mem_size = inum*elsize;
  out      = (void *)calloc((size_t)inum, elsize);

  return (out);
}

void genvectorfree(void *var)
{
  free(var);
}

void **genmatrix(int jnum, int inum, size_t elsize)
{
  void **out;
  int mem_size;

  mem_size = jnum*sizeof(void *);
  out      = (void **)malloc(mem_size);

  mem_size = jnum*inum*elsize;
  out[0]    = (void *)calloc(jnum*inum, elsize);

  for (int i = 1; i < jnum; i++) {
    out[i] = (void *)((char *)out[i-1] + inum*elsize);
  }

  return (out);
}

void genmatrixfree(void **var)
{
  free(var[0]);
  free(var);
}
