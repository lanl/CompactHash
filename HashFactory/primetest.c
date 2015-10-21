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

#include <math.h>
#include <stdio.h>
#define LCG_A 2147483629
#define LCG_C 2147483587
#define LCG_M 2147483647
#define PRIME_NUM_CHECKS 20 

int modularPow(int base, int exponent, int modulus){
    int result = 1;
    while (exponent){
        if (exponent & 1)
            result = ((long long int)result * base) % modulus;
        exponent >>= 1;
        base = ((long long int)base * base) % modulus;
    }
    return result;
}

int largestProthPrimeUnder(int N){
  if(N < 4){
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
  while(p > 3){
    //check if a proth number is prime
    for(i = 0; i < PRIME_NUM_CHECKS; i++){
      a = rand();
      if(modularPow(a, (p - 1) / 2, p) == p - 1){
        return p;
      }
    }
    //determine the next proth number
    if(p - 1 == s * s / 4){
      s /= 2;
    }
    p -= s;
  }
  return 3;
}

int smallestProthPrimeAbove(int N){
  if(N < 4){
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
  while(1){
    //determine the next proth number
    if(p - 1 == s * s){
      s *= 2;
    }
    p += s;
    //check if a proth number is prime
    for(i = 0; i < PRIME_NUM_CHECKS; i++){
      a = rand();
      if(modularPow(a, (p - 1) / 2, p) == p - 1){
        return p;
      }
    }
  }
  return 3;
}

int main(){
  printf("%d\n", largestProthPrimeUnder(2));
  printf("%d\n", smallestProthPrimeAbove(2));
  printf("%d\n", largestProthPrimeUnder(3));
  printf("%d\n", smallestProthPrimeAbove(3));
  printf("%d\n", largestProthPrimeUnder(4));
  printf("%d\n", smallestProthPrimeAbove(4));
  printf("%d\n", largestProthPrimeUnder(5));
  printf("%d\n", smallestProthPrimeAbove(5));
  printf("%d\n", largestProthPrimeUnder(11));
  printf("%d\n", smallestProthPrimeAbove(11));
  printf("%d\n", largestProthPrimeUnder(13));
  printf("%d\n", smallestProthPrimeAbove(13));
  printf("%d\n", largestProthPrimeUnder(2000024));
  printf("%d\n", smallestProthPrimeAbove(1993728));
  int p = 3;
  int m = 2;
  int c = 1;
  int i;
  for(i = 0; i < 0; i++){
     printf("%d %d\n", i, p);
     if(c == 0){
       m *= 2;
       c = i + 2;
     }
     p += m;
     c--;
  }
  return 0;
}
